#!/usr/bin/env Rscript
###############################################################################
# run_1gb_c1c2c3.R -- 1GB benchmark for C1, C2, C3 with time_df=24
#
# Runs all four pipeline stages: crossbasis, glm, crosspred, crossreduce.
# Merges results into timing_results.csv, substage_timings.csv,
# and scaling_analysis.csv.
###############################################################################

library(data.table)
library(splines)
pkgload::load_all(".", quiet = TRUE)

TIME_DF <- 24L
SCALE   <- "1gb"

configs <- list(
  C1 = list(
    argvar  = list(fun = "lin"),
    arglag  = list(fun = "poly", degree = 4),
    lag     = c(0L, 15L),
    cb_cols = 5L,
    desc    = "Minimal/DLM: lin x poly(4), lag 0-15"
  ),
  C2 = list(
    argvar  = list(fun = "ns", df = 5),
    arglag  = list(fun = "ns", df = 4),
    lag     = c(0L, 21L),
    cb_cols = 20L,
    desc    = "Typical epi: ns(5) x ns(4), lag 0-21"
  ),
  C3 = list(
    argvar  = list(fun = "bs", df = 6),
    arglag  = list(fun = "ns", df = 4),
    lag     = c(0L, 40L),
    cb_cols = 24L,
    desc    = "Extended lag: bs(6) x ns(4), lag 0-40"
  )
)

# -- Load data ----------------------------------------------------------------
cat("Loading scale_1gb.rds...\n")
t0 <- proc.time()
dt <- readRDS("benchmarks/data/scale_1gb.rds")
cat(sprintf("%.1fs  (%s rows)\n\n", (proc.time() - t0)[3],
            format(nrow(dt), big.mark = ",")))

n_rows   <- nrow(dt)
x        <- dt$temp
death    <- dt$death
time_var <- dt$time
dow      <- dt$dow
group    <- dt$city
rm(dt); gc(verbose = FALSE)

df_model <- data.frame(death = death, time = time_var, dow = dow)
rm(death, time_var, dow); gc(verbose = FALSE)

# -- Results accumulators -----------------------------------------------------
new_timing   <- list()
new_substage <- list()
t_idx <- 0L
s_idx <- 0L

# -- Loop over configs --------------------------------------------------------
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  sep <- paste(rep("-", 60), collapse = "")
  cat(sprintf("%s\n%s: %s\n%s\n", sep, cfg_name, cfg$desc, sep))

  # ---- crossbasis -----------------------------------------------------------
  cat("  crossbasis... ")
  gc(verbose = FALSE)
  t0 <- proc.time()
  cb <- crossbasis(x, lag = cfg$lag, argvar = cfg$argvar, arglag = cfg$arglag,
                   group = group)
  cb_sec <- (proc.time() - t0)[3]
  cat(sprintf("%.3fs  (%d cols, %.0f MB)\n", cb_sec, ncol(cb),
              object.size(cb) / 1e6))

  t_idx <- t_idx + 1L
  new_timing[[t_idx]] <- data.frame(
    scale = SCALE, config = cfg_name, stage = "crossbasis",
    n_rows = n_rows, cb_cols = ncol(cb), median_time_sec = cb_sec,
    is_extrapolated = FALSE, stringsAsFactors = FALSE
  )

  # ---- crossbasis sub-stages ------------------------------------------------
  cat("  profiling crossbasis internals...\n")
  basisvar <- do.call(onebasis, modifyList(cfg$argvar, list(x = as.numeric(x))))
  al <- cfg$arglag
  if ((is.null(al$fun) || "intercept" %in% names(formals(al$fun))) &&
      sum(pmatch(names(al), "intercept", nomatch = 0)) == 0)
    al$intercept <- TRUE
  al$cen <- NULL
  basislag <- do.call(onebasis, modifyList(al, list(x = seqlag(cfg$lag))))
  nv <- ncol(basisvar); nl <- ncol(basislag)

  gc(verbose = FALSE); t0 <- proc.time()
  lag_mats <- vector("list", nv)
  for (v in seq_len(nv))
    lag_mats[[v]] <- as.matrix(tsModel::Lag(basisvar[, v], seqlag(cfg$lag),
                                            group = group))
  lag_mat_sec <- (proc.time() - t0)[3]

  gc(verbose = FALSE); t0 <- proc.time()
  cb_tmp <- matrix(0, nrow = n_rows, ncol = nv * nl)
  for (v in seq_len(nv))
    for (l in seq_len(nl))
      cb_tmp[, nl * (v - 1) + l] <- lag_mats[[v]] %*% basislag[, l]
  multiply_sec <- (proc.time() - t0)[3]
  rm(cb_tmp, lag_mats, basisvar, basislag); gc(verbose = FALSE)

  cat(sprintf("    lag_matrix:    %.3fs\n", lag_mat_sec))
  cat(sprintf("    multiply_loop: %.3fs\n", multiply_sec))

  for (sub in list(list("lag_matrix", lag_mat_sec),
                   list("multiply_loop", multiply_sec))) {
    s_idx <- s_idx + 1L
    new_substage[[s_idx]] <- data.frame(
      scale = SCALE, config = cfg_name, substage = sub[[1]],
      time_seconds = sub[[2]], stringsAsFactors = FALSE
    )
  }

  # ---- glm ------------------------------------------------------------------
  total_cols <- ncol(cb) + TIME_DF + 6L
  cat(sprintf("  glm (time_df=%d, ~%.1f GB model matrix)... ",
              TIME_DF, n_rows * total_cols * 8 / 1e9))
  gc(verbose = FALSE)
  t0 <- proc.time()
  model <- tryCatch(
    glm(death ~ cb + ns(time, df = TIME_DF) + dow,
        data = df_model, family = quasipoisson()),
    error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
  )
  glm_sec <- (proc.time() - t0)[3]

  if (!is.null(model)) {
    cat(sprintf("%.3fs  (converged=%s, iter=%d)\n",
                glm_sec, model$converged, model$iter))
  } else {
    glm_sec <- NA_real_
  }

  t_idx <- t_idx + 1L
  new_timing[[t_idx]] <- data.frame(
    scale = SCALE, config = cfg_name, stage = "glm",
    n_rows = n_rows, cb_cols = ncol(cb), median_time_sec = glm_sec,
    is_extrapolated = FALSE, stringsAsFactors = FALSE
  )

  # ---- crosspred + crossreduce (only if model succeeded) --------------------
  if (!is.null(model)) {
    cat("  crosspred... ")
    gc(verbose = FALSE)
    t0 <- proc.time()
    pred <- tryCatch(
      crosspred(cb, model,
                at  = seq(attr(cb, "range")[1], attr(cb, "range")[2],
                          length.out = 50),
                cumul = TRUE),
      error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
    )
    pred_sec <- (proc.time() - t0)[3]
    if (!is.null(pred)) cat(sprintf("%.3fs\n", pred_sec)) else pred_sec <- NA_real_

    t_idx <- t_idx + 1L
    new_timing[[t_idx]] <- data.frame(
      scale = SCALE, config = cfg_name, stage = "crosspred",
      n_rows = n_rows, cb_cols = ncol(cb), median_time_sec = pred_sec,
      is_extrapolated = FALSE, stringsAsFactors = FALSE
    )

    cat("  crossreduce... ")
    gc(verbose = FALSE)
    t0 <- proc.time()
    cr <- tryCatch(
      crossreduce(cb, model, type = "overall"),
      error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
    )
    cr_sec <- (proc.time() - t0)[3]
    if (!is.null(cr)) cat(sprintf("%.3fs\n", cr_sec)) else cr_sec <- NA_real_

    t_idx <- t_idx + 1L
    new_timing[[t_idx]] <- data.frame(
      scale = SCALE, config = cfg_name, stage = "crossreduce",
      n_rows = n_rows, cb_cols = ncol(cb), median_time_sec = cr_sec,
      is_extrapolated = FALSE, stringsAsFactors = FALSE
    )

    rm(pred, cr, model)
  }

  rm(cb); gc(verbose = FALSE)
  cat("\n")
}

# -- Merge with existing results ----------------------------------------------
results_dir   <- "benchmarks/results"
timing_file   <- file.path(results_dir, "timing_results.csv")
substage_file <- file.path(results_dir, "substage_timings.csv")

existing_timing   <- fread(timing_file)
existing_substage <- fread(substage_file)

new_timing_dt   <- rbindlist(new_timing)
new_substage_dt <- rbindlist(new_substage)

# Remove old 1gb rows for C1/C2/C3 (both actuals and extrapolated)
updated_timing <- rbind(
  existing_timing[!(scale == SCALE & config %in% names(configs))],
  new_timing_dt,
  fill = TRUE
)
updated_timing[is.na(is_extrapolated), is_extrapolated := FALSE]

updated_substage <- rbind(
  existing_substage[!(scale == SCALE & config %in% names(configs))],
  new_substage_dt
)

fwrite(updated_timing,   timing_file)
fwrite(updated_substage, substage_file)
fwrite(updated_timing[!is.na(median_time_sec)],
       file.path(results_dir, "scaling_analysis.csv"))

# -- Print final table --------------------------------------------------------
cat("=== Updated 1GB rows ===\n")
print(updated_timing[scale == SCALE,
                      .(scale, config, stage, n_rows, cb_cols,
                        median_time_sec, is_extrapolated)])

cat("\n=== Scaling exponents (crossbasis + glm) ===\n")
for (cfg_name in names(configs)) {
  for (stg in c("crossbasis", "glm")) {
    sub <- updated_timing[config == cfg_name & stage == stg &
                          !is_extrapolated & !is.na(median_time_sec)]
    sub <- sub[order(n_rows)]
    if (nrow(sub) >= 2) {
      fit   <- lm(log10(median_time_sec) ~ log10(n_rows), data = sub)
      alpha <- coef(fit)[2]
      cat(sprintf("  %s / %-12s  alpha=%.2f\n", cfg_name, stg, alpha))
    }
  }
}

cat("\nDone.\n")
