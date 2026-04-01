#!/usr/bin/env Rscript
###############################################################################
# run_glm_backends.R -- Compare GLM backends at 1GB scale
#
# Tests glm (base), bam, bam+discrete (mgcv), and fastglm for the GLM stage
# on C1, C2, C3 configs with time_df=24.
#
# crosspred/crossreduce are intentionally skipped for non-base backends.
# Results saved to benchmarks/results/glm_backends.csv.
###############################################################################

library(data.table)
library(splines)
library(mgcv)
pkgload::load_all(".", quiet = TRUE)

has_fastglm <- requireNamespace("fastglm", quietly = TRUE)
if (has_fastglm) {
  library(fastglm)
  cat("fastglm: available\n")
} else {
  cat("fastglm: NOT installed -- skipping\n")
}

n_cores  <- max(1L, parallel::detectCores() - 1L)
TIME_DF  <- 24L
SCALE    <- "1gb"

configs <- list(
  C1 = list(
    argvar  = list(fun = "lin"),
    arglag  = list(fun = "poly", degree = 4),
    lag     = c(0L, 15L),
    desc    = "Minimal/DLM: lin x poly(4), lag 0-15"
  ),
  C2 = list(
    argvar  = list(fun = "ns", df = 5),
    arglag  = list(fun = "ns", df = 4),
    lag     = c(0L, 21L),
    desc    = "Typical epi: ns(5) x ns(4), lag 0-21"
  ),
  C3 = list(
    argvar  = list(fun = "bs", df = 6),
    arglag  = list(fun = "ns", df = 4),
    lag     = c(0L, 40L),
    desc    = "Extended lag: bs(6) x ns(4), lag 0-40"
  )
)

# -- Load data ----------------------------------------------------------------
cat("\nLoading scale_1gb.rds...\n")
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

# -- Results accumulator ------------------------------------------------------
results <- list()
r_idx   <- 0L

run_backend <- function(label, expr_fn) {
  gc(verbose = FALSE)
  t0  <- proc.time()
  out <- tryCatch(expr_fn(), error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e))); NULL
  })
  sec <- (proc.time() - t0)[3]
  if (!is.null(out)) cat(sprintf("%.3fs\n", sec)) else sec <- NA_real_
  rm(out); gc(verbose = FALSE)
  sec
}

# -- Loop over configs --------------------------------------------------------
for (cfg_name in names(configs)) {
  cfg <- configs[[cfg_name]]
  sep <- paste(rep("-", 60), collapse = "")
  cat(sprintf("%s\n%s: %s\n%s\n", sep, cfg_name, cfg$desc, sep))

  # Build crossbasis once per config (reused across all backends)
  cat("  crossbasis... ")
  gc(verbose = FALSE)
  t0 <- proc.time()
  cb <- crossbasis(x, lag = cfg$lag, argvar = cfg$argvar, arglag = cfg$arglag,
                   group = group)
  cb_sec <- (proc.time() - t0)[3]
  cat(sprintf("%.3fs  (%d cols, %.0f MB)\n", cb_sec, ncol(cb),
              object.size(cb) / 1e6))
  cb_cols <- ncol(cb)

  # ---- glm (base R) ---------------------------------------------------------
  cat(sprintf("  glm      (base)... "))
  sec_glm <- run_backend("glm", function()
    glm(death ~ cb + ns(time, df = TIME_DF) + dow,
        data = df_model, family = quasipoisson()))

  # ---- bam ------------------------------------------------------------------
  cat(sprintf("  bam      (%d threads)... ", n_cores))
  sec_bam <- run_backend("bam", function()
    bam(death ~ cb + ns(time, df = TIME_DF) + dow,
        data = df_model, family = quasipoisson(),
        nthreads = n_cores))

  # ---- bam + discrete -------------------------------------------------------
  cat(sprintf("  bam_d    (%d threads, discrete)... ", n_cores))
  sec_bam_d <- run_backend("bam_d", function()
    bam(death ~ cb + ns(time, df = TIME_DF) + dow,
        data = df_model, family = quasipoisson(),
        discrete = TRUE, nthreads = n_cores))

  # ---- fastglm --------------------------------------------------------------
  sec_fg <- NA_real_
  if (has_fastglm) {
    cat("  fastglm  (includes model.matrix)... ")
    sec_fg <- run_backend("fastglm", function() {
      X <- model.matrix(~ cb + ns(time, df = TIME_DF) + dow, data = df_model)
      fastglm(X, df_model$death, family = quasipoisson())
    })
  }

  # Store results
  for (row in list(
    list("glm",          sec_glm),
    list("bam",          sec_bam),
    list("bam_discrete", sec_bam_d),
    list("fastglm",      sec_fg)
  )) {
    r_idx <- r_idx + 1L
    results[[r_idx]] <- data.frame(
      scale   = SCALE,
      config  = cfg_name,
      backend = row[[1]],
      n_rows  = n_rows,
      cb_cols = cb_cols,
      time_sec = row[[2]],
      stringsAsFactors = FALSE
    )
  }

  rm(cb); gc(verbose = FALSE)
  cat("\n")
}

# -- Save results -------------------------------------------------------------
results_dt <- rbindlist(results)
out_file   <- "benchmarks/results/glm_backends.csv"
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
fwrite(results_dt, out_file)
cat(sprintf("Saved to %s\n\n", out_file))

# -- Print comparison table ---------------------------------------------------
cat("=== GLM Backend Comparison (1GB, time_df=24) ===\n\n")

wide <- dcast(results_dt, config + cb_cols ~ backend, value.var = "time_sec")
# Add speedup columns vs glm baseline
for (col in c("bam", "bam_discrete", "fastglm")) {
  if (col %in% names(wide)) {
    wide[[paste0(col, "_speedup")]] <- round(wide$glm / wide[[col]], 2)
  }
}
print(wide)

# Also print a human-readable summary
cat("\n=== Speedup vs base glm ===\n")
for (cfg_name in names(configs)) {
  sub <- results_dt[config == cfg_name]
  base_t <- sub[backend == "glm", time_sec]
  cat(sprintf("\n%s (cb_cols=%d):\n", cfg_name, sub$cb_cols[1]))
  for (i in seq_len(nrow(sub))) {
    be  <- sub$backend[i]
    t   <- sub$time_sec[i]
    if (is.na(t)) {
      cat(sprintf("  %-14s  NA\n", be))
    } else {
      spd <- if (be == "glm") "baseline" else sprintf("%.2fx faster", base_t / t)
      cat(sprintf("  %-14s  %7.2fs  %s\n", be, t, spd))
    }
  }
}

cat("\nDone.\n")
