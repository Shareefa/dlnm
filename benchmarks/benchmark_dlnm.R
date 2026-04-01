#!/usr/bin/env Rscript
###############################################################################
# benchmark_dlnm.R -- Benchmark DLNM pipeline at multiple scales/configs
#
# Benchmarks four stages independently:
#   1. crossbasis()  -- cross-basis matrix construction
#   2. glm()         -- model fitting (quasipoisson)
#   3. crosspred()   -- predictions with cumul=TRUE
#   4. crossreduce() -- overall cumulative reduction
#
# Also profiles internal sub-functions within crossbasis() and crosspred().
#
# Usage: Rscript benchmarks/benchmark_dlnm.R [scales] [configs]
#   scales:  comma-separated, e.g. "10mb,100mb" (default: 10mb,100mb,1gb)
#   configs: comma-separated, e.g. "C1,C2" (default: all)
#
# Prerequisites:
#   - Run benchmarks/generate_data.R first to create datasets
#   - Packages: dlnm (local), bench, data.table
###############################################################################

# -- Dependencies -------------------------------------------------------------
required_pkgs <- c("bench", "data.table")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

library(bench)
library(data.table)
library(splines)

# Load dlnm from local source
pkgload::load_all(".", quiet = TRUE)

# -- Parse CLI args -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
all_scales <- c("10mb", "100mb", "1gb")
all_configs <- c("C1", "C2", "C3", "C4", "C5")

if (length(args) >= 1 && nchar(args[1]) > 0) {
  scales_to_run <- tolower(unlist(strsplit(args[1], ",")))
} else {
  scales_to_run <- all_scales
}
if (length(args) >= 2 && nchar(args[2]) > 0) {
  configs_to_run <- toupper(unlist(strsplit(args[2], ",")))
} else {
  configs_to_run <- all_configs
}

cat("Scales:", paste(scales_to_run, collapse = ", "), "\n")
cat("Configs:", paste(configs_to_run, collapse = ", "), "\n\n")

# -- Output directory ---------------------------------------------------------
results_dir <- "benchmarks/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -- Model configurations -----------------------------------------------------
# Each config defines: argvar, arglag, lag, n_pred_values, description
configs <- list(
  C1 = list(
    argvar = list(fun = "lin"),
    arglag = list(fun = "poly", degree = 4),
    lag = c(0, 15),
    cb_cols = 5,
    desc = "Minimal/DLM: lin x poly(4), lag 0-15"
  ),
  C2 = list(
    argvar = list(fun = "ns", df = 5),
    arglag = list(fun = "ns", df = 4),
    lag = c(0, 21),
    cb_cols = 20,
    desc = "Typical epi: ns(5) x ns(4), lag 0-21"
  ),
  C3 = list(
    argvar = list(fun = "bs", df = 6),
    arglag = list(fun = "ns", df = 4),
    lag = c(0, 40),
    cb_cols = 24,
    desc = "Extended lag: bs(6) x ns(4), lag 0-40"
  ),
  C4 = list(
    argvar = list(fun = "ps", df = 10),
    arglag = list(fun = "ps", df = 5),
    lag = c(0, 30),
    cb_cols = 50,
    desc = "Penalized: ps(10) x ps(5), lag 0-30"
  ),
  C5 = list(
    argvar = list(fun = "ps", df = 15),
    arglag = list(fun = "ps", df = 8),
    lag = c(0, 60),
    cb_cols = 120,
    desc = "Stress test: ps(15) x ps(8), lag 0-60"
  )
)

# -- Helper: safe timing with GC tracking ------------------------------------
time_stage <- function(expr, label, n_iter = 3) {
  # Force GC before measurement
  gc(verbose = FALSE, reset = TRUE)

  # Use bench::mark for reliable timing
  result <- tryCatch({
    bm <- bench::mark(
      expr,
      iterations = n_iter,
      check = FALSE,
      memory = TRUE,
      filter_gc = FALSE,
      min_time = 0
    )
    list(
      median_sec  = as.numeric(bm$median),
      mem_alloc   = as.numeric(bm$mem_alloc) / 1e6,  # MB
      n_iters     = bm$n_itr,
      n_gc        = bm$n_gc,
      min_sec     = as.numeric(bm$min),
      max_sec     = as.numeric(bm$max)
    )
  }, error = function(e) {
    # Fallback: simple system.time
    gc(verbose = FALSE, reset = TRUE)
    mem_before <- gc(verbose = FALSE)[2, 2]  # Vcells MB
    times <- system.time(expr)
    mem_after <- gc(verbose = FALSE)[2, 2]
    list(
      median_sec  = times[3],
      mem_alloc   = max(0, mem_after - mem_before),
      n_iters     = 1,
      n_gc        = NA,
      min_sec     = times[3],
      max_sec     = times[3]
    )
  })

  cat(sprintf("    %-25s %8.3fs  (mem: %.1f MB)\n",
              label, result$median_sec, result$mem_alloc))
  result$label <- label
  result
}

# -- Helper: profile crossbasis internals -------------------------------------
profile_crossbasis_internals <- function(x, lag, argvar, arglag, group = NULL) {
  timings <- list()

  # Stage 1: onebasis for variable space
  gc(verbose = FALSE)
  t0 <- proc.time()
  basisvar <- do.call(onebasis, modifyList(argvar, list(x = as.numeric(x))))
  timings$onebasis_var <- (proc.time() - t0)[3]

  # Stage 2: onebasis for lag space
  gc(verbose = FALSE)
  t0 <- proc.time()
  basislag <- do.call(onebasis, modifyList(
    if (length(arglag) == 0L || diff(lag) == 0L)
      list(fun = "strata", df = 1, intercept = TRUE)
    else {
      al <- arglag
      if ((is.null(al$fun) || "intercept" %in% names(formals(al$fun))) &&
          sum(pmatch(names(al), "intercept", nomatch = 0)) == 0)
        al$intercept <- TRUE
      al$cen <- NULL
      al
    },
    list(x = seqlag(lag))
  ))
  timings$onebasis_lag <- (proc.time() - t0)[3]

  # Stage 3: Lag matrix construction + cross-basis multiply
  nv <- ncol(basisvar)
  nl <- ncol(basislag)
  n <- nrow(as.matrix(x))

  gc(verbose = FALSE)
  t0 <- proc.time()
  lag_matrices <- vector("list", nv)
  for (v in seq_len(nv)) {
    lag_matrices[[v]] <- as.matrix(tsModel::Lag(basisvar[, v], seqlag(lag), group = group))
  }
  timings$lag_matrix <- (proc.time() - t0)[3]

  # Stage 4: Matrix multiply loop
  gc(verbose = FALSE)
  t0 <- proc.time()
  cb <- matrix(0, nrow = n, ncol = nv * nl)
  for (v in seq_len(nv)) {
    mat <- lag_matrices[[v]]
    for (l in seq_len(nl)) {
      cb[, nl * (v - 1) + l] <- mat %*% basislag[, l]
    }
  }
  timings$multiply_loop <- (proc.time() - t0)[3]

  timings
}

# -- Helper: profile crosspred internals --------------------------------------
profile_crosspred_internals <- function(cb, model, n_at = 50) {
  timings <- list()

  coef_full <- coef(model)
  vcov_full <- vcov(model)
  name <- "cb"

  # Extract relevant coefficients
  cond <- paste0(name, "[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
  ind <- grep(cond, names(coef_full))
  coef_cb <- coef_full[ind]
  vcov_cb <- vcov_full[ind, ind, drop = FALSE]

  # Set up prediction values
  range_x <- attr(cb, "range")
  at <- seq(range_x[1], range_x[2], length.out = n_at)
  lag <- attr(cb, "lag")
  predlag <- seqlag(lag)

  # Stage 1: mkXpred (includes tensor product)
  gc(verbose = FALSE)
  t0 <- proc.time()
  Xpred <- mkXpred("cb", cb, at, at, predlag, cen = NULL)
  timings$mkXpred <- (proc.time() - t0)[3]

  # Stage 2: Xpred %*% coef (prediction)
  gc(verbose = FALSE)
  t0 <- proc.time()
  matfit <- matrix(Xpred %*% coef_cb, length(at), length(predlag))
  timings$prediction <- (proc.time() - t0)[3]

  # Stage 3: SE quadratic form  rowSums((Xpred %*% vcov) * Xpred)
  gc(verbose = FALSE)
  t0 <- proc.time()
  matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% vcov_cb) * Xpred))),
                  length(at), length(predlag))
  timings$se_quadform <- (proc.time() - t0)[3]

  # Re-create for overall/cumulative (integer lags only)
  predlag_int <- seqlag(lag)
  Xpred2 <- mkXpred("cb", cb, at, at, predlag_int, cen = NULL)

  # Stage 4: Cumulative accumulation loop
  gc(verbose = FALSE)
  t0 <- proc.time()
  Xpredall <- 0
  cumfit <- cumse <- matrix(0, length(at), length(predlag_int))
  for (i in seq_along(predlag_int)) {
    idx <- seq(length(at)) + length(at) * (i - 1)
    Xpredall <- Xpredall + Xpred2[idx, , drop = FALSE]
    cumfit[, i] <- Xpredall %*% coef_cb
    cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% vcov_cb) * Xpredall)))
  }
  timings$cumulative_loop <- (proc.time() - t0)[3]

  timings
}

# -- Load dataset helper ------------------------------------------------------
load_dataset <- function(scale_name) {
  path <- file.path("benchmarks/data", paste0("scale_", scale_name, ".rds"))
  if (!file.exists(path)) {
    stop("Dataset not found: ", path, "\n  Run benchmarks/generate_data.R first")
  }
  cat(sprintf("Loading %s... ", path))
  t0 <- proc.time()
  dt <- readRDS(path)
  cat(sprintf("%.1fs (%s rows)\n", (proc.time() - t0)[3],
              format(nrow(dt), big.mark = ",")))
  dt
}

# -- Determine iteration count based on scale ---------------------------------
get_n_iter <- function(scale_name, stage) {
  base <- switch(scale_name,
    "10mb"  = 5,
    "100mb" = 3,
    "1gb"   = 1,
    1
  )
  # Fewer iterations for slower stages
  if (stage %in% c("glm", "crosspred")) base <- max(1, base - 1)
  base
}

# -- Main benchmark loop ------------------------------------------------------
timing_results <- list()
substage_results <- list()
result_idx <- 0
sub_idx <- 0

for (scale_name in scales_to_run) {
  sep <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n", sep))
  cat(sprintf("SCALE: %s\n", toupper(scale_name)))
  cat(sprintf("%s\n", sep))

  # Load dataset
  dt <- tryCatch(load_dataset(scale_name), error = function(e) {
    cat("  SKIPPED:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(dt)) next

  n_rows <- nrow(dt)
  n_years <- length(unique(dt$year))

  for (config_name in configs_to_run) {
    cfg <- configs[[config_name]]
    cat(sprintf("\n--- Config %s: %s ---\n", config_name, cfg$desc))

    # Check if this config is feasible for this scale
    # Cross-basis matrix size: n_rows * cb_cols * 8 bytes
    cb_mem_gb <- n_rows * cfg$cb_cols * 8 / 1e9
    if (cb_mem_gb > 32) {
      cat(sprintf("  SKIPPED: cross-basis would need %.1f GB (>32 GB limit)\n", cb_mem_gb))
      next
    }
    cat(sprintf("  Estimated CB matrix: %.2f GB\n", cb_mem_gb))

    # Impute NA in pm10 (required for crossbasis)
    x <- dt$temp
    death <- dt$death
    time_var <- dt$time
    dow <- dt$dow
    group <- dt$city

    # ---- Stage 1: crossbasis() ----
    cat("  Stage 1: crossbasis()\n")
    n_iter <- get_n_iter(scale_name, "crossbasis")

    gc(verbose = FALSE)
    cb <- NULL
    cb_time <- tryCatch({
      t0 <- proc.time()
      cb <- crossbasis(x, lag = cfg$lag, argvar = cfg$argvar, arglag = cfg$arglag,
                       group = group)
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs  (cols: %d)\n", "crossbasis", elapsed, ncol(cb)))

      # Run additional iterations for timing
      if (n_iter > 1) {
        times <- numeric(n_iter - 1)
        for (it in seq_len(n_iter - 1)) {
          gc(verbose = FALSE)
          t0 <- proc.time()
          crossbasis(x, lag = cfg$lag, argvar = cfg$argvar, arglag = cfg$arglag,
                     group = group)
          times[it] <- (proc.time() - t0)[3]
        }
        median(c(elapsed, times))
      } else {
        elapsed
      }
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      NA_real_
    })

    if (is.null(cb) || is.na(cb_time)) {
      cat("  Skipping remaining stages (crossbasis failed)\n")
      next
    }

    result_idx <- result_idx + 1
    timing_results[[result_idx]] <- data.frame(
      scale = scale_name, config = config_name, stage = "crossbasis",
      n_rows = n_rows, cb_cols = ncol(cb),
      median_time_sec = cb_time,
      stringsAsFactors = FALSE
    )

    # ---- Sub-function profiling: crossbasis internals ----
    cat("  Profiling crossbasis internals:\n")
    cb_subs <- tryCatch({
      profile_crossbasis_internals(x, cfg$lag, cfg$argvar, cfg$arglag, group = group)
    }, error = function(e) {
      cat(sprintf("    Profile ERROR: %s\n", conditionMessage(e)))
      NULL
    })
    if (!is.null(cb_subs)) {
      for (nm in names(cb_subs)) {
        cat(sprintf("    %-25s %8.3fs\n", nm, cb_subs[[nm]]))
        sub_idx <- sub_idx + 1
        substage_results[[sub_idx]] <- data.frame(
          scale = scale_name, config = config_name, substage = nm,
          time_seconds = cb_subs[[nm]], stringsAsFactors = FALSE
        )
      }
    }

    # ---- Stage 2: glm() ----
    cat("  Stage 2: glm()\n")
    n_iter_glm <- get_n_iter(scale_name, "glm")

    # Build formula: death ~ cb + ns(time, 7*n_years) + dow
    # For multi-city: include city as factor strata
    model <- NULL
    glm_time <- tryCatch({
      df_model <- data.frame(
        death = death,
        time = time_var,
        dow = dow,
        city = group
      )
      # Number of df for time spline: 7 per year
      time_df <- 7 * n_years

      # Cap time_df to avoid extreme memory usage
      time_df <- min(time_df, 500)

      gc(verbose = FALSE)
      t0 <- proc.time()
      model <- glm(death ~ cb + ns(time, df = time_df) + dow,
                    data = df_model, family = quasipoisson())
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs\n", "glm", elapsed))
      elapsed
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      NA_real_
    })

    if (is.null(model)) {
      cat("  Skipping prediction stages (model fit failed)\n")

      result_idx <- result_idx + 1
      timing_results[[result_idx]] <- data.frame(
        scale = scale_name, config = config_name, stage = "glm",
        n_rows = n_rows, cb_cols = ncol(cb),
        median_time_sec = glm_time,
        stringsAsFactors = FALSE
      )
      next
    }

    result_idx <- result_idx + 1
    timing_results[[result_idx]] <- data.frame(
      scale = scale_name, config = config_name, stage = "glm",
      n_rows = n_rows, cb_cols = ncol(cb),
      median_time_sec = glm_time,
      stringsAsFactors = FALSE
    )

    # ---- Stage 3: crosspred() ----
    cat("  Stage 3: crosspred()\n")
    pred <- NULL
    pred_time <- tryCatch({
      gc(verbose = FALSE)
      t0 <- proc.time()
      pred <- crosspred(cb, model, at = seq(
        attr(cb, "range")[1], attr(cb, "range")[2], length.out = 50),
        cumul = TRUE)
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs\n", "crosspred", elapsed))
      elapsed
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      NA_real_
    })

    result_idx <- result_idx + 1
    timing_results[[result_idx]] <- data.frame(
      scale = scale_name, config = config_name, stage = "crosspred",
      n_rows = n_rows, cb_cols = ncol(cb),
      median_time_sec = pred_time,
      stringsAsFactors = FALSE
    )

    # ---- Sub-function profiling: crosspred internals ----
    if (!is.null(model)) {
      cat("  Profiling crosspred internals:\n")
      cp_subs <- tryCatch({
        profile_crosspred_internals(cb, model, n_at = 50)
      }, error = function(e) {
        cat(sprintf("    Profile ERROR: %s\n", conditionMessage(e)))
        NULL
      })
      if (!is.null(cp_subs)) {
        for (nm in names(cp_subs)) {
          cat(sprintf("    %-25s %8.3fs\n", nm, cp_subs[[nm]]))
          sub_idx <- sub_idx + 1
          substage_results[[sub_idx]] <- data.frame(
            scale = scale_name, config = config_name, substage = nm,
            time_seconds = cp_subs[[nm]], stringsAsFactors = FALSE
          )
        }
      }
    }

    # ---- Stage 4: crossreduce() ----
    cat("  Stage 4: crossreduce()\n")
    reduce_time <- tryCatch({
      gc(verbose = FALSE)
      t0 <- proc.time()
      cr <- crossreduce(cb, model, type = "overall")
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs\n", "crossreduce", elapsed))
      elapsed
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      NA_real_
    })

    result_idx <- result_idx + 1
    timing_results[[result_idx]] <- data.frame(
      scale = scale_name, config = config_name, stage = "crossreduce",
      n_rows = n_rows, cb_cols = ncol(cb),
      median_time_sec = reduce_time,
      stringsAsFactors = FALSE
    )

    # Clean up per-config
    rm(cb, model, pred)
    gc(verbose = FALSE)
    cat("\n")
  }

  # Clean up per-scale
  rm(dt)
  gc(verbose = FALSE)
}

# -- Save results -------------------------------------------------------------
cat("\n=== Saving Results ===\n")

if (length(timing_results) > 0) {
  timing_df <- rbindlist(timing_results)
  timing_file <- file.path(results_dir, "timing_results.csv")
  fwrite(timing_df, timing_file)
  cat("Timing results:", timing_file, "\n")
  print(timing_df)
}

if (length(substage_results) > 0) {
  substage_df <- rbindlist(substage_results)
  substage_file <- file.path(results_dir, "substage_timings.csv")
  fwrite(substage_df, substage_file)
  cat("\nSub-stage timings:", substage_file, "\n")
  print(substage_df)
}

# -- Scaling analysis ---------------------------------------------------------
if (length(timing_results) > 0) {
  cat("\n=== Scaling Analysis ===\n")
  scaling_df <- timing_df[!is.na(median_time_sec)]
  scaling_file <- file.path(results_dir, "scaling_analysis.csv")
  fwrite(scaling_df, scaling_file)

  # Print scaling summary: for each (config, stage), show time vs n_rows
  for (cfg_name in unique(scaling_df$config)) {
    for (stg in unique(scaling_df$stage)) {
      subset <- scaling_df[config == cfg_name & stage == stg]
      if (nrow(subset) >= 2) {
        cat(sprintf("\n%s / %s:\n", cfg_name, stg))
        for (i in seq_len(nrow(subset))) {
          cat(sprintf("  %6s: %8.3fs  (%s rows)\n",
                      subset$scale[i], subset$median_time_sec[i],
                      format(subset$n_rows[i], big.mark = ",")))
        }
        # Estimate scaling exponent: time ~ n^alpha
        if (nrow(subset) >= 2) {
          log_n <- log10(subset$n_rows)
          log_t <- log10(pmax(1e-6, subset$median_time_sec))
          fit <- lm(log_t ~ log_n)
          alpha <- coef(fit)[2]
          cat(sprintf("  Scaling exponent (alpha): %.2f  (linear=1.0, quadratic=2.0)\n", alpha))
        }
      }
    }
  }
  cat("\nScaling analysis:", scaling_file, "\n")
}

cat("\nBenchmark complete.\n")
