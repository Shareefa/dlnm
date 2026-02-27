#!/usr/bin/env Rscript
###############################################################################
# run_optimized_benchmarks.R -- Benchmark Rust-optimized DLNM pipeline
#
# Runs the same four pipeline stages as benchmark_dlnm.R but with the
# Rust-optimized crossbasis() and crosspred() functions:
#   1. crossbasis()  -- Rust P1 fused kernel
#   2. glm()         -- model fitting (quasipoisson)
#   3. crosspred()   -- Rust P2 quadratic form SE
#   4. crossreduce() -- overall cumulative reduction
#
# Also profiles internal sub-stages and computes speedup vs baseline.
#
# Usage: Rscript benchmarks/run_optimized_benchmarks.R [scales] [configs]
#   scales:  comma-separated, e.g. "10mb,100mb" (default: 10mb,100mb,1gb)
#   configs: comma-separated, e.g. "C1,C2" (default: see per-scale defaults)
#
# Scale-specific config defaults:
#   10mb, 100mb: C1,C2,C3,C4,C5
#   1gb:         C1,C2,C3
#   10gb:        C1,C2 (crossbasis-only if GLM won't fit)
#
# Prerequisites:
#   - Rust crate built (cd src/rust && cargo build --release)
#   - Benchmark data exists at required scales
#   - Baseline timing_results.csv exists
###############################################################################

# -- Dependencies -------------------------------------------------------------
required_pkgs <- c("data.table")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

library(data.table)
library(splines)

# Load dlnm from local source (includes Rust optimizations)
pkgload::load_all(".", quiet = TRUE)

# -- Verify Rust backend is available -----------------------------------------
rust_available <- tryCatch({
  res <- fused_crossbasis(
    matrix(as.numeric(1:6), nrow = 3, ncol = 2),
    matrix(as.numeric(1:2), nrow = 2, ncol = 1),
    0L, 1L, 1L, 3L
  )
  !is.null(res)
}, error = function(e) FALSE)

if (!rust_available) {
  stop("Rust backend is NOT available. Build with: cd src/rust && cargo build --release")
}
cat("Rust backend: AVAILABLE\n\n")

# -- Parse CLI args -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
all_scales <- c("10mb", "100mb", "1gb", "10gb")
all_configs <- c("C1", "C2", "C3", "C4", "C5")

# Default configs per scale
default_configs <- list(
  "10mb"  = c("C1", "C2", "C3", "C4", "C5"),
  "100mb" = c("C1", "C2", "C3", "C4", "C5"),
  "1gb"   = c("C1", "C2", "C3"),
  "10gb"  = c("C1", "C2")
)

if (length(args) >= 1 && nchar(args[1]) > 0) {
  scales_to_run <- tolower(unlist(strsplit(args[1], ",")))
} else {
  scales_to_run <- c("10mb", "100mb", "1gb")
}

user_configs <- NULL
if (length(args) >= 2 && nchar(args[2]) > 0) {
  user_configs <- toupper(unlist(strsplit(args[2], ",")))
}

cat("Scales:", paste(scales_to_run, collapse = ", "), "\n")
if (!is.null(user_configs)) {
  cat("Configs (user override):", paste(user_configs, collapse = ", "), "\n")
}
cat("\n")

# -- Output directory ---------------------------------------------------------
results_dir <- "benchmarks/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -- Model configurations -----------------------------------------------------
configs <- list(
  C1 = list(
    argvar = list(fun = "lin"),
    arglag = list(fun = "poly", degree = 4),
    lag = c(0L, 15L),
    cb_cols = 5L,
    desc = "Minimal/DLM: lin x poly(4), lag 0-15"
  ),
  C2 = list(
    argvar = list(fun = "ns", df = 5),
    arglag = list(fun = "ns", df = 4),
    lag = c(0L, 21L),
    cb_cols = 20L,
    desc = "Typical epi: ns(5) x ns(4), lag 0-21"
  ),
  C3 = list(
    argvar = list(fun = "bs", df = 6),
    arglag = list(fun = "ns", df = 4),
    lag = c(0L, 40L),
    cb_cols = 24L,
    desc = "Extended lag: bs(6) x ns(4), lag 0-40"
  ),
  C4 = list(
    argvar = list(fun = "ps", df = 10),
    arglag = list(fun = "ps", df = 5),
    lag = c(0L, 30L),
    cb_cols = 50L,
    desc = "Penalized: ps(10) x ps(5), lag 0-30"
  ),
  C5 = list(
    argvar = list(fun = "ps", df = 15),
    arglag = list(fun = "ps", df = 8),
    lag = c(0L, 60L),
    cb_cols = 120L,
    desc = "Stress test: ps(15) x ps(8), lag 0-60"
  )
)

# -- Helper: determine GLM time_df per scale ----------------------------------
get_time_df <- function(scale_name, n_years) {
  base_df <- 7 * n_years
  switch(scale_name,
    "10mb"  = min(base_df, 500),
    "100mb" = min(base_df, 500),
    "1gb"   = 24L,
    "10gb"  = 24L,
    min(base_df, 500)
  )
}

# -- Helper: load dataset (RDS or parquet for 10GB) ---------------------------
load_dataset <- function(scale_name) {
  if (scale_name == "10gb") {
    parquet_dir <- "benchmarks/data/scale_10gb"
    if (!dir.exists(parquet_dir)) {
      stop("10GB parquet directory not found: ", parquet_dir)
    }
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop("arrow package required for 10GB scale")
    }
    cat(sprintf("Loading %s via arrow... ", parquet_dir))
    t0 <- proc.time()
    dt <- as.data.table(arrow::read_parquet(
      list.files(parquet_dir, pattern = "\\.parquet$", full.names = TRUE)[1]
    ))
    # Check total row count
    ds <- arrow::open_dataset(parquet_dir)
    total_rows <- ds$num_rows
    cat(sprintf("%.1fs (chunk: %s rows, total dataset: %s rows)\n",
                (proc.time() - t0)[3],
                format(nrow(dt), big.mark = ","),
                format(total_rows, big.mark = ",")))
    attr(dt, "total_rows") <- total_rows
    attr(dt, "is_10gb") <- TRUE
    return(dt)
  }

  path <- file.path("benchmarks/data", paste0("scale_", scale_name, ".rds"))
  if (!file.exists(path)) {
    stop("Dataset not found: ", path, "\n  Run benchmarks/generate_data.R first")
  }
  cat(sprintf("Loading %s... ", path))
  t0 <- proc.time()
  dt <- readRDS(path)
  cat(sprintf("%.1fs (%s rows)\n", (proc.time() - t0)[3],
              format(nrow(dt), big.mark = ",")))
  attr(dt, "is_10gb") <- FALSE
  dt
}

# -- Helper: load 10GB dataset (only needed columns, minimal memory) ----------
load_10gb_full <- function() {
  parquet_dir <- "benchmarks/data/scale_10gb"
  if (!dir.exists(parquet_dir)) {
    stop("10GB parquet directory not found: ", parquet_dir)
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("arrow package required for 10GB scale")
  }
  cat("Loading 10GB dataset via arrow (only needed columns)...\n")
  t0 <- proc.time()
  ds <- arrow::open_dataset(parquet_dir)
  # Only select columns needed for benchmarking to minimize memory
  dt <- as.data.table(
    ds |>
      dplyr::select(temp, death, time, year, dow, city) |>
      dplyr::collect()
  )
  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  Loaded in %.1fs (%s rows, %.1f GB)\n", elapsed,
              format(nrow(dt), big.mark = ","),
              object.size(dt) / 1e9))
  dt
}

# -- Profile crossbasis internals (with Rust) ---------------------------------
profile_crossbasis_internals <- function(x, lag, argvar, arglag, group = NULL) {
  timings <- list()

  # Stage 1: onebasis for variable space
  gc(verbose = FALSE)
  t0 <- proc.time()
  basisvar <- do.call(onebasis, modifyList(argvar, list(x = as.numeric(x))))
  timings$onebasis_var <- (proc.time() - t0)[3]

  # Stage 2: onebasis for lag space
  gc(verbose = FALSE)
  al <- arglag
  if (length(al) == 0L || diff(lag) == 0L) {
    al <- list(fun = "strata", df = 1, intercept = TRUE)
  } else {
    if ((is.null(al$fun) || "intercept" %in% names(formals(al$fun))) &&
        sum(pmatch(names(al), "intercept", nomatch = 0)) == 0)
      al$intercept <- TRUE
    al$cen <- NULL
  }
  t0 <- proc.time()
  basislag <- do.call(onebasis, modifyList(al, list(x = seqlag(lag))))
  timings$onebasis_lag <- (proc.time() - t0)[3]

  # Stage 3: Rust fused kernel timing (combined lag_matrix + multiply)
  nv <- ncol(basisvar)
  nl <- ncol(basislag)
  n <- nrow(as.matrix(x))

  gc(verbose = FALSE)
  t0 <- proc.time()
  if (!is.null(group)) {
    grle <- rle(as.integer(group))
    gends <- cumsum(grle$lengths)
    gstarts <- c(1L, gends[-length(gends)] + 1L)
  } else {
    gstarts <- 1L
    gends <- as.integer(n)
  }
  rust_cb <- fused_crossbasis(
    unclass(basisvar), unclass(basislag),
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(gstarts), as.integer(gends)
  )
  timings$rust_fused_kernel <- (proc.time() - t0)[3]

  timings
}

# -- Profile crosspred internals (with Rust) ----------------------------------
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

  # Stage 3: Rust SE quadratic form
  gc(verbose = FALSE)
  t0 <- proc.time()
  matse_vec <- as.numeric(quad_form_se(Xpred, vcov_cb))
  matse <- matrix(matse_vec, length(at), length(predlag))
  timings$rust_se_quadform <- (proc.time() - t0)[3]

  # Stage 4: Rust cumulative SE
  predlag_int <- seqlag(lag)
  Xpred2 <- mkXpred("cb", cb, at, at, predlag_int, cen = NULL)

  gc(verbose = FALSE)
  t0 <- proc.time()
  n_lag <- length(predlag_int)
  res <- cumulative_quad_form_se(Xpred2, vcov_cb, as.integer(length(at)),
                                  as.integer(n_lag))
  timings$rust_cumulative_se <- (proc.time() - t0)[3]

  timings
}

# -- Main benchmark loop ------------------------------------------------------
timing_results <- list()
substage_results <- list()
result_idx <- 0L
sub_idx <- 0L

for (scale_name in scales_to_run) {
  sep <- paste(rep("=", 70), collapse = "")
  cat(sprintf("\n%s\n", sep))
  cat(sprintf("SCALE: %s (Rust-optimized)\n", toupper(scale_name)))
  cat(sprintf("%s\n", sep))

  # Determine configs to run for this scale
  configs_to_run <- if (!is.null(user_configs)) {
    user_configs
  } else {
    default_configs[[scale_name]]
  }
  if (is.null(configs_to_run)) configs_to_run <- all_configs

  cat("Configs:", paste(configs_to_run, collapse = ", "), "\n")

  # Load dataset
  if (scale_name == "10gb") {
    dt <- tryCatch(load_10gb_full(), error = function(e) {
      cat("  LOAD ERROR:", conditionMessage(e), "\n")
      NULL
    })
  } else {
    dt <- tryCatch(load_dataset(scale_name), error = function(e) {
      cat("  SKIPPED:", conditionMessage(e), "\n")
      NULL
    })
  }
  if (is.null(dt)) next

  n_rows <- nrow(dt)
  n_years <- length(unique(dt$year))
  time_df <- get_time_df(scale_name, n_years)

  x <- dt$temp
  death <- dt$death
  time_var <- dt$time
  dow <- dt$dow
  group <- dt$city

  # Build model data frame once per scale
  df_model <- data.frame(death = death, time = time_var, dow = dow)

  for (config_name in configs_to_run) {
    cfg <- configs[[config_name]]
    cat(sprintf("\n--- Config %s: %s ---\n", config_name, cfg$desc))

    # Check if this config is feasible for this scale
    cb_mem_gb <- as.numeric(n_rows) * cfg$cb_cols * 8 / 1e9
    if (cb_mem_gb > 28) {
      cat(sprintf("  SKIPPED: cross-basis would need %.1f GB (>28 GB limit)\n",
                  cb_mem_gb))
      next
    }
    cat(sprintf("  Estimated CB matrix: %.2f GB\n", cb_mem_gb))

    # For 10GB, we limit to crossbasis-only if GLM model matrix is too large
    total_model_cols <- cfg$cb_cols + time_df + 6L  # cb + time_spline + dow
    model_matrix_gb <- as.numeric(n_rows) * total_model_cols * 8 / 1e9
    crossbasis_only <- (scale_name == "10gb" && model_matrix_gb > 28)
    if (crossbasis_only) {
      cat(sprintf("  GLM model matrix: %.1f GB -- crossbasis-only mode\n",
                  model_matrix_gb))
    }

    # ---- Stage 1: crossbasis() (Rust-optimized) ----
    cat("  Stage 1: crossbasis() [Rust]\n")

    gc(verbose = FALSE)
    cb <- NULL
    cb_time <- tryCatch({
      t0 <- proc.time()
      cb <- crossbasis(x, lag = cfg$lag, argvar = cfg$argvar, arglag = cfg$arglag,
                       group = group)
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs  (cols: %d)\n", "crossbasis", elapsed,
                  ncol(cb)))

      # Run 1 additional iteration for 10mb/100mb scales for more stable timing
      if (scale_name %in% c("10mb", "100mb")) {
        gc(verbose = FALSE)
        t1 <- proc.time()
        crossbasis(x, lag = cfg$lag, argvar = cfg$argvar, arglag = cfg$arglag,
                   group = group)
        elapsed2 <- (proc.time() - t1)[3]
        median(c(elapsed, elapsed2))
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

    result_idx <- result_idx + 1L
    timing_results[[result_idx]] <- data.frame(
      scale = scale_name, config = config_name, stage = "crossbasis",
      n_rows = n_rows, cb_cols = ncol(cb),
      median_time_sec = cb_time,
      stringsAsFactors = FALSE
    )

    # ---- Sub-function profiling: crossbasis internals ----
    cat("  Profiling crossbasis internals:\n")
    cb_subs <- tryCatch({
      profile_crossbasis_internals(x, cfg$lag, cfg$argvar, cfg$arglag,
                                    group = group)
    }, error = function(e) {
      cat(sprintf("    Profile ERROR: %s\n", conditionMessage(e)))
      NULL
    })
    if (!is.null(cb_subs)) {
      for (nm in names(cb_subs)) {
        cat(sprintf("    %-25s %8.3fs\n", nm, cb_subs[[nm]]))
        sub_idx <- sub_idx + 1L
        substage_results[[sub_idx]] <- data.frame(
          scale = scale_name, config = config_name, substage = nm,
          time_seconds = cb_subs[[nm]], stringsAsFactors = FALSE
        )
      }
    }

    # Skip GLM/crosspred/crossreduce for crossbasis-only mode
    if (crossbasis_only) {
      cat("  Skipping GLM/crosspred/crossreduce (crossbasis-only mode)\n")
      rm(cb)
      gc(verbose = FALSE)
      next
    }

    # ---- Stage 2: glm() ----
    cat(sprintf("  Stage 2: glm() [time_df=%d]\n", time_df))

    model <- NULL
    glm_time <- tryCatch({
      gc(verbose = FALSE)
      t0 <- proc.time()
      model <- glm(death ~ cb + ns(time, df = time_df) + dow,
                    data = df_model, family = quasipoisson())
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs  (converged=%s, iter=%d)\n",
                  "glm", elapsed, model$converged, model$iter))
      elapsed
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      NA_real_
    })

    result_idx <- result_idx + 1L
    timing_results[[result_idx]] <- data.frame(
      scale = scale_name, config = config_name, stage = "glm",
      n_rows = n_rows, cb_cols = ncol(cb),
      median_time_sec = glm_time,
      stringsAsFactors = FALSE
    )

    if (is.null(model)) {
      cat("  Skipping prediction stages (model fit failed)\n")
      rm(cb)
      gc(verbose = FALSE)
      next
    }

    # ---- Stage 3: crosspred() (Rust-optimized SE) ----
    cat("  Stage 3: crosspred() [Rust SE]\n")
    pred <- NULL
    pred_time <- tryCatch({
      gc(verbose = FALSE)
      t0 <- proc.time()
      pred <- crosspred(cb, model,
                        at = seq(attr(cb, "range")[1], attr(cb, "range")[2],
                                 length.out = 50),
                        cumul = TRUE)
      elapsed <- (proc.time() - t0)[3]
      cat(sprintf("    %-25s %8.3fs\n", "crosspred", elapsed))
      elapsed
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", conditionMessage(e)))
      NA_real_
    })

    result_idx <- result_idx + 1L
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
          sub_idx <- sub_idx + 1L
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

    result_idx <- result_idx + 1L
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
  rm(dt, df_model, x, death, time_var, dow, group)
  gc(verbose = FALSE)
}

# -- Save timing results ------------------------------------------------------
cat("\n=== Saving Results ===\n")

timing_file <- file.path(results_dir, "optimized_timing_results.csv")
substage_file <- file.path(results_dir, "optimized_substage_timings.csv")

if (length(timing_results) > 0) {
  new_timing_df <- rbindlist(timing_results)

  # Merge with existing results (replace rows for re-run scale/config combos)
  if (file.exists(timing_file)) {
    existing <- fread(timing_file)
    # Remove old rows for scales+configs we just ran
    keys <- unique(new_timing_df[, .(scale, config)])
    for (i in seq_len(nrow(keys))) {
      existing <- existing[!(scale == keys$scale[i] & config == keys$config[i])]
    }
    timing_df <- rbind(existing, new_timing_df, fill = TRUE)
  } else {
    timing_df <- new_timing_df
  }

  fwrite(timing_df, timing_file)
  cat("Optimized timing results:", timing_file, "\n")
  print(timing_df)
} else {
  timing_df <- if (file.exists(timing_file)) fread(timing_file) else data.table()
}

if (length(substage_results) > 0) {
  new_substage_df <- rbindlist(substage_results)

  # Merge with existing
  if (file.exists(substage_file)) {
    existing_sub <- fread(substage_file)
    keys <- unique(new_substage_df[, .(scale, config)])
    for (i in seq_len(nrow(keys))) {
      existing_sub <- existing_sub[!(scale == keys$scale[i] & config == keys$config[i])]
    }
    substage_df <- rbind(existing_sub, new_substage_df, fill = TRUE)
  } else {
    substage_df <- new_substage_df
  }

  fwrite(substage_df, substage_file)
  cat("\nOptimized sub-stage timings:", substage_file, "\n")
  print(substage_df)
}

# -- Compute speedup vs baseline ----------------------------------------------
baseline_file <- file.path(results_dir, "timing_results.csv")
if (file.exists(baseline_file) && length(timing_results) > 0) {
  cat("\n=== Speedup Comparison ===\n")

  baseline_df <- fread(baseline_file)

  # Merge baseline and optimized on (scale, config, stage)
  comp <- merge(
    baseline_df[, .(scale, config, stage, r_time_sec = median_time_sec)],
    timing_df[, .(scale, config, stage, rust_time_sec = median_time_sec)],
    by = c("scale", "config", "stage"),
    all = FALSE
  )

  if (nrow(comp) > 0) {
    comp[, speedup_factor := round(r_time_sec / rust_time_sec, 2)]
    comp <- comp[order(scale, config, stage)]

    speedup_file <- file.path(results_dir, "speedup_comparison.csv")
    fwrite(comp, speedup_file)
    cat("Speedup comparison:", speedup_file, "\n")
    print(comp)

    # Highlight key result: C2 1GB crossbasis speedup
    c2_1gb <- comp[scale == "1gb" & config == "C2" & stage == "crossbasis"]
    if (nrow(c2_1gb) > 0) {
      cat(sprintf("\n*** KEY RESULT: C2 1GB crossbasis ***\n"))
      cat(sprintf("    R baseline:     %.3fs\n", c2_1gb$r_time_sec))
      cat(sprintf("    Rust optimized: %.3fs\n", c2_1gb$rust_time_sec))
      cat(sprintf("    Speedup:        %.2fx\n", c2_1gb$speedup_factor))
      if (c2_1gb$rust_time_sec < 21) {
        cat("    STATUS: PASS (< 21 seconds target)\n")
      } else {
        cat("    STATUS: WARN (>= 21 seconds target)\n")
      }
    }
  } else {
    cat("No matching baseline rows for speedup comparison.\n")
  }
} else {
  cat("\nBaseline timing_results.csv not found -- skipping speedup comparison.\n")
}

# -- Scaling analysis ---------------------------------------------------------
if (length(timing_results) > 0) {
  cat("\n=== Scaling Analysis (Optimized) ===\n")
  for (cfg_name in unique(timing_df$config)) {
    for (stg in unique(timing_df$stage)) {
      subset <- timing_df[config == cfg_name & stage == stg &
                           !is.na(median_time_sec)]
      if (nrow(subset) >= 2) {
        subset <- subset[order(n_rows)]
        cat(sprintf("\n%s / %s:\n", cfg_name, stg))
        for (i in seq_len(nrow(subset))) {
          cat(sprintf("  %6s: %8.3fs  (%s rows)\n",
                      subset$scale[i], subset$median_time_sec[i],
                      format(subset$n_rows[i], big.mark = ",")))
        }
        if (nrow(subset) >= 2) {
          log_n <- log10(subset$n_rows)
          log_t <- log10(pmax(1e-6, subset$median_time_sec))
          fit <- lm(log_t ~ log_n)
          alpha <- coef(fit)[2]
          cat(sprintf("  Scaling exponent (alpha): %.2f\n", alpha))
        }
      }
    }
  }
}

cat("\nOptimized benchmark complete.\n")
