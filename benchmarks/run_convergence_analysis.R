#!/usr/bin/env Rscript
###############################################################################
# run_convergence_analysis.R -- Compare single-stage DLNM vs two-stage
#                               DLNM + mixmeta on truth data with known DGP
#
# Runs BOTH approaches on the SAME truth data file, then compares each
# against the known TRUE exposure-lag-response surface.
#
# === APPROACHES ===
#
# 1) SINGLE-STAGE optimized DLNM (Rust-enabled):
#    crossbasis(group=city) -> glm -> crosspred -> extract estimated RR surface
#
# 2) TWO-STAGE unoptimized DLNM + mixmeta (R fallback):
#    Per-city: crossbasis -> glm -> crossreduce
#    Pool with mixmeta (REML random-effects)
#    Extract pooled RR surface
#
# === METRICS ===
#
# Statistical accuracy (vs TRUE surface):
#   - Bias = estimated_RR - true_RR at each of 50 temperature points
#   - MSE = mean(bias^2) across temperatures
#   - RMSE = sqrt(MSE)
#   - Relative bias = bias / |true_RR| at 5th, 50th, 95th percentile temps
#   - CI coverage = proportion of 50 temps where true RR in 95% CI
#
# Numerical convergence:
#   - GLM iteration count, convergence flag
#   - Total timing
#
# === OUTPUTS ===
#
# benchmarks/results/convergence_metrics.csv
#   Columns: approach, config, scale, bias, mse, rmse, rel_bias, ci_coverage
#
# benchmarks/results/convergence_details.csv
#   Columns: approach, config, temperature, estimated_rr, true_rr, se,
#            ci_lower, ci_upper
#
# benchmarks/results/numerical_convergence.csv
#   Columns: approach, config, scale, iterations, converged, time_sec
#
# Usage:
#   Rscript benchmarks/run_convergence_analysis.R [scale] [config]
#     scale:  10mb, 100mb, 1gb (default: 10mb)
#     config: C2 (default) — must match truth data DGP config
#
###############################################################################

# -- Dependencies -------------------------------------------------------------
required_pkgs <- c("data.table", "splines", "jsonlite", "mixmeta")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

library(data.table)
library(splines)
library(jsonlite)
library(mixmeta)

# Load dlnm from local source (includes Rust optimizations if available)
pkgload::load_all(".", quiet = TRUE)

# -- Parse CLI args -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

scale_name <- if (length(args) >= 1 && nchar(args[1]) > 0) {
  tolower(args[1])
} else {
  "10mb"
}

config_name <- if (length(args) >= 2 && nchar(args[2]) > 0) {
  toupper(args[2])
} else {
  "C2"
}

valid_scales <- c("10mb", "100mb", "1gb")
if (!scale_name %in% valid_scales) {
  stop("Unknown scale: ", scale_name, ". Choose from: ", paste(valid_scales, collapse = ", "))
}

cat(sprintf("=== Convergence Analysis: scale=%s, config=%s ===\n\n",
            scale_name, config_name))

# -- Output directory ---------------------------------------------------------
results_dir <- "benchmarks/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -- Load truth data and specification ----------------------------------------
truth_data_path <- file.path("benchmarks/data", paste0("truth_", scale_name, ".rds"))
truth_spec_path <- "benchmarks/data/truth_spec.json"

if (!file.exists(truth_data_path)) {
  stop("Truth data not found: ", truth_data_path,
       "\nRun: Rscript benchmarks/generate_truth_data.R ", scale_name)
}
if (!file.exists(truth_spec_path)) {
  stop("Truth spec not found: ", truth_spec_path,
       "\nRun: Rscript benchmarks/generate_truth_data.R")
}

cat(sprintf("Loading truth data: %s\n", truth_data_path))
t_load <- proc.time()
dt <- readRDS(truth_data_path)
load_time <- (proc.time() - t_load)[3]
cat(sprintf("  Loaded in %.1fs: %s rows, %d cities\n",
            load_time, format(nrow(dt), big.mark = ","),
            length(unique(dt$city))))

cat(sprintf("Loading truth spec: %s\n", truth_spec_path))
spec <- fromJSON(truth_spec_path)

# SAME truth data file used for BOTH approaches
cat(sprintf("BOTH approaches will use: %s\n\n", truth_data_path))

# -- Extract true coefficients from spec --------------------------------------
beta_true <- unlist(spec$true_coefficients$beta_true)
n_var_basis <- spec$crossbasis_config$n_var_basis
n_lag_basis <- spec$crossbasis_config$n_lag_basis

# Cross-basis configuration (C2: ns(5) x ns(4), lag 0-21)
cb_argvar <- list(fun = spec$crossbasis_config$argvar$fun,
                  df = spec$crossbasis_config$argvar$df)
cb_arglag <- list(fun = spec$crossbasis_config$arglag$fun,
                  df = spec$crossbasis_config$arglag$df)
cb_lag <- unlist(spec$crossbasis_config$lag)

# Time trend
n_time_basis <- 7  # as specified in DGP

cat(sprintf("True beta (%d coefficients): %s\n",
            length(beta_true),
            paste(sprintf("%.4f", beta_true[1:min(5, length(beta_true))]), collapse = ", ")))
cat(sprintf("Cross-basis: %s(df=%d) x %s(df=%d), lag %d-%d\n\n",
            cb_argvar$fun, cb_argvar$df, cb_arglag$fun, cb_arglag$df,
            cb_lag[1], cb_lag[2]))

# -- Compute TRUE overall cumulative RR surface --------------------------------
cat("=== COMPUTING TRUE RR SURFACE ===\n\n")

temp_range <- range(dt$temp)
pred_temps <- seq(temp_range[1], temp_range[2], length.out = 50)
cen_temp <- 20  # Center near the optimum temperature

# Build variable basis at prediction temperatures
pred_var_basis <- do.call(onebasis, c(list(x = pred_temps), cb_argvar))
cen_var_basis <- do.call(onebasis, c(list(x = cen_temp), cb_argvar))
pred_var_cen <- scale(pred_var_basis, center = cen_var_basis, scale = FALSE)

# Build lag basis and compute lag sum for overall reduction
lag_seq <- seqlag(cb_lag)
pred_lag_basis <- do.call(onebasis, c(list(x = lag_seq), cb_arglag))
lag_sum <- colSums(pred_lag_basis)

# Overall reduced coefficient: for each var basis v,
# reduced_coef_v = sum_l(beta_true[v,l] * sum_j(basislag[j,l]))
reduced_coef_true <- numeric(n_var_basis)
for (v in seq_len(n_var_basis)) {
  for (l in seq_len(n_lag_basis)) {
    reduced_coef_true[v] <- reduced_coef_true[v] +
      beta_true[(v - 1) * n_lag_basis + l] * lag_sum[l]
  }
}

true_log_rr <- as.numeric(unclass(pred_var_cen) %*% reduced_coef_true)
true_rr <- exp(true_log_rr)

cat(sprintf("True RR range: %.4f - %.4f (centered at %.1fC)\n",
            min(true_rr), max(true_rr), cen_temp))
cat(sprintf("Prediction temperatures: %.1f to %.1fC (%d points)\n\n",
            min(pred_temps), max(pred_temps), length(pred_temps)))

# -- Identify key percentile temperatures ------------------------------------
p05_idx <- which.min(abs(pred_temps - quantile(dt$temp, 0.05)))
p50_idx <- which.min(abs(pred_temps - quantile(dt$temp, 0.50)))
p95_idx <- which.min(abs(pred_temps - quantile(dt$temp, 0.95)))

cat(sprintf("Percentile temperatures: P5=%.1fC, P50=%.1fC, P95=%.1fC\n\n",
            pred_temps[p05_idx], pred_temps[p50_idx], pred_temps[p95_idx]))

# ===========================================================================
# APPROACH 1: SINGLE-STAGE OPTIMIZED DLNM (Rust-enabled)
# ===========================================================================
cat("=== APPROACH 1: SINGLE-STAGE DLNM ===\n\n")

t_single <- proc.time()

# Step 1: Build crossbasis with group (uses Rust fused kernel if available)
cat("  Building crossbasis (group=city)... ")
t_cb <- proc.time()
cb_single <- crossbasis(dt$temp, lag = cb_lag,
                        argvar = cb_argvar, arglag = cb_arglag,
                        group = dt$city)
cb_time <- (proc.time() - t_cb)[3]
cat(sprintf("%.2fs\n", cb_time))

# Step 2: Fit GLM
cat("  Fitting GLM... ")
t_glm <- proc.time()
model_single <- glm(death ~ cb_single + ns(time, df = n_time_basis) + dow,
                    data = dt, family = quasipoisson())
glm_time <- (proc.time() - t_glm)[3]
cat(sprintf("%.2fs (converged: %s, iterations: %d)\n",
            glm_time, model_single$converged, model_single$iter))

single_iterations <- model_single$iter
single_converged <- model_single$converged

# Step 3: Predictions with crosspred
cat("  Computing crosspred... ")
t_cp <- proc.time()
cp_single <- crosspred(cb_single, model_single, at = pred_temps, cen = cen_temp)
cp_time <- (proc.time() - t_cp)[3]
cat(sprintf("%.2fs\n", cp_time))

single_total_time <- (proc.time() - t_single)[3]
cat(sprintf("  Total single-stage time: %.2fs\n\n", single_total_time))

# Extract estimated RR and SE
single_rr <- cp_single$allRRfit
single_log_rr <- log(single_rr)
single_se <- cp_single$allse  # SE of log(RR)
z <- qnorm(0.975)
single_ci_lower <- exp(single_log_rr - z * single_se)
single_ci_upper <- exp(single_log_rr + z * single_se)

# ===========================================================================
# APPROACH 2: TWO-STAGE DLNM + MIXMETA (R fallback)
# ===========================================================================
cat("=== APPROACH 2: TWO-STAGE DLNM + MIXMETA ===\n\n")

t_two <- proc.time()
cities <- unique(dt$city)
n_cities <- length(cities)

city_coefs <- list()
city_vcovs <- list()
city_iters <- rep(0L, n_cities)
city_converged <- rep(FALSE, n_cities)
names(city_iters) <- as.character(cities)
names(city_converged) <- as.character(cities)

# Stage 1: Per-city DLNM
cat(sprintf("  Stage 1: Fitting %d per-city models...\n", n_cities))
t_stage1 <- proc.time()

for (i in seq_along(cities)) {
  city_name <- cities[i]
  city_dt <- dt[dt$city == city_name, ]
  n_years <- length(unique(city_dt$year))
  time_df <- min(7 * n_years, 100)

  # Build crossbasis (single time series, no group)
  cb_city <- crossbasis(city_dt$temp, lag = cb_lag,
                        argvar = cb_argvar, arglag = cb_arglag)

  # Fit GLM
  model_city <- glm(death ~ cb_city + ns(time, df = time_df) + dow,
                    data = city_dt, family = quasipoisson())

  city_iters[i] <- model_city$iter
  city_converged[i] <- model_city$converged

  # Crossreduce to overall coefficients + vcov
  cr <- crossreduce(cb_city, model_city, type = "overall")
  city_coefs[[city_name]] <- coef(cr)
  city_vcovs[[city_name]] <- cr$vcov

  # Progress
  if (i <= 3 || i %% max(1, n_cities %/% 10) == 0 || i == n_cities) {
    cat(sprintf("    City %d/%d: iter=%d, converged=%s\n",
                i, n_cities, model_city$iter, model_city$converged))
  }
}

stage1_time <- (proc.time() - t_stage1)[3]
n_converged_cities <- sum(city_converged)
cat(sprintf("  Stage 1 done: %d/%d converged, %.2fs\n\n",
            n_converged_cities, n_cities, stage1_time))

# Stage 2: Pool with mixmeta
cat("  Stage 2: Mixmeta pooling... ")
t_stage2 <- proc.time()

coef_matrix <- do.call(rbind, city_coefs)
vcov_list <- city_vcovs

pooled <- tryCatch({
  mixmeta(coef_matrix ~ 1, S = vcov_list, method = "reml")
}, error = function(e) {
  cat(sprintf("REML failed: %s, trying ML\n", conditionMessage(e)))
  tryCatch({
    mixmeta(coef_matrix ~ 1, S = vcov_list, method = "ml")
  }, error = function(e2) {
    cat(sprintf("ML also failed: %s\n", conditionMessage(e2)))
    NULL
  })
})

stage2_time <- (proc.time() - t_stage2)[3]
two_total_time <- (proc.time() - t_two)[3]

if (is.null(pooled)) {
  stop("Mixmeta pooling failed")
}

cat(sprintf("%.2fs\n", stage2_time))
cat(sprintf("  Total two-stage time: %.2fs\n\n", two_total_time))

# Extract pooled RR surface
# Build variable basis at prediction temps (centered)
pooled_coef <- coef(pooled)
pooled_vcov <- vcov(pooled)

pooled_log_rr <- as.numeric(unclass(pred_var_cen) %*% pooled_coef)
pooled_se <- sqrt(pmax(0, rowSums((unclass(pred_var_cen) %*% pooled_vcov) *
                                    unclass(pred_var_cen))))
two_rr <- exp(pooled_log_rr)
two_ci_lower <- exp(pooled_log_rr - z * pooled_se)
two_ci_upper <- exp(pooled_log_rr + z * pooled_se)

# Two-stage convergence: use median city iterations, all converged flag
two_iterations <- as.integer(median(city_iters))
two_converged <- all(city_converged)

# ===========================================================================
# COMPUTE METRICS: Compare both approaches against TRUE surface
# ===========================================================================
cat("=== COMPUTING CONVERGENCE & ACCURACY METRICS ===\n\n")

compute_metrics <- function(est_rr, est_ci_lower, est_ci_upper, est_se,
                            true_rr, approach_name, config, scale,
                            pred_temps, p05_idx, p50_idx, p95_idx) {
  # Bias at each temperature
  bias <- est_rr - true_rr

  # MSE and RMSE
  mse <- mean(bias^2)
  rmse <- sqrt(mse)

  # Mean bias across all temperatures
  mean_bias <- mean(bias)

  # Relative bias at percentile temperatures
  # Relative bias = bias / |true_value|
  rel_bias_p05 <- bias[p05_idx] / abs(true_rr[p05_idx])
  rel_bias_p50 <- bias[p50_idx] / abs(true_rr[p50_idx])
  rel_bias_p95 <- bias[p95_idx] / abs(true_rr[p95_idx])
  # Average relative bias across the three percentiles
  avg_rel_bias <- mean(c(rel_bias_p05, rel_bias_p50, rel_bias_p95))

  # CI coverage: proportion of temperature points where true RR falls in 95% CI
  ci_coverage <- mean(true_rr >= est_ci_lower & true_rr <= est_ci_upper)

  cat(sprintf("  %s:\n", approach_name))
  cat(sprintf("    Mean bias:      %+.6f\n", mean_bias))
  cat(sprintf("    MSE:            %.6f\n", mse))
  cat(sprintf("    RMSE:           %.6f\n", rmse))
  cat(sprintf("    Rel bias (P5):  %+.4f\n", rel_bias_p05))
  cat(sprintf("    Rel bias (P50): %+.4f\n", rel_bias_p50))
  cat(sprintf("    Rel bias (P95): %+.4f\n", rel_bias_p95))
  cat(sprintf("    Avg rel bias:   %+.4f\n", avg_rel_bias))
  cat(sprintf("    CI coverage:    %.2f (%.0f%% of %d temps)\n\n",
              ci_coverage, ci_coverage * 100, length(pred_temps)))

  # Return summary row
  list(
    metrics = data.frame(
      approach = approach_name,
      config = config,
      scale = scale,
      bias = mean_bias,
      mse = mse,
      rmse = rmse,
      rel_bias = avg_rel_bias,
      ci_coverage = ci_coverage,
      stringsAsFactors = FALSE
    ),
    details = data.frame(
      approach = approach_name,
      config = config,
      scale = scale,
      temperature = pred_temps,
      estimated_rr = est_rr,
      true_rr = true_rr,
      se = est_se,
      ci_lower = est_ci_lower,
      ci_upper = est_ci_upper,
      stringsAsFactors = FALSE
    )
  )
}

# Single-stage metrics
single_results <- compute_metrics(
  est_rr = single_rr,
  est_ci_lower = single_ci_lower,
  est_ci_upper = single_ci_upper,
  est_se = single_se,
  true_rr = true_rr,
  approach_name = "single_stage",
  config = config_name,
  scale = scale_name,
  pred_temps = pred_temps,
  p05_idx = p05_idx,
  p50_idx = p50_idx,
  p95_idx = p95_idx
)

# Two-stage metrics
two_results <- compute_metrics(
  est_rr = two_rr,
  est_ci_lower = two_ci_lower,
  est_ci_upper = two_ci_upper,
  est_se = pooled_se,
  true_rr = true_rr,
  approach_name = "two_stage",
  config = config_name,
  scale = scale_name,
  pred_temps = pred_temps,
  p05_idx = p05_idx,
  p50_idx = p50_idx,
  p95_idx = p95_idx
)

# ===========================================================================
# SAVE RESULTS TO CSV
# ===========================================================================
cat("=== SAVING RESULTS ===\n\n")

# -- convergence_metrics.csv --------------------------------------------------
metrics_df <- rbind(single_results$metrics, two_results$metrics)
metrics_file <- file.path(results_dir, "convergence_metrics.csv")

if (file.exists(metrics_file)) {
  existing <- fread(metrics_file)
  # Remove old rows for same scale+config
  existing <- existing[!(scale == scale_name & config == config_name)]
  metrics_out <- rbind(existing, as.data.table(metrics_df), fill = TRUE)
} else {
  metrics_out <- as.data.table(metrics_df)
}
fwrite(metrics_out, metrics_file)
cat(sprintf("Saved: %s (%d rows)\n", metrics_file, nrow(metrics_out)))

# -- convergence_details.csv --------------------------------------------------
details_df <- rbind(single_results$details, two_results$details)
details_file <- file.path(results_dir, "convergence_details.csv")

if (file.exists(details_file)) {
  existing <- fread(details_file)
  existing <- existing[!(scale == scale_name & config == config_name)]
  details_out <- rbind(existing, as.data.table(details_df), fill = TRUE)
} else {
  details_out <- as.data.table(details_df)
}
fwrite(details_out, details_file)
cat(sprintf("Saved: %s (%d rows)\n", details_file, nrow(details_out)))

# -- numerical_convergence.csv ------------------------------------------------
conv_df <- data.frame(
  approach = c("single_stage", "two_stage"),
  config = config_name,
  scale = scale_name,
  iterations = c(single_iterations, two_iterations),
  converged = c(single_converged, two_converged),
  time_sec = c(single_total_time, two_total_time),
  stringsAsFactors = FALSE
)
conv_file <- file.path(results_dir, "numerical_convergence.csv")

if (file.exists(conv_file)) {
  existing <- fread(conv_file)
  existing <- existing[!(scale == scale_name & config == config_name)]
  conv_out <- rbind(existing, as.data.table(conv_df), fill = TRUE)
} else {
  conv_out <- as.data.table(conv_df)
}
fwrite(conv_out, conv_file)
cat(sprintf("Saved: %s (%d rows)\n", conv_file, nrow(conv_out)))

# -- Summary Table ------------------------------------------------------------
cat("\n=== SUMMARY TABLE ===\n\n")
cat(sprintf("%-14s | %-6s | %-5s | %10s | %10s | %10s | %10s | %10s\n",
            "Approach", "Config", "Scale", "Bias", "MSE", "RMSE",
            "Rel.Bias", "CI.Cover"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (i in seq_len(nrow(metrics_df))) {
  row <- metrics_df[i, ]
  cat(sprintf("%-14s | %-6s | %-5s | %+10.6f | %10.6f | %10.6f | %+10.4f | %10.4f\n",
              row$approach, row$config, row$scale,
              row$bias, row$mse, row$rmse, row$rel_bias, row$ci_coverage))
}

cat("\n=== NUMERICAL CONVERGENCE ===\n\n")
cat(sprintf("%-14s | %-6s | %-5s | %10s | %10s | %10s\n",
            "Approach", "Config", "Scale", "Iterations", "Converged", "Time(s)"))
cat(paste(rep("-", 75), collapse = ""), "\n")
for (i in seq_len(nrow(conv_df))) {
  row <- conv_df[i, ]
  cat(sprintf("%-14s | %-6s | %-5s | %10d | %10s | %10.2f\n",
              row$approach, row$config, row$scale,
              row$iterations, row$converged, row$time_sec))
}

cat(sprintf("\n=== Convergence Analysis Complete (scale=%s, config=%s) ===\n",
            scale_name, config_name))
