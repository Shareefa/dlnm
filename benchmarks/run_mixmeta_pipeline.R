#!/usr/bin/env Rscript
###############################################################################
# run_mixmeta_pipeline.R -- Two-Stage DLNM + Mixmeta Pipeline
#
# Implements the full two-stage workflow:
#   Stage 1: Per-city DLNM (crossbasis -> glm -> crossreduce)
#   Stage 2: Pool with mixmeta (multivariate random-effects meta-analysis)
#
# Outputs:
#   - benchmarks/results/mixmeta_timing.csv    (per-city and total timing)
#   - benchmarks/results/mixmeta_estimates.csv  (pooled coefficients + heterogeneity)
#
# Usage:
#   Rscript benchmarks/run_mixmeta_pipeline.R [scale] [config]
#     scale:  10mb or 100mb (default: 10mb)
#     config: C1, C2, C3, C4, C5 (default: C2)
#
# Notes:
#   - Forces R fallback (unoptimized) for baseline timing
#   - Supports appending results across multiple runs
###############################################################################

# -- Dependencies -------------------------------------------------------------
required_pkgs <- c("mixmeta", "data.table", "splines")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

library(mixmeta)
library(data.table)
library(splines)

# Load dlnm from local source
# Note: Uses whatever backend is available (Rust if compiled, R fallback otherwise).
# Per-city fitting calls crossbasis without group= since each city is a single
# time series, so the pipeline works correctly with either backend.
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

cat(sprintf("=== Mixmeta Pipeline: scale=%s, config=%s ===\n\n", scale_name, config_name))

# -- Model configurations -----------------------------------------------------
configs <- list(
  C1 = list(
    argvar = list(fun = "lin"),
    arglag = list(fun = "poly", degree = 4),
    lag = c(0, 15),
    desc = "Minimal/DLM: lin x poly(4), lag 0-15"
  ),
  C2 = list(
    argvar = list(fun = "ns", df = 5),
    arglag = list(fun = "ns", df = 4),
    lag = c(0, 21),
    desc = "Typical epi: ns(5) x ns(4), lag 0-21"
  ),
  C3 = list(
    argvar = list(fun = "bs", df = 6),
    arglag = list(fun = "ns", df = 4),
    lag = c(0, 40),
    desc = "Extended lag: bs(6) x ns(4), lag 0-40"
  ),
  C4 = list(
    argvar = list(fun = "ps", df = 10),
    arglag = list(fun = "ps", df = 5),
    lag = c(0, 30),
    desc = "Penalized: ps(10) x ps(5), lag 0-30"
  ),
  C5 = list(
    argvar = list(fun = "ps", df = 15),
    arglag = list(fun = "ps", df = 8),
    lag = c(0, 60),
    desc = "Stress test: ps(15) x ps(8), lag 0-60"
  )
)

if (!config_name %in% names(configs)) {
  stop("Unknown config: ", config_name,
       ". Choose from: ", paste(names(configs), collapse = ", "))
}

cfg <- configs[[config_name]]
cat(sprintf("Config %s: %s\n", config_name, cfg$desc))

# -- Output directory ---------------------------------------------------------
results_dir <- "benchmarks/results"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# -- Load dataset -------------------------------------------------------------
data_path <- file.path("benchmarks/data", paste0("scale_", scale_name, ".rds"))
if (!file.exists(data_path)) {
  stop("Dataset not found: ", data_path,
       "\nRun benchmarks/generate_data.R first")
}

cat(sprintf("\nLoading %s... ", data_path))
t_load <- proc.time()
dt <- readRDS(data_path)
cat(sprintf("%.1fs (%s rows)\n", (proc.time() - t_load)[3],
            format(nrow(dt), big.mark = ",")))

cities <- unique(dt$city)
n_cities <- length(cities)
cat(sprintf("Cities: %d\n\n", n_cities))

# -- Stage 1: Per-City DLNM --------------------------------------------------
cat("=== STAGE 1: Per-City DLNM ===\n\n")

city_coefs <- list()
city_vcovs <- list()
city_timings <- list()
city_converged <- logical(n_cities)
names(city_converged) <- cities

total_stage1_start <- proc.time()

for (i in seq_along(cities)) {
  city_name <- cities[i]
  city_start <- proc.time()

  city_dt <- dt[dt$city == city_name, ]
  n_years <- length(unique(city_dt$year))
  # 7 df per year for time trend, capped at 100
  time_df <- min(7 * n_years, 100)

  # Step 1a: Build cross-basis for this city (single time series, no group)
  t_cb <- proc.time()
  cb <- crossbasis(city_dt$temp, lag = cfg$lag,
                   argvar = cfg$argvar, arglag = cfg$arglag)
  cb_time <- (proc.time() - t_cb)[3]

  # Step 1b: Fit GLM
  t_glm <- proc.time()
  model <- glm(death ~ cb + ns(time, df = time_df) + dow,
               data = city_dt, family = quasipoisson())
  glm_time <- (proc.time() - t_glm)[3]

  city_converged[i] <- model$converged

  # Step 1c: Crossreduce to overall coefficients + vcov
  t_cr <- proc.time()
  cr <- crossreduce(cb, model, type = "overall")
  cr_time <- (proc.time() - t_cr)[3]

  city_coefs[[city_name]] <- coef(cr)
  city_vcovs[[city_name]] <- cr$vcov

  city_total <- (proc.time() - city_start)[3]

  city_timings[[i]] <- data.frame(
    scale = scale_name,
    config = config_name,
    city = city_name,
    n_rows = nrow(city_dt),
    cb_time_sec = cb_time,
    glm_time_sec = glm_time,
    cr_time_sec = cr_time,
    total_time_sec = city_total,
    converged = model$converged,
    stringsAsFactors = FALSE
  )

  # Progress indicator every 10 cities or at start/end
  if (i <= 3 || i %% 10 == 0 || i == n_cities) {
    cat(sprintf("  City %d/%d (%s): cb=%.3fs, glm=%.3fs, cr=%.3fs, total=%.3fs%s\n",
                i, n_cities, city_name, cb_time, glm_time, cr_time, city_total,
                if (!model$converged) " [NOT CONVERGED]" else ""))
  }
}

total_stage1_time <- (proc.time() - total_stage1_start)[3]
n_converged <- sum(city_converged)

cat(sprintf("\nStage 1 complete: %d/%d cities converged, total time: %.2fs\n\n",
            n_converged, n_cities, total_stage1_time))

# -- Stage 2: Pool with Mixmeta ----------------------------------------------
cat("=== STAGE 2: Mixmeta Pooling ===\n\n")

# Assemble coefficient matrix and vcov list
coef_matrix <- do.call(rbind, city_coefs)
vcov_list <- city_vcovs

cat(sprintf("Coefficient matrix: %d cities x %d coefficients\n",
            nrow(coef_matrix), ncol(coef_matrix)))

# Fit multivariate random-effects meta-analysis
t_pool <- proc.time()
pooled <- tryCatch({
  mixmeta(coef_matrix ~ 1, S = vcov_list, method = "reml")
}, error = function(e) {
  cat(sprintf("  mixmeta REML failed: %s\n", conditionMessage(e)))
  cat("  Trying method='ml'...\n")
  tryCatch({
    mixmeta(coef_matrix ~ 1, S = vcov_list, method = "ml")
  }, error = function(e2) {
    cat(sprintf("  mixmeta ML also failed: %s\n", conditionMessage(e2)))
    NULL
  })
})
pool_time <- (proc.time() - t_pool)[3]

if (is.null(pooled)) {
  stop("Mixmeta pooling failed for both REML and ML methods")
}

cat(sprintf("Pooling time: %.3fs\n", pool_time))

# Extract pooled estimates
pooled_coef <- coef(pooled)
pooled_vcov <- vcov(pooled)
pooled_se <- sqrt(diag(pooled_vcov))
pooled_summary <- summary(pooled)

cat(sprintf("Pooled coefficients: %s\n",
            paste(sprintf("%.4f", pooled_coef), collapse = ", ")))
cat(sprintf("Pooled SEs: %s\n",
            paste(sprintf("%.4f", pooled_se), collapse = ", ")))

# -- Heterogeneity Statistics -------------------------------------------------
cat("\n=== HETEROGENEITY STATISTICS ===\n\n")

# Between-study variance (Psi)
psi <- pooled$Psi
cat("Between-study variance-covariance (Psi, diagonal):\n")
cat(sprintf("  %s\n", paste(sprintf("%.6f", diag(psi)), collapse = ", ")))

# Tau-squared (diagonal of Psi = between-study variances for each coefficient)
tau2 <- diag(psi)
cat(sprintf("\nTau-squared (between-city variance per coefficient):\n"))
for (j in seq_along(tau2)) {
  cat(sprintf("  coef %d: tau2 = %.6f, tau = %.4f\n", j, tau2[j], sqrt(tau2[j])))
}

# Cochran Q test for heterogeneity
q_test <- tryCatch({
  qtest(pooled)
}, error = function(e) {
  cat(sprintf("  Q-test error: %s\n", conditionMessage(e)))
  NULL
})

i_squared <- NULL
q_values <- NULL
if (!is.null(q_test)) {
  cat(sprintf("\nCochran Q test:\n"))
  q_values <- q_test$Q
  p_values <- q_test$pvalue
  df_values <- q_test$df

  # Compute I-squared manually: I2 = max(0, (Q - df) / Q * 100)
  # For multivariate models, compute per-dimension
  # Skip the ".all" entry (first element), use per-outcome values
  if (length(q_values) > 1) {
    # Per-outcome Q/df (skip .all)
    per_outcome_q <- q_values[-1]
    per_outcome_df <- df_values[-1]
    i_squared <- pmax(0, (per_outcome_q - per_outcome_df) / per_outcome_q * 100)
  } else {
    # Univariate case
    i_squared <- pmax(0, (q_values - df_values) / q_values * 100)
    per_outcome_q <- q_values
    per_outcome_df <- df_values
  }

  # Print overall Q
  cat(sprintf("  Overall: Q = %.2f, df = %d, p = %.4f\n",
              q_values[1], df_values[1], p_values[1]))

  # Print per-dimension Q and I-squared
  outcome_names <- names(per_outcome_q)
  if (is.null(outcome_names)) outcome_names <- seq_along(per_outcome_q)
  for (j in seq_along(per_outcome_q)) {
    cat(sprintf("  Coef %s: Q = %.2f, df = %d, p = %.4f, I-squared = %.1f%%\n",
                outcome_names[j], per_outcome_q[j], per_outcome_df[j],
                p_values[j + if (length(q_values) > 1) 1 else 0],
                i_squared[j]))
  }
}

# -- Prediction Intervals (Pooled Exposure-Response) --------------------------
cat("\n=== POOLED EXPOSURE-RESPONSE ESTIMATES ===\n\n")

# Reconstruct the pooled exposure-response curve
# Use onebasis with the same argvar to create prediction basis
temp_range <- range(dt$temp, na.rm = TRUE)
pred_temps <- seq(temp_range[1], temp_range[2], length.out = 50)

# Build the variable basis at prediction temperatures
pred_basis <- do.call(onebasis, c(list(x = pred_temps), cfg$argvar))

# Center at median temperature
med_temp <- median(dt$temp, na.rm = TRUE)
cen_basis <- do.call(onebasis, c(list(x = med_temp), cfg$argvar))
pred_basis_cen <- scale(pred_basis, center = cen_basis, scale = FALSE)

# Compute pooled overall log-RR and SE
pooled_logRR <- as.vector(pred_basis_cen %*% pooled_coef)
pooled_logRR_se <- sqrt(pmax(0, rowSums((pred_basis_cen %*% pooled_vcov) * pred_basis_cen)))

# Compute RR and 95% CI
z <- qnorm(0.975)
pooled_RR <- exp(pooled_logRR)
pooled_RR_low <- exp(pooled_logRR - z * pooled_logRR_se)
pooled_RR_high <- exp(pooled_logRR + z * pooled_logRR_se)

# Print key temperature points
key_temps <- c(temp_range[1], quantile(dt$temp, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), temp_range[2])
key_temps <- sort(unique(round(key_temps, 1)))

cat("Pooled cumulative RR at key temperatures (centered at median):\n")
for (kt in key_temps) {
  idx <- which.min(abs(pred_temps - kt))
  cat(sprintf("  %6.1f°C: RR = %6.3f (95%% CI: %6.3f - %6.3f)\n",
              pred_temps[idx], pooled_RR[idx], pooled_RR_low[idx], pooled_RR_high[idx]))
}

# -- Total Timing -------------------------------------------------------------
total_time <- total_stage1_time + pool_time
cat(sprintf("\n=== TOTAL TIME: %.2fs (Stage1: %.2fs, Stage2: %.2fs) ===\n\n",
            total_time, total_stage1_time, pool_time))

# -- Save Timing Results ------------------------------------------------------
cat("=== SAVING RESULTS ===\n\n")

# Per-city timing
city_timing_df <- rbindlist(city_timings)

# Summary timing (stage-level)
timing_summary <- data.frame(
  scale = scale_name,
  config = config_name,
  stage = c("stage1_per_city", "stage2_mixmeta_pool", "total"),
  n_cities = n_cities,
  time_sec = c(total_stage1_time, pool_time, total_time),
  n_converged = c(n_converged, NA, n_converged),
  stringsAsFactors = FALSE
)

# Append to existing file or create new
timing_file <- file.path(results_dir, "mixmeta_timing.csv")
if (file.exists(timing_file)) {
  existing <- fread(timing_file)
  # Remove old rows for same scale+config
  existing <- existing[!(scale == scale_name & config == config_name)]
  timing_out <- rbind(existing, as.data.table(timing_summary), fill = TRUE)
} else {
  timing_out <- as.data.table(timing_summary)
}
fwrite(timing_out, timing_file)
cat(sprintf("Timing saved to: %s\n", timing_file))

# Also save per-city timing
percity_file <- file.path(results_dir, "mixmeta_percity_timing.csv")
if (file.exists(percity_file)) {
  existing_pc <- fread(percity_file)
  existing_pc <- existing_pc[!(scale == scale_name & config == config_name)]
  percity_out <- rbind(existing_pc, as.data.table(city_timing_df), fill = TRUE)
} else {
  percity_out <- as.data.table(city_timing_df)
}
fwrite(percity_out, percity_file)
cat(sprintf("Per-city timing saved to: %s\n", percity_file))

# -- Save Estimates -----------------------------------------------------------
# Pooled coefficient estimates with heterogeneity statistics
estimates_rows <- list()
for (j in seq_along(pooled_coef)) {
  row <- data.frame(
    scale = scale_name,
    config = config_name,
    coefficient_index = j,
    pooled_coef = pooled_coef[j],
    pooled_se = pooled_se[j],
    tau2 = tau2[j],
    tau = sqrt(tau2[j]),
    I2 = if (!is.null(i_squared) && j <= length(i_squared)) i_squared[j] else NA_real_,
    Q = if (!is.null(q_values) && length(q_values) > 1 && j < length(q_values)) q_values[j + 1] else if (!is.null(q_values) && j <= length(q_values)) q_values[j] else NA_real_,
    stringsAsFactors = FALSE
  )
  estimates_rows[[j]] <- row
}
estimates_df <- rbindlist(estimates_rows)

# Also add pooled RR at key temperature points
rr_rows <- list()
for (ki in seq_along(key_temps)) {
  kt <- key_temps[ki]
  idx <- which.min(abs(pred_temps - kt))
  rr_rows[[ki]] <- data.frame(
    scale = scale_name,
    config = config_name,
    temperature = pred_temps[idx],
    pooled_RR = pooled_RR[idx],
    pooled_RR_low = pooled_RR_low[idx],
    pooled_RR_high = pooled_RR_high[idx],
    stringsAsFactors = FALSE
  )
}
rr_df <- rbindlist(rr_rows)

# Save coefficient estimates
estimates_file <- file.path(results_dir, "mixmeta_estimates.csv")
if (file.exists(estimates_file)) {
  existing_est <- fread(estimates_file)
  existing_est <- existing_est[!(scale == scale_name & config == config_name)]
  estimates_out <- rbind(existing_est, as.data.table(estimates_df), fill = TRUE)
} else {
  estimates_out <- as.data.table(estimates_df)
}
fwrite(estimates_out, estimates_file)
cat(sprintf("Pooled estimates saved to: %s\n", estimates_file))

# Save RR estimates
rr_file <- file.path(results_dir, "mixmeta_rr_estimates.csv")
if (file.exists(rr_file)) {
  existing_rr <- fread(rr_file)
  existing_rr <- existing_rr[!(scale == scale_name & config == config_name)]
  rr_out <- rbind(existing_rr, as.data.table(rr_df), fill = TRUE)
} else {
  rr_out <- as.data.table(rr_df)
}
fwrite(rr_out, rr_file)
cat(sprintf("RR estimates saved to: %s\n", rr_file))

cat("\n=== Pipeline Complete ===\n")
