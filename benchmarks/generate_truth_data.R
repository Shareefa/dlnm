#!/usr/bin/env Rscript
###############################################################################
# generate_truth_data.R -- Generate synthetic data with a KNOWN true
#                          exposure-lag-response surface for convergence testing
#
# === DATA GENERATING PROCESS (DGP) SPECIFICATION ===
#
# Model: death_it ~ Poisson(mu_it)
#   log(mu_it) = intercept + cb_it %*% beta_true_i + ns(time_it, df=time_df) %*% gamma_true + dow_effects
#
# Where:
#   i = city index, t = day index within city
#   cb_it = crossbasis(temp_it, lag=c(0,21), argvar=list(fun="ns",df=5), arglag=list(fun="ns",df=4))
#   This produces 20 cross-basis columns (5 var basis x 4 lag basis)
#
# TRUE COEFFICIENTS:
#   intercept_true = 3.5 (baseline ~33 deaths/day = exp(3.5))
#
#   beta_true (20 values, C2: 5 var x 4 lag basis):
#     The coefficients are chosen to create a U-shaped exposure-response
#     (high risk at cold and hot extremes, lowest risk at ~20C) that
#     decays with increasing lag (strongest at lag 0, weak by lag 21).
#     Arranged as: for var basis v=1..5, lag basis l=1..4:
#       beta_true[(v-1)*4 + l]
#
#   gamma_true (7 values for ns(time, df=7)):
#     Smooth seasonal confounding trend
#     Values: c(0.05, -0.03, 0.04, -0.02, 0.03, -0.01, 0.02)
#
#   dow_effects (6 values, reference = Sunday):
#     Mon=-0.02, Tue=-0.03, Wed=-0.025, Thu=-0.02, Fri=-0.01, Sat=0.01
#
# BETWEEN-CITY HETEROGENEITY:
#   beta_true_city = beta_true + N(0, tau^2 * I)
#   where tau = 0.005 (small between-city variance)
#   intercept_city = intercept_true + N(0, 0.10^2)
#   (cities have different baseline mortality rates)
#
# TEMPERATURE GENERATION:
#   temp_it = mean_temp_city + amplitude_city * sin(2*pi*(doy-80)/365) + noise
#   mean_temp_city ~ N(15, 5^2)    (different climates)
#   amplitude_city ~ N(12, 2^2)    (different seasonal amplitudes)
#   noise ~ N(0, 3^2)              (daily variability)
#
# SCALE CONFIGURATION:
#   10MB:  ~24 cities   x 5114 days = ~123K rows
#   100MB: ~235 cities  x 5114 days = ~1.2M rows
#   1GB:   ~2350 cities x 5114 days = ~12M rows
#
###############################################################################

# -- Dependencies -------------------------------------------------------------
required_pkgs <- c("data.table", "splines", "jsonlite")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(data.table)
library(splines)
library(jsonlite)

# Load dlnm from source
pkgload::load_all(".", quiet = TRUE)

# -- Parse CLI args -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
valid_scales <- c("10mb", "100mb", "1gb")
if (length(args) >= 1 && nchar(args[1]) > 0) {
  requested <- tolower(args[1])
  if (!requested %in% valid_scales) {
    stop("Unknown scale: ", requested, ". Choose from: ", paste(valid_scales, collapse = ", "))
  }
  scales_to_generate <- requested
} else {
  scales_to_generate <- "10mb"
}

cat("=== Truth Data Generation ===\n")
cat("Scales:", paste(scales_to_generate, collapse = ", "), "\n\n")

# -- TRUE DGP PARAMETERS (MATHEMATICAL SPECIFICATION) -------------------------
set.seed(42)

# Intercept: exp(3.5) ~ 33 deaths/day baseline
intercept_true <- 3.5

# Cross-basis configuration (C2: ns(5) x ns(4), lag 0-21)
cb_argvar <- list(fun = "ns", df = 5)
cb_arglag <- list(fun = "ns", df = 4)
cb_lag <- c(0, 21)
n_var_basis <- 5
n_lag_basis <- 4
n_cb_coefs <- n_var_basis * n_lag_basis  # 20

# TRUE beta coefficients for the cross-basis (20 values)
# Structure: beta_true[(v-1)*n_lag_basis + l] for var basis v, lag basis l
#
# Design rationale:
#   - Variable basis columns capture different parts of the temperature range
#     (cold edge, cool, moderate, warm, hot edge due to ns spline placement)
#   - Lag basis columns capture different lag patterns (immediate, short, medium, long)
#   - We want: U-shape in temperature (high at extremes, low at moderate)
#     and decay with lag (strongest immediate effect)
#   - Effect sizes are chosen to be large enough for reliable recovery:
#     the cumulative log-RR at temperature extremes should be ~0.1-0.3
#     (RR of ~1.1-1.35), matching realistic temperature-mortality associations
#
# The coefficients are arranged as v1.l1, v1.l2, v1.l3, v1.l4, v2.l1, ...
beta_true <- c(
  # v1 (cold end of temp spectrum): strong immediate effect, moderate decay
   0.08,  0.05,  0.03,  0.015,
  # v2 (cool-moderate): smaller effects
   0.04,  0.025,  0.015,  0.008,
  # v3 (moderate - near optimum): near zero (minimum risk region)
  -0.01, -0.005, -0.003,  0.001,
  # v4 (warm): moderate positive effects
   0.05,  0.03,  0.02,  0.01,
  # v5 (hot end): strong immediate effect, faster decay
   0.10,  0.06,  0.03,  0.015
)
names(beta_true) <- paste0("v", rep(1:n_var_basis, each = n_lag_basis),
                           ".l", rep(1:n_lag_basis, n_var_basis))

# Time trend coefficients (ns with df=7, ~1 per year for 14-year span)
n_time_basis <- 7
gamma_true <- c(0.05, -0.03, 0.04, -0.02, 0.03, -0.01, 0.02)
names(gamma_true) <- paste0("time_ns", 1:n_time_basis)

# Day-of-week effects (6 coefficients, Sunday = reference)
dow_effects_true <- c(
  Monday = -0.020,
  Tuesday = -0.030,
  Wednesday = -0.025,
  Thursday = -0.020,
  Friday = -0.010,
  Saturday = 0.010
)

# Between-city heterogeneity
tau_beta <- 0.005  # SD for random effects on beta (small to ensure recovery)
tau_intercept <- 0.10  # SD for random effects on intercept

# Temperature generation parameters
temp_mean_global <- 15     # Global mean temperature (C)
temp_mean_sd <- 5          # Between-city SD of mean temperature
temp_amplitude_mean <- 12  # Mean seasonal amplitude
temp_amplitude_sd <- 2     # Between-city SD of seasonal amplitude
temp_noise_sd <- 3         # Daily temperature noise SD

# Time series parameters
n_days <- 5114  # 14 years (1987-2000), matching chicagoNMMAPS
start_date <- as.Date("1987-01-01")

# -- Scale definitions --------------------------------------------------------
scale_config <- list(
  "10mb"  = list(n_cities = 24),
  "100mb" = list(n_cities = 235),
  "1gb"   = list(n_cities = 2350)
)

# -- Data output directory ----------------------------------------------------
data_dir <- "benchmarks/data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -- Helper: generate one city's data -----------------------------------------
generate_city <- function(city_id, city_seed,
                          intercept_true, beta_true, gamma_true, dow_effects_true,
                          tau_beta, tau_intercept,
                          cb_argvar, cb_arglag, cb_lag,
                          n_days, start_date,
                          temp_mean_global, temp_mean_sd,
                          temp_amplitude_mean, temp_amplitude_sd,
                          temp_noise_sd, n_time_basis) {
  set.seed(city_seed)

  # -- City-specific parameters --
  intercept_city <- intercept_true + rnorm(1, 0, tau_intercept)
  beta_city <- beta_true + rnorm(length(beta_true), 0, tau_beta)
  # gamma and dow are shared across cities (confounders)

  # -- Generate dates and time variables --
  dates <- start_date + 0:(n_days - 1)
  doy <- as.integer(format(dates, "%j"))
  dow <- factor(weekdays(dates),
                levels = c("Sunday", "Monday", "Tuesday", "Wednesday",
                           "Thursday", "Friday", "Saturday"))
  time_idx <- seq_len(n_days)
  year <- as.integer(format(dates, "%Y"))

  # -- Generate temperature --
  mean_temp_city <- rnorm(1, temp_mean_global, temp_mean_sd)
  amplitude_city <- rnorm(1, temp_amplitude_mean, temp_amplitude_sd)
  temp <- mean_temp_city + amplitude_city * sin(2 * pi * (doy - 80) / 365) +
    rnorm(n_days, 0, temp_noise_sd)

  # -- Build crossbasis matrix --
  cb <- crossbasis(temp, lag = cb_lag, argvar = cb_argvar, arglag = cb_arglag)

  # -- Time trend basis --
  time_basis <- ns(time_idx, df = n_time_basis)

  # -- Compute linear predictor --
  # log(mu) = intercept + cb %*% beta + time_basis %*% gamma + dow_effects
  cb_effect <- as.numeric(unclass(cb) %*% beta_city)
  # The first lag[2] rows have NA in the crossbasis (incomplete lag window).
  # For data generation, set those to 0 (baseline effect = no temperature effect).
  cb_effect[is.na(cb_effect)] <- 0

  log_mu <- intercept_city + cb_effect +
    as.numeric(time_basis %*% gamma_true)

  # Add day-of-week effects
  dow_numeric <- as.integer(dow)  # 1=Sunday, 2=Monday, ...
  for (d in 2:7) {
    log_mu[dow_numeric == d] <- log_mu[dow_numeric == d] + dow_effects_true[d - 1]
  }

  # -- Generate death counts --
  mu <- exp(log_mu)
  # Cap extremely large mu to avoid numerical issues
  mu <- pmin(mu, 1e6)
  death <- rpois(n_days, lambda = mu)

  # -- Assemble data.table --
  dt <- data.table(
    city = city_id,
    date = dates,
    time = time_idx,
    year = year,
    doy = doy,
    dow = dow,
    temp = round(temp, 2),
    death = death
  )

  list(
    data = dt,
    beta_city = beta_city,
    intercept_city = intercept_city,
    mean_temp_city = mean_temp_city,
    amplitude_city = amplitude_city
  )
}

# -- Generate for each scale --------------------------------------------------
for (scale_name in scales_to_generate) {
  cfg <- scale_config[[scale_name]]
  n_cities <- cfg$n_cities

  cat(sprintf("=== Generating %s: %d cities x %d days = %s rows ===\n",
              toupper(scale_name), n_cities, n_days,
              format(n_cities * n_days, big.mark = ",")))

  t0 <- proc.time()

  all_data <- vector("list", n_cities)
  city_params <- vector("list", n_cities)

  for (i in seq_len(n_cities)) {
    city_id <- sprintf("city_%05d", i)
    city_seed <- 1000 + i

    result <- generate_city(
      city_id = city_id,
      city_seed = city_seed,
      intercept_true = intercept_true,
      beta_true = beta_true,
      gamma_true = gamma_true,
      dow_effects_true = dow_effects_true,
      tau_beta = tau_beta,
      tau_intercept = tau_intercept,
      cb_argvar = cb_argvar,
      cb_arglag = cb_arglag,
      cb_lag = cb_lag,
      n_days = n_days,
      start_date = start_date,
      temp_mean_global = temp_mean_global,
      temp_mean_sd = temp_mean_sd,
      temp_amplitude_mean = temp_amplitude_mean,
      temp_amplitude_sd = temp_amplitude_sd,
      temp_noise_sd = temp_noise_sd,
      n_time_basis = n_time_basis
    )

    all_data[[i]] <- result$data
    city_params[[i]] <- list(
      city = city_id,
      intercept = result$intercept_city,
      beta = result$beta_city,
      mean_temp = result$mean_temp_city,
      amplitude = result$amplitude_city
    )

    if (i <= 3 || i %% max(1, n_cities %/% 10) == 0 || i == n_cities) {
      cat(sprintf("  City %d/%d (%s): mean_temp=%.1fC, deaths mean=%.1f\n",
                  i, n_cities, city_id, result$mean_temp_city,
                  mean(result$data$death)))
    }
  }

  gen_time <- (proc.time() - t0)[3]
  cat(sprintf("  Data generation: %.1fs\n", gen_time))

  # Combine all cities
  dt <- rbindlist(all_data)
  dt[, city := as.factor(city)]

  # Save data
  outfile <- file.path(data_dir, paste0("truth_", scale_name, ".rds"))
  t1 <- proc.time()
  saveRDS(dt, outfile)
  save_time <- (proc.time() - t1)[3]

  file_size <- file.size(outfile)
  cat(sprintf("  Saved: %s (%.1f MB) in %.1fs\n", outfile, file_size / 1e6, save_time))
  cat(sprintf("  Rows: %s, Cols: %d, Cities: %d\n",
              format(nrow(dt), big.mark = ","), ncol(dt),
              length(unique(dt$city))))
  cat(sprintf("  Death mean: %.1f, range: %d-%d\n",
              mean(dt$death), min(dt$death), max(dt$death)))
  cat(sprintf("  Temp mean: %.1fC, range: %.1f-%.1fC\n",
              mean(dt$temp), min(dt$temp), max(dt$temp)))

  # -- Save truth specification (JSON) -- only once per run ---
  # Build the truth spec with ALL true parameters
  truth_spec <- list(
    dgp_description = paste(
      "Poisson GLM with cross-basis (DLNM) for temperature-lag-mortality.",
      "Model: death ~ Poisson(mu), log(mu) = intercept + cb %*% beta_true + ns(time,df=7) %*% gamma + dow_effects.",
      "Cross-basis: ns(df=5) x ns(df=4), lag 0-21 (C2 configuration).",
      "Between-city heterogeneity: beta_city = beta_true + N(0, tau_beta^2), intercept_city = intercept_true + N(0, tau_intercept^2)."
    ),
    model_formula = "death ~ crossbasis(temp, lag=c(0,21), argvar=list(fun='ns',df=5), arglag=list(fun='ns',df=4)) + ns(time, df=7) + dow",
    family = "quasipoisson (data generated as Poisson)",
    crossbasis_config = list(
      argvar = cb_argvar,
      arglag = cb_arglag,
      lag = cb_lag,
      n_var_basis = n_var_basis,
      n_lag_basis = n_lag_basis,
      n_cb_coefs = n_cb_coefs
    ),
    true_coefficients = list(
      intercept = intercept_true,
      beta_true = as.list(beta_true),
      beta_true_names = names(beta_true),
      gamma_true = as.list(gamma_true),
      gamma_true_names = names(gamma_true),
      dow_effects = as.list(dow_effects_true),
      dow_effects_names = names(dow_effects_true)
    ),
    between_city_heterogeneity = list(
      tau_beta = tau_beta,
      tau_intercept = tau_intercept,
      description = "beta_city = beta_true + N(0, tau_beta^2*I), intercept_city = intercept_true + N(0, tau_intercept^2)"
    ),
    temperature_generation = list(
      mean_global = temp_mean_global,
      mean_sd = temp_mean_sd,
      amplitude_mean = temp_amplitude_mean,
      amplitude_sd = temp_amplitude_sd,
      noise_sd = temp_noise_sd,
      formula = "temp = mean_temp_city + amplitude_city * sin(2*pi*(doy-80)/365) + N(0, noise_sd^2)"
    ),
    time_series = list(
      n_days = n_days,
      start_date = as.character(start_date),
      end_date = as.character(start_date + n_days - 1)
    ),
    scale = scale_name,
    n_cities = n_cities,
    n_rows = nrow(dt),
    random_seed_base = 42,
    city_seed_formula = "city_seed = 1000 + city_index",
    city_specific_params = lapply(city_params, function(p) {
      list(city = p$city, intercept = p$intercept,
           beta = as.list(p$beta), mean_temp = p$mean_temp,
           amplitude = p$amplitude)
    }),
    beta_design_rationale = paste(
      "Coefficients ordered as v1.l1, v1.l2, ..., v5.l4.",
      "v1 (cold edge): strong effect decaying with lag.",
      "v2 (cool): moderate effect.",
      "v3 (moderate/optimum): near-zero or slightly protective.",
      "v4 (warm): moderate positive effect.",
      "v5 (hot edge): strongest immediate effect, fast decay.",
      "This creates a U-shaped cumulative RR curve (min near 20C, high at extremes)."
    )
  )

  spec_file <- file.path(data_dir, "truth_spec.json")
  writeLines(toJSON(truth_spec, pretty = TRUE, auto_unbox = TRUE), spec_file)
  cat(sprintf("  Truth spec saved: %s\n", spec_file))

  # -- Verification: recover the true surface on a small subset ----------------
  if (scale_name == "10mb") {
    cat("\n=== VERIFICATION: Recovering True Surface ===\n\n")

    # Fit DLNM on the generated data
    t_verify <- proc.time()

    # Build crossbasis with group (same as benchmark pipeline)
    cb_verify <- crossbasis(dt$temp, lag = cb_lag,
                            argvar = cb_argvar, arglag = cb_arglag,
                            group = dt$city)

    # Fit GLM
    model_verify <- glm(death ~ cb_verify + ns(time, df = n_time_basis) + dow,
                        data = dt, family = quasipoisson())

    cat(sprintf("  GLM converged: %s (iterations: %d)\n",
                model_verify$converged, model_verify$iter))

    # Extract estimated overall cumulative RR via crosspred + compare with truth
    # Use a grid of temperatures for comparison
    temp_range <- range(dt$temp)
    pred_temps <- seq(temp_range[1], temp_range[2], length.out = 50)
    cen_temp <- 20  # Center near the optimum

    # Estimated RR
    cp <- crosspred(cb_verify, model_verify, at = pred_temps, cen = cen_temp)
    estimated_rr <- cp$allRRfit

    # True RR: compute from beta_true using the same basis functions
    # Build variable basis at prediction temps
    pred_var_basis <- do.call(onebasis, c(list(x = pred_temps), cb_argvar))
    cen_var_basis <- do.call(onebasis, c(list(x = cen_temp), cb_argvar))
    pred_var_cen <- scale(pred_var_basis, center = cen_var_basis, scale = FALSE)

    # Build lag basis
    lag_seq <- seqlag(cb_lag)
    pred_lag_basis <- do.call(onebasis, c(list(x = lag_seq), cb_arglag))

    # Overall cumulative effect: sum over all lags
    # Reduction matrix M: for each var basis col, multiply by sum of lag basis
    lag_sum <- colSums(pred_lag_basis)  # sum basislag over all lags
    # Overall reduced coef: for each var basis v, reduced_coef_v = sum_l(beta[v,l] * sum_j(basislag[j,l]))
    # where j runs over lags
    reduced_coef <- numeric(n_var_basis)
    for (v in seq_len(n_var_basis)) {
      for (l in seq_len(n_lag_basis)) {
        reduced_coef[v] <- reduced_coef[v] + beta_true[(v - 1) * n_lag_basis + l] * lag_sum[l]
      }
    }
    true_log_rr <- as.numeric(unclass(pred_var_cen) %*% reduced_coef)
    true_rr <- exp(true_log_rr)

    # Compute correlation
    corr <- cor(estimated_rr, true_rr)
    cat(sprintf("  Correlation (estimated vs true cumulative RR): %.4f\n", corr))
    cat(sprintf("  Verification %s (threshold: 0.9)\n",
                if (corr > 0.9) "PASSED" else "FAILED"))

    # Print a few comparison points
    cat("\n  Temperature | True RR | Estimated RR | Diff\n")
    cat("  -----------|---------|--------------|-----\n")
    for (idx in c(1, 10, 20, 25, 30, 40, 50)) {
      if (idx <= length(pred_temps)) {
        cat(sprintf("  %10.1fC | %7.4f | %12.4f | %+.4f\n",
                    pred_temps[idx], true_rr[idx], estimated_rr[idx],
                    estimated_rr[idx] - true_rr[idx]))
      }
    }

    verify_time <- (proc.time() - t_verify)[3]
    cat(sprintf("\n  Verification time: %.1fs\n", verify_time))
  }

  rm(dt, all_data, city_params)
  gc(verbose = FALSE)
  cat("\n")
}

cat("=== Truth Data Generation Complete ===\n")
