###############################################################################
# test-truth-data.R -- Tests for truth data generation script
#
# Verifies: data structure, truth_spec.json, true surface recoverability,
# and DGP consistency.
###############################################################################

library(testthat)
library(data.table)
library(splines)
library(jsonlite)

# Load dlnm from source
pkgload::load_all(".", quiet = TRUE)

# -- Helper: generate truth data at minimal scale for testing -----------------
# We run the DGP inline (matching generate_truth_data.R) rather than calling
# the script, so tests are self-contained and fast.
generate_small_truth <- function(n_cities = 3, n_days = 365) {
  set.seed(42)

  intercept_true <- 3.5
  cb_argvar <- list(fun = "ns", df = 5)
  cb_arglag <- list(fun = "ns", df = 4)
  cb_lag <- c(0, 21)
  n_var_basis <- 5
  n_lag_basis <- 4

  beta_true <- c(
    0.08, 0.05, 0.03, 0.015,
    0.04, 0.025, 0.015, 0.008,
    -0.01, -0.005, -0.003, 0.001,
    0.05, 0.03, 0.02, 0.01,
    0.10, 0.06, 0.03, 0.015
  )

  gamma_true <- c(0.05, -0.03, 0.04, -0.02, 0.03, -0.01, 0.02)
  dow_effects_true <- c(-0.020, -0.030, -0.025, -0.020, -0.010, 0.010)
  tau_beta <- 0.005
  tau_intercept <- 0.10

  start_date <- as.Date("1987-01-01")
  n_time_basis <- 7

  all_data <- list()
  for (i in seq_len(n_cities)) {
    set.seed(1000 + i)
    city_id <- sprintf("city_%05d", i)
    intercept_city <- intercept_true + rnorm(1, 0, tau_intercept)
    beta_city <- beta_true + rnorm(length(beta_true), 0, tau_beta)

    dates <- start_date + 0:(n_days - 1)
    doy <- as.integer(format(dates, "%j"))
    dow <- factor(weekdays(dates),
                  levels = c("Sunday", "Monday", "Tuesday", "Wednesday",
                             "Thursday", "Friday", "Saturday"))
    time_idx <- seq_len(n_days)

    mean_temp <- rnorm(1, 15, 5)
    amp <- rnorm(1, 12, 2)
    temp <- mean_temp + amp * sin(2 * pi * (doy - 80) / 365) + rnorm(n_days, 0, 3)

    cb <- crossbasis(temp, lag = cb_lag, argvar = cb_argvar, arglag = cb_arglag)
    time_basis <- ns(time_idx, df = n_time_basis)

    cb_effect <- as.numeric(unclass(cb) %*% beta_city)
    cb_effect[is.na(cb_effect)] <- 0

    log_mu <- intercept_city + cb_effect + as.numeric(time_basis %*% gamma_true)

    dow_num <- as.integer(dow)
    for (d in 2:7) {
      log_mu[dow_num == d] <- log_mu[dow_num == d] + dow_effects_true[d - 1]
    }

    mu <- exp(pmin(log_mu, log(1e6)))
    death <- rpois(n_days, lambda = mu)

    all_data[[i]] <- data.table(
      city = city_id,
      date = dates,
      time = time_idx,
      year = as.integer(format(dates, "%Y")),
      doy = doy,
      dow = dow,
      temp = round(temp, 2),
      death = death
    )
  }

  dt <- rbindlist(all_data)
  dt[, city := as.factor(city)]

  list(
    data = dt,
    beta_true = beta_true,
    gamma_true = gamma_true,
    dow_effects_true = dow_effects_true,
    intercept_true = intercept_true,
    cb_argvar = cb_argvar,
    cb_arglag = cb_arglag,
    cb_lag = cb_lag
  )
}

# ============================================================================
# TESTS
# ============================================================================

test_that("truth data has expected structure: death, temp, time, dow, city columns", {
  truth <- generate_small_truth()
  dt <- truth$data

  # Required columns
  expected_cols <- c("death", "temp", "time", "dow", "city")
  expect_true(all(expected_cols %in% names(dt)))

  # Column types
  expect_true(is.integer(dt$death))
  expect_true(is.numeric(dt$temp))
  expect_true(is.integer(dt$time))
  expect_true(is.factor(dt$dow))
  expect_true(is.factor(dt$city))

  # No NAs in key columns
  expect_equal(sum(is.na(dt$death)), 0)
  expect_equal(sum(is.na(dt$temp)), 0)
  expect_equal(sum(is.na(dt$city)), 0)
})

test_that("truth data has expected dimensions", {
  truth <- generate_small_truth(n_cities = 3, n_days = 365)
  dt <- truth$data

  expect_equal(nrow(dt), 3 * 365)
  expect_equal(length(unique(dt$city)), 3)
})

test_that("death counts are positive integers", {
  truth <- generate_small_truth()
  dt <- truth$data

  expect_true(all(dt$death >= 0))
  expect_true(all(dt$death == as.integer(dt$death)))
  expect_true(mean(dt$death) > 10)  # reasonable death rate
  expect_true(mean(dt$death) < 200)
})

test_that("temperature has realistic seasonal pattern", {
  truth <- generate_small_truth(n_cities = 1, n_days = 365)
  dt <- truth$data

  # Temperature should vary seasonally
  expect_true(sd(dt$temp) > 3)     # non-trivial variation
  expect_true(max(dt$temp) - min(dt$temp) > 20)  # reasonable range
})

test_that("true beta coefficients have expected count", {
  truth <- generate_small_truth()
  expect_equal(length(truth$beta_true), 20)  # 5 var x 4 lag = 20
})

test_that("true surface is recoverable (correlation > 0.9) from pooled model", {
  truth <- generate_small_truth(n_cities = 10, n_days = 5114)
  dt <- truth$data

  # Fit DLNM on truth data
  cb <- crossbasis(dt$temp, lag = truth$cb_lag,
                   argvar = truth$cb_argvar, arglag = truth$cb_arglag,
                   group = dt$city)

  model <- glm(death ~ cb + ns(time, df = 7) + dow,
               data = dt, family = quasipoisson())

  expect_true(model$converged)

  # Get estimated overall cumulative RR
  temp_range <- range(dt$temp)
  pred_temps <- seq(temp_range[1], temp_range[2], length.out = 50)
  cen_temp <- 20

  cp <- crosspred(cb, model, at = pred_temps, cen = cen_temp)
  estimated_rr <- cp$allRRfit

  # Compute true overall cumulative RR
  pred_var_basis <- do.call(onebasis, c(list(x = pred_temps), truth$cb_argvar))
  cen_var_basis <- do.call(onebasis, c(list(x = cen_temp), truth$cb_argvar))
  pred_var_cen <- scale(pred_var_basis, center = cen_var_basis, scale = FALSE)

  lag_seq <- seqlag(truth$cb_lag)
  pred_lag_basis <- do.call(onebasis, c(list(x = lag_seq), truth$cb_arglag))
  lag_sum <- colSums(pred_lag_basis)

  n_var <- 5
  n_lag <- 4
  reduced_coef <- numeric(n_var)
  for (v in seq_len(n_var)) {
    for (l in seq_len(n_lag)) {
      reduced_coef[v] <- reduced_coef[v] + truth$beta_true[(v - 1) * n_lag + l] * lag_sum[l]
    }
  }
  true_rr <- exp(as.numeric(unclass(pred_var_cen) %*% reduced_coef))

  corr <- cor(estimated_rr, true_rr)
  expect_gt(corr, 0.9)
})

test_that("truth_spec.json exists and has required fields (if generated)", {
  # Locate from either tests/testthat/ or project root
  spec_file <- "benchmarks/data/truth_spec.json"
  if (!file.exists(spec_file)) spec_file <- "../../benchmarks/data/truth_spec.json"
  skip_if(!file.exists(spec_file), "truth_spec.json not generated yet")

  spec <- fromJSON(spec_file)

  # Check top-level fields
  expect_true("dgp_description" %in% names(spec))
  expect_true("model_formula" %in% names(spec))
  expect_true("crossbasis_config" %in% names(spec))
  expect_true("true_coefficients" %in% names(spec))
  expect_true("between_city_heterogeneity" %in% names(spec))
  expect_true("temperature_generation" %in% names(spec))

  # Check true coefficients
  coefs <- spec$true_coefficients
  expect_true("intercept" %in% names(coefs))
  expect_true("beta_true" %in% names(coefs))
  expect_true("gamma_true" %in% names(coefs))
  expect_true("dow_effects" %in% names(coefs))

  # Check beta_true has 20 coefficients
  beta <- unlist(coefs$beta_true)
  expect_equal(length(beta), 20)
})

test_that("truth_10mb.rds exists with expected properties (if generated)", {
  data_file <- "benchmarks/data/truth_10mb.rds"
  if (!file.exists(data_file)) data_file <- "../../benchmarks/data/truth_10mb.rds"
  skip_if(!file.exists(data_file), "truth_10mb.rds not generated yet")

  dt <- readRDS(data_file)

  # Expected columns
  expected_cols <- c("city", "date", "time", "year", "doy", "dow", "temp", "death")
  expect_true(all(expected_cols %in% names(dt)))

  # Dimensions
  expect_gt(nrow(dt), 100000)  # 24 cities x 5114 = 122,736
  expect_equal(ncol(dt), 8)

  # No NAs
  expect_equal(sum(is.na(dt$death)), 0)
  expect_equal(sum(is.na(dt$temp)), 0)

  # Multiple cities
  expect_gte(length(unique(dt$city)), 20)
})
