###############################################################################
# test-mixmeta-pipeline.R — Tests for the two-stage DLNM + mixmeta pipeline
###############################################################################

# Helper: resolve paths relative to project root (testthat CWD is tests/testthat/)
project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
data_path_10mb <- file.path(project_root, "benchmarks", "data", "scale_10mb.rds")
pipeline_script <- file.path(project_root, "benchmarks", "run_mixmeta_pipeline.R")
results_dir <- file.path(project_root, "benchmarks", "results")

test_that("mixmeta package loads and can perform basic meta-analysis", {
  skip_if_not_installed("mixmeta")
  library(mixmeta)

  # Create dummy two-outcome meta-analysis
  set.seed(42)
  n_studies <- 5
  coefs <- matrix(rnorm(n_studies * 2), ncol = 2)
  vcovs <- lapply(seq_len(n_studies), function(i) {
    S <- diag(runif(2, 0.01, 0.1))
    S
  })

  fit <- mixmeta(coefs ~ 1, S = vcovs, method = "reml")
  expect_s3_class(fit, "mixmeta")
  expect_length(coef(fit), 2)
  expect_true(!any(is.na(coef(fit))))
})

test_that("two-stage pipeline produces pooled estimates at 10MB", {
  skip_if_not_installed("mixmeta")
  skip_if(!file.exists(data_path_10mb),
          "10MB benchmark data not found")

  pkgload::load_all(project_root, quiet = TRUE)
  library(mixmeta)
  library(splines)

  dt <- readRDS(data_path_10mb)
  cities <- unique(dt$city)
  expect_true(length(cities) >= 2)

  # Use C2 config: ns(5) x ns(4), lag 0-21
  argvar <- list(fun = "ns", df = 5)
  arglag <- list(fun = "ns", df = 4)
  lag_range <- c(0, 21)

  # Stage 1: per-city DLNM
  city_coefs <- list()
  city_vcovs <- list()
  # Test on first 3 cities for speed
  test_cities <- cities[1:min(3, length(cities))]

  for (city_name in test_cities) {
    city_dt <- dt[dt$city == city_name, ]
    n_years <- length(unique(city_dt$year))
    time_df <- min(7 * n_years, 100)

    cb <- crossbasis(city_dt$temp, lag = lag_range,
                     argvar = argvar, arglag = arglag)
    model <- glm(death ~ cb + ns(time, df = time_df) + dow,
                 data = city_dt, family = quasipoisson())

    # crossreduce to get overall coefficients
    cr <- crossreduce(cb, model, type = "overall")
    city_coefs[[city_name]] <- coef(cr)
    city_vcovs[[city_name]] <- cr$vcov
  }

  # Stage 2: Pool with mixmeta
  coef_matrix <- do.call(rbind, city_coefs)
  vcov_list <- city_vcovs

  pooled <- mixmeta(coef_matrix ~ 1, S = vcov_list, method = "reml")
  expect_s3_class(pooled, "mixmeta")
  expect_true(!any(is.na(coef(pooled))))
  expect_equal(length(coef(pooled)), ncol(coef_matrix))
})

test_that("heterogeneity statistics are computed and in plausible ranges", {
  skip_if_not_installed("mixmeta")
  skip_if(!file.exists(data_path_10mb),
          "10MB benchmark data not found")

  pkgload::load_all(project_root, quiet = TRUE)
  library(mixmeta)
  library(splines)

  dt <- readRDS(data_path_10mb)
  cities <- unique(dt$city)

  argvar <- list(fun = "ns", df = 5)
  arglag <- list(fun = "ns", df = 4)
  lag_range <- c(0, 21)

  city_coefs <- list()
  city_vcovs <- list()
  test_cities <- cities[1:min(5, length(cities))]

  for (city_name in test_cities) {
    city_dt <- dt[dt$city == city_name, ]
    n_years <- length(unique(city_dt$year))
    time_df <- min(7 * n_years, 100)

    cb <- crossbasis(city_dt$temp, lag = lag_range,
                     argvar = argvar, arglag = arglag)
    model <- glm(death ~ cb + ns(time, df = time_df) + dow,
                 data = city_dt, family = quasipoisson())
    cr <- crossreduce(cb, model, type = "overall")
    city_coefs[[city_name]] <- coef(cr)
    city_vcovs[[city_name]] <- cr$vcov
  }

  coef_matrix <- do.call(rbind, city_coefs)
  vcov_list <- city_vcovs

  pooled <- mixmeta(coef_matrix ~ 1, S = vcov_list, method = "reml")
  summ <- summary(pooled)

  # Check that Psi (between-study variance) exists
  expect_true(!is.null(pooled$Psi))

  # Cochran Q test should be available
  q_stat <- qtest(pooled)
  expect_true(!is.null(q_stat))
  expect_true(is.numeric(q_stat$Q))
  expect_true(all(q_stat$Q >= 0))

  # Compute I-squared manually: I2 = max(0, (Q - df) / Q * 100)
  # For multivariate models, use per-outcome Q (skip first .all entry)
  per_q <- q_stat$Q[-1]
  per_df <- q_stat$df[-1]
  i_squared <- pmax(0, (per_q - per_df) / per_q * 100)
  expect_true(is.numeric(i_squared))
  expect_true(all(i_squared >= 0 & i_squared <= 100))
})

test_that("timing CSV has expected columns", {
  skip_if_not_installed("mixmeta")
  skip_if(!file.exists(data_path_10mb),
          "10MB benchmark data not found")

  # Run the pipeline script on 10MB from the project root
  # Use sh -c to cd first since system2 inherits testthat's working directory
  result <- system2("sh",
    args = c("-c", shQuote(paste0(
      "cd ", shQuote(project_root),
      " && Rscript ", shQuote(pipeline_script), " 10mb C2"))),
    stdout = TRUE, stderr = TRUE)
  exit_code <- attr(result, "status")
  if (is.null(exit_code)) exit_code <- 0
  expect_equal(exit_code, 0,
    info = paste("Pipeline failed:", paste(tail(result, 20), collapse = "\n")))

  # Check timing CSV
  timing_file <- file.path(results_dir, "mixmeta_timing.csv")
  expect_true(file.exists(timing_file))
  timing_df <- read.csv(timing_file)
  expected_cols <- c("scale", "config", "stage", "n_cities", "time_sec")
  for (col in expected_cols) {
    expect_true(col %in% names(timing_df),
      info = paste("Missing column:", col))
  }
  # Should have rows for 10mb
  expect_true(any(timing_df$scale == "10mb"))
})

test_that("estimates CSV has pooled coefficients", {
  # This test relies on the pipeline having been run by the previous test
  estimates_file <- file.path(results_dir, "mixmeta_estimates.csv")
  skip_if(!file.exists(estimates_file),
          "Run the pipeline first to create estimates CSV")

  est_df <- read.csv(estimates_file)
  expect_true(nrow(est_df) > 0)
  # Should have scale, config, and coefficient columns
  expect_true("scale" %in% names(est_df))
  expect_true("config" %in% names(est_df))
  # Should have pooled coefficient values and heterogeneity stats
  expect_true("pooled_coef" %in% names(est_df))
  expect_true("pooled_se" %in% names(est_df))
  expect_true("tau2" %in% names(est_df))
  # Coefficients should be finite numbers
  expect_true(all(is.finite(est_df$pooled_coef)))
})
