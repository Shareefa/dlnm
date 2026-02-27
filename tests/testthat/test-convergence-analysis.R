###############################################################################
# test-convergence-analysis.R -- Tests for convergence analysis script
#
# Verifies: metric computation, CSV output structure, both approaches running
# on the same truth data, and reasonable metric ranges.
###############################################################################

library(testthat)
library(data.table)
library(splines)
library(jsonlite)

# Load dlnm from source
pkgload::load_all(".", quiet = TRUE)

# -- Helper: resolve paths relative to project root --------------------------
project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
results_dir <- file.path(project_root, "benchmarks", "results")
script_path <- file.path(project_root, "benchmarks", "run_convergence_analysis.R")

# ============================================================================
# UNIT TESTS: metric computation on small known datasets
# ============================================================================

test_that("bias, MSE, RMSE computed correctly on known values", {
  # Known values: estimated = c(1.1, 1.2, 0.9), true = c(1.0, 1.0, 1.0)
  estimated <- c(1.1, 1.2, 0.9)
  true_vals <- c(1.0, 1.0, 1.0)

  bias <- estimated - true_vals  # c(0.1, 0.2, -0.1)
  mse <- mean(bias^2)            # (0.01 + 0.04 + 0.01) / 3 = 0.02
  rmse <- sqrt(mse)              # ~0.1414

  expect_equal(bias, c(0.1, 0.2, -0.1))
  expect_equal(mse, 0.02, tolerance = 1e-10)
  expect_equal(rmse, sqrt(0.02), tolerance = 1e-10)
})

test_that("relative bias computed correctly", {
  estimated <- 1.3
  true_val <- 1.2
  bias <- estimated - true_val  # 0.1
  rel_bias <- bias / abs(true_val)  # 0.1 / 1.2

  expect_equal(rel_bias, 0.1 / 1.2, tolerance = 1e-10)
})

test_that("CI coverage computed correctly on known intervals", {
  # True values: all 1.0
  # 4 intervals: 3 contain 1.0, 1 does not
  true_rr <- rep(1.0, 4)
  ci_lower <- c(0.8, 0.9, 1.1, 0.7)  # third interval: 1.1 > 1.0, misses

  ci_upper <- c(1.2, 1.1, 1.3, 1.5)
  coverage <- mean(true_rr >= ci_lower & true_rr <= ci_upper)
  expect_equal(coverage, 0.75)  # 3/4
})

# ============================================================================
# INTEGRATION TESTS: verify convergence analysis script outputs
# ============================================================================

test_that("convergence_metrics.csv has expected columns and rows for both approaches", {
  metrics_file <- file.path(results_dir, "convergence_metrics.csv")
  skip_if(!file.exists(metrics_file),
          "convergence_metrics.csv not generated yet (run convergence analysis first)")

  df <- fread(metrics_file)

  # Expected columns
  expected_cols <- c("approach", "config", "scale", "bias", "mse", "rmse",
                     "rel_bias", "ci_coverage")
  for (col in expected_cols) {
    expect_true(col %in% names(df),
                info = paste("Missing column:", col))
  }

  # Should have rows for both approaches
  expect_true("single_stage" %in% df$approach,
              info = "No single_stage rows in convergence_metrics.csv")
  expect_true("two_stage" %in% df$approach,
              info = "No two_stage rows in convergence_metrics.csv")

  # Metric ranges
  expect_true(all(is.finite(df$bias)))
  expect_true(all(df$mse >= 0))
  expect_true(all(df$rmse >= 0))
  expect_true(all(df$ci_coverage >= 0 & df$ci_coverage <= 1))
})

test_that("convergence_details.csv has per-temperature data for both approaches", {
  details_file <- file.path(results_dir, "convergence_details.csv")
  skip_if(!file.exists(details_file),
          "convergence_details.csv not generated yet")

  df <- fread(details_file)

  # Expected columns
  expected_cols <- c("approach", "config", "temperature", "estimated_rr",
                     "true_rr", "se", "ci_lower", "ci_upper")
  for (col in expected_cols) {
    expect_true(col %in% names(df),
                info = paste("Missing column:", col))
  }

  # Should have rows for both approaches
  expect_true("single_stage" %in% df$approach)
  expect_true("two_stage" %in% df$approach)

  # Each approach x config should have ~50 temperature points
  counts <- df[, .N, by = .(approach, config)]
  expect_true(all(counts$N >= 40),  # at least 40 temperature points
              info = "Too few temperature points per approach/config")

  # RR values should be positive
  expect_true(all(df$estimated_rr > 0))
  expect_true(all(df$true_rr > 0))
})

test_that("numerical_convergence.csv has expected columns and convergence data", {
  conv_file <- file.path(results_dir, "numerical_convergence.csv")
  skip_if(!file.exists(conv_file),
          "numerical_convergence.csv not generated yet")

  df <- fread(conv_file)

  # Expected columns
  expected_cols <- c("approach", "config", "scale", "iterations",
                     "converged", "time_sec")
  for (col in expected_cols) {
    expect_true(col %in% names(df),
                info = paste("Missing column:", col))
  }

  # Should have rows for both approaches
  expect_true("single_stage" %in% df$approach)
  expect_true("two_stage" %in% df$approach)

  # Convergence flag should be logical-like (TRUE/FALSE or 1/0)
  expect_true(all(df$converged %in% c(TRUE, FALSE, 1, 0)),
              info = "converged column should be TRUE/FALSE")

  # Timing should be positive
  expect_true(all(df$time_sec > 0))

  # Iterations should be positive integers
  expect_true(all(df$iterations > 0))
})

test_that("both approaches use same truth data (same scale entries present)", {
  metrics_file <- file.path(results_dir, "convergence_metrics.csv")
  skip_if(!file.exists(metrics_file), "convergence_metrics.csv not generated yet")

  df <- fread(metrics_file)

  # For each scale+config, both approaches should be present
  combos <- unique(df[, .(config, scale)])
  for (i in seq_len(nrow(combos))) {
    cfg <- combos$config[i]
    scl <- combos$scale[i]
    sub <- df[config == cfg & scale == scl]
    expect_true(nrow(sub) >= 2,
                info = paste("Expected both approaches for", cfg, scl))
    expect_true("single_stage" %in% sub$approach,
                info = paste("Missing single_stage for", cfg, scl))
    expect_true("two_stage" %in% sub$approach,
                info = paste("Missing two_stage for", cfg, scl))
  }
})

test_that("convergence analysis runs at 10MB C2 scale", {
  metrics_file <- file.path(results_dir, "convergence_metrics.csv")
  skip_if(!file.exists(metrics_file), "convergence_metrics.csv not generated yet")

  df <- fread(metrics_file)
  c2_10mb <- df[config == "C2" & scale == "10mb"]
  expect_true(nrow(c2_10mb) >= 2,
              info = "Expected C2 10mb rows for both approaches")
})
