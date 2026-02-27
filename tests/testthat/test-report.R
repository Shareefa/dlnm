###############################################################################
# test-report.R -- Tests for the final comprehensive markdown report
#
# Verifies: report file existence, all 7 required sections, formatted tables,
# PNG plot references, CSV-sourced numeric values, reproducible pipeline docs.
###############################################################################

library(testthat)
library(data.table)

# -- Helper: resolve paths relative to project root --------------------------
project_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
report_dir <- file.path(project_root, "benchmarks", "results", "report")
results_dir <- file.path(project_root, "benchmarks", "results")
report_file <- file.path(report_dir, "report.md")

# ============================================================================
# REPORT FILE EXISTENCE
# ============================================================================

test_that("report.md exists in benchmarks/results/report/", {
  expect_true(file.exists(report_file),
              info = "benchmarks/results/report/report.md does not exist")
})

# ============================================================================
# REPORT SECTION HEADERS
# ============================================================================

test_that("report has all 7 required sections", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- readLines(report_file, warn = FALSE)
  full_text <- paste(content, collapse = "\n")

  # Check for section headers (case-insensitive)
  expect_true(grepl("Introduction", full_text, ignore.case = TRUE),
              info = "Missing Introduction section")
  expect_true(grepl("Optimization Methodology", full_text, ignore.case = TRUE),
              info = "Missing Optimization Methodology section")
  expect_true(grepl("Benchmark Results", full_text, ignore.case = TRUE),
              info = "Missing Benchmark Results section")
  expect_true(grepl("Memory.*Scaling", full_text, ignore.case = TRUE),
              info = "Missing Memory and Scaling Analysis section")
  expect_true(grepl("Mixmeta Comparison", full_text, ignore.case = TRUE),
              info = "Missing Mixmeta Comparison section")
  expect_true(grepl("Convergence Analysis", full_text, ignore.case = TRUE),
              info = "Missing Convergence Analysis section")
  expect_true(grepl("Conclusions", full_text, ignore.case = TRUE),
              info = "Missing Conclusions section")
})

# ============================================================================
# FORMATTED TABLES
# ============================================================================

test_that("report includes formatted pipe-separated tables", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- readLines(report_file, warn = FALSE)

  # Count lines containing pipe characters (table rows)
  pipe_lines <- sum(grepl("\\|", content))
  expect_true(pipe_lines > 5,
              info = paste("Expected more than 5 pipe-separated table rows, found",
                           pipe_lines))
})

test_that("report includes at least one timing comparison table", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  # Should contain timing data references like seconds or speedup
  expect_true(grepl("speedup|time.*sec|median.*time", content, ignore.case = TRUE),
              info = "No timing comparison data found in report")
})

# ============================================================================
# PNG PLOT REFERENCES
# ============================================================================

test_that("report embeds PNG plots with markdown image syntax", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- readLines(report_file, warn = FALSE)

  # Count PNG references using ![Caption](filename.png) syntax
  png_refs <- sum(grepl("\\.png", content))
  expect_true(png_refs >= 4,
              info = paste("Expected at least 4 PNG references, found", png_refs))

  # Verify actual markdown image syntax
  img_refs <- sum(grepl("!\\[.*\\]\\(.*\\.png\\)", content))
  expect_true(img_refs >= 4,
              info = paste("Expected at least 4 ![Caption](file.png) refs, found",
                           img_refs))
})

test_that("referenced PNG files exist in report directory", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- readLines(report_file, warn = FALSE)

  # Extract PNG filenames from ![...](filename.png) patterns
  img_matches <- regmatches(content,
                            gregexpr("!\\[.*?\\]\\(([^)]+\\.png)\\)", content))
  img_matches <- unlist(img_matches)

  # Extract just the filenames
  png_files <- gsub("!\\[.*?\\]\\(([^)]+\\.png)\\)", "\\1", img_matches)
  png_files <- unique(png_files)

  expect_true(length(png_files) >= 4,
              info = "Expected at least 4 unique PNG files referenced")

  # Check that each referenced PNG exists

  for (pf in png_files) {
    full_path <- file.path(report_dir, pf)
    expect_true(file.exists(full_path),
                info = paste("Referenced PNG not found:", pf))
  }
})

# ============================================================================
# CONVERGENCE ANALYSIS CONTENT
# ============================================================================

test_that("report includes convergence analysis with quantitative values", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  # Check for quantitative convergence data
  expect_true(grepl("bias", content, ignore.case = TRUE),
              info = "No bias values in convergence section")
  expect_true(grepl("MSE|mse", content),
              info = "No MSE values in convergence section")
  expect_true(grepl("coverage", content, ignore.case = TRUE),
              info = "No CI coverage values in convergence section")
  expect_true(grepl("RMSE|rmse", content, ignore.case = TRUE),
              info = "No RMSE values in convergence section")
})

# ============================================================================
# MEMORY AND SCALING ANALYSIS
# ============================================================================

test_that("report includes memory and scaling analysis", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  expect_true(grepl("memory", content, ignore.case = TRUE),
              info = "No memory analysis in report")
  expect_true(grepl("scaling", content, ignore.case = TRUE),
              info = "No scaling analysis in report")
})

# ============================================================================
# CSV-SOURCED NUMERIC VALUES
# ============================================================================

test_that("at least one speedup value matches speedup_comparison.csv", {
  skip_if(!file.exists(report_file), "report.md not generated yet")

  speedup_file <- file.path(results_dir, "speedup_comparison.csv")
  skip_if(!file.exists(speedup_file), "speedup_comparison.csv not available")

  content <- paste(readLines(report_file, warn = FALSE), collapse = "\n")
  df <- fread(speedup_file)

  # Extract speedup factors and check if any appear in the report
  speedups <- df$speedup_factor
  found <- FALSE
  for (sp in speedups) {
    # Check for the value with at least 1 decimal place
    sp_str <- sprintf("%.2f", sp)
    if (grepl(sp_str, content, fixed = TRUE)) {
      found <- TRUE
      break
    }
    # Also try 1 decimal place
    sp_str1 <- sprintf("%.1f", sp)
    if (grepl(sp_str1, content, fixed = TRUE)) {
      found <- TRUE
      break
    }
  }
  expect_true(found,
              info = "No speedup values from CSV found in report")
})

test_that("at least one convergence metric matches convergence_metrics.csv", {
  skip_if(!file.exists(report_file), "report.md not generated yet")

  metrics_file <- file.path(results_dir, "convergence_metrics.csv")
  skip_if(!file.exists(metrics_file), "convergence_metrics.csv not available")

  content <- paste(readLines(report_file, warn = FALSE), collapse = "\n")
  df <- fread(metrics_file)

  # Check for RMSE or CI coverage values
  found <- FALSE
  for (i in seq_len(nrow(df))) {
    # Check RMSE (4 decimal places)
    rmse_str <- sprintf("%.4f", df$rmse[i])
    if (grepl(rmse_str, content, fixed = TRUE)) {
      found <- TRUE
      break
    }
    # Check CI coverage (2 decimal places as %)
    cov_str <- sprintf("%.2f", df$ci_coverage[i])
    if (grepl(cov_str, content, fixed = TRUE)) {
      found <- TRUE
      break
    }
  }
  expect_true(found,
              info = "No convergence metrics from CSV found in report")
})

# ============================================================================
# REPRODUCIBLE PIPELINE
# ============================================================================

test_that("reproducible pipeline documentation exists", {
  skip_if(!file.exists(report_file), "report.md not generated yet")
  content <- paste(readLines(report_file, warn = FALSE), collapse = "\n")

  # Check for pipeline/reproduce section in report
  has_pipeline_section <- grepl("reproducib|pipeline|reproduce", content,
                                ignore.case = TRUE)
  expect_true(has_pipeline_section,
              info = "No reproducible pipeline documentation found in report")

  # Check for actual commands (should mention Rscript or cargo or bash)
  has_commands <- grepl("Rscript|cargo build|bash", content)
  expect_true(has_commands,
              info = "Pipeline documentation should include executable commands")
})
