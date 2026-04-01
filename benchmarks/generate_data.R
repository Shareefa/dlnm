#!/usr/bin/env Rscript
###############################################################################
# generate_data.R -- Generate scaled DLNM benchmark datasets
#
# Creates multi-city panel datasets by replicating chicagoNMMAPS with
# realistic perturbations. Scales: 10MB, 100MB, 1GB, 10GB.
#
# Usage: Rscript benchmarks/generate_data.R [scales]
#   scales: comma-separated list, e.g. "10mb,100mb" (default: all)
#
# Output:
#   benchmarks/data/scale_10mb.rds
#   benchmarks/data/scale_100mb.rds
#   benchmarks/data/scale_1gb.rds
#   benchmarks/data/scale_10gb/chunk_001.parquet ... (requires arrow)
###############################################################################

# -- Dependencies -------------------------------------------------------------
required_pkgs <- c("data.table")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
library(data.table)

# -- Parse CLI args -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
all_scales <- c("10mb", "100mb", "1gb", "10gb")
if (length(args) > 0) {
  requested <- tolower(unlist(strsplit(args[1], ",")))
  invalid <- setdiff(requested, all_scales)
  if (length(invalid) > 0) stop("Unknown scales: ", paste(invalid, collapse = ", "))
  scales_to_generate <- requested
} else {
  scales_to_generate <- all_scales
}

# -- Load source data ---------------------------------------------------------
load("data/chicagoNMMAPS.rda")
orig <- as.data.table(chicagoNMMAPS)
n_orig <- nrow(orig)  # 5114

cat("Source: chicagoNMMAPS -", n_orig, "rows x", ncol(orig), "cols\n")
cat("Scales to generate:", paste(scales_to_generate, collapse = ", "), "\n\n")

# -- Compute correlation structure of continuous vars -------------------------
env_cols <- c("temp", "dptp", "rhum", "pm10", "o3")
env_mat <- as.matrix(orig[, ..env_cols])
# Use pairwise complete obs for correlation (pm10 has NAs)
cor_mat <- cor(env_mat, use = "pairwise.complete.obs")
# Cholesky decomposition for correlated noise generation
chol_L <- chol(cor_mat)
# Column standard deviations for scaling noise
env_sds <- apply(env_mat, 2, sd, na.rm = TRUE)

# -- Scale definitions --------------------------------------------------------
# chicagoNMMAPS: 5114 rows x 15 cols (with city column), ~10.3 MB for 23 cities
# Using fixed city counts from the plan to hit target in-memory sizes
scale_config <- list(
  "10mb"  = list(n_cities = 24,    method = "memory"),
  "100mb" = list(n_cities = 235,   method = "memory"),
  "1gb"   = list(n_cities = 2350,  method = "memory"),
  "10gb"  = list(n_cities = 23500, method = "chunked", chunk_cities = 500)
)

# -- Per-city perturbation function -------------------------------------------
perturb_city <- function(dt, city_id, chol_L, env_sds) {
  # Deep copy to avoid modifying original
  city <- copy(dt)
  n <- nrow(city)

  # 1. Add correlated noise to environmental variables
  #    Generate uncorrelated standard normal, then apply Cholesky factor
  noise_raw <- matrix(rnorm(n * 5), nrow = n, ncol = 5)
  noise_corr <- noise_raw %*% chol_L
  # Scale noise to ~20% of original SD
  noise_scaled <- sweep(noise_corr, 2, env_sds * 0.2, "*")

  city[, temp := temp + noise_scaled[, 1]]
  city[, dptp := dptp + noise_scaled[, 2]]
  city[, rhum := pmin(100, pmax(0, rhum + noise_scaled[, 3]))]
  city[, pm10 := pmax(0, pm10 + noise_scaled[, 4])]
  city[, o3   := pmax(0, o3 + noise_scaled[, 5])]

  # 2. Shift temperature mean to simulate different climates
  temp_shift <- rnorm(1, mean = 0, sd = 5)
  city[, temp := temp + temp_shift]

  # 3. Scale pollution levels to simulate different pollution environments
  poll_scale <- runif(1, min = 0.5, max = 2.0)
  city[, pm10 := pm10 * poll_scale]
  city[, o3   := o3 * poll_scale]

  # 4. Regenerate mortality counts from perturbed Poisson process
  death_scale <- runif(1, min = 0.5, max = 3.0)
  base_rate <- dt$death * death_scale
  city[, death := rpois(n, lambda = pmax(1, base_rate))]
  city[, cvd   := rpois(n, lambda = pmax(1, dt$cvd * death_scale))]
  city[, resp  := rpois(n, lambda = pmax(1, dt$resp * death_scale))]

  # 5. Add city identifier
  city[, city := city_id]

  # 6. Preserve date/time structure (doy, dow, month already present)
  return(city)
}

# -- Generate a batch of cities -----------------------------------------------
generate_batch <- function(n_cities, start_id, orig, chol_L, env_sds) {
  set.seed(42 + start_id)  # Reproducible per batch
  cities <- vector("list", n_cities)
  for (i in seq_len(n_cities)) {
    city_id <- sprintf("city_%05d", start_id + i - 1)
    cities[[i]] <- perturb_city(orig, city_id, chol_L, env_sds)
  }
  rbindlist(cities)
}

# -- Output directory ---------------------------------------------------------
data_dir <- "benchmarks/data"
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -- Generate each requested scale --------------------------------------------
for (scale_name in scales_to_generate) {
  cfg <- scale_config[[scale_name]]
  n_cities <- cfg$n_cities
  total_rows <- n_cities * n_orig

  cat(sprintf("=== %s: %d cities, ~%s rows ===\n",
              toupper(scale_name), n_cities, format(total_rows, big.mark = ",")))

  if (cfg$method == "memory") {
    # --- In-memory generation, save as .rds ---
    t0 <- proc.time()
    dt <- generate_batch(n_cities, 1, orig, chol_L, env_sds)
    dt[, city := as.factor(city)]
    gen_time <- (proc.time() - t0)[3]

    outfile <- file.path(data_dir, paste0("scale_", scale_name, ".rds"))
    t1 <- proc.time()
    saveRDS(dt, outfile)
    save_time <- (proc.time() - t1)[3]

    file_size <- file.size(outfile)
    cat(sprintf("  Generated in %.1fs, saved in %.1fs\n", gen_time, save_time))
    cat(sprintf("  File: %s (%.1f MB)\n", outfile,  file_size / 1e6))
    cat(sprintf("  In-memory: %s\n", format(object.size(dt), units = "auto")))

    # Quick validation
    cat(sprintf("  Rows: %s, Cols: %d, Cities: %d\n",
                format(nrow(dt), big.mark = ","), ncol(dt), length(unique(dt$city))))
    cat(sprintf("  temp mean=%.2f (orig=%.2f), death mean=%.1f (orig=%.1f)\n",
                mean(dt$temp, na.rm = TRUE), mean(orig$temp, na.rm = TRUE),
                mean(dt$death), mean(orig$death)))
    rm(dt); gc(verbose = FALSE)

  } else {
    # --- Chunked generation for 10GB, save as parquet ---
    if (!requireNamespace("arrow", quietly = TRUE)) {
      message("  Skipping 10gb: 'arrow' package not installed.")
      message("  Install with: install.packages('arrow')")
      next
    }
    library(arrow)

    chunk_dir <- file.path(data_dir, paste0("scale_", scale_name))
    dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

    chunk_size <- cfg$chunk_cities
    n_chunks <- ceiling(n_cities / chunk_size)
    total_size <- 0
    t0 <- proc.time()

    for (chunk_i in seq_len(n_chunks)) {
      start_id <- (chunk_i - 1) * chunk_size + 1
      end_id <- min(chunk_i * chunk_size, n_cities)
      n_this <- end_id - start_id + 1

      cat(sprintf("\r  Chunk %d/%d (cities %d-%d)...",
                  chunk_i, n_chunks, start_id, end_id))

      dt_chunk <- generate_batch(n_this, start_id, orig, chol_L, env_sds)
      dt_chunk[, city := as.factor(city)]

      chunk_file <- file.path(chunk_dir, sprintf("chunk_%03d.parquet", chunk_i))
      write_parquet(dt_chunk, chunk_file)
      total_size <- total_size + file.size(chunk_file)

      rm(dt_chunk); gc(verbose = FALSE)
    }

    gen_time <- (proc.time() - t0)[3]
    cat(sprintf("\n  Generated %d chunks in %.1fs\n", n_chunks, gen_time))
    cat(sprintf("  Total size: %.1f GB\n", total_size / 1e9))
    cat(sprintf("  Total rows: ~%s\n", format(total_rows, big.mark = ",")))
  }
  cat("\n")
}

# -- Verification summary -----------------------------------------------------
cat("=== Verification Summary ===\n")
cat("Original chicagoNMMAPS correlation structure (temp, dptp, rhum, pm10, o3):\n")
print(round(cor_mat, 3))

# Check the smallest generated .rds file for correlation preservation
rds_files <- list.files(data_dir, pattern = "^scale_.*\\.rds$", full.names = TRUE)
if (length(rds_files) > 0) {
  check_file <- rds_files[which.min(file.size(rds_files))]
  dt_check <- readRDS(check_file)
  gen_cor <- cor(as.matrix(dt_check[, ..env_cols]), use = "pairwise.complete.obs")
  cat(sprintf("\nGenerated %s correlation structure:\n", basename(check_file)))
  print(round(gen_cor, 3))
  cat("\nMax absolute correlation difference:",
      round(max(abs(cor_mat - gen_cor)), 4), "\n")
  rm(dt_check); gc(verbose = FALSE)
}

cat("\nDone.\n")
