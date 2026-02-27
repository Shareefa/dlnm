# Environment

Environment variables, external dependencies, and setup notes.

**What belongs here:** Required env vars, external API keys/services, dependency quirks, platform-specific notes.
**What does NOT belong here:** Service ports/commands (use `.factory/services.yaml`).

---

## System

- **OS:** macOS ARM (aarch64-apple-darwin)
- **CPU:** 12 cores (Apple Silicon)
- **RAM:** 32 GB
- **R:** 4.5.2
- **Rust:** 1.87.0 (cargo 1.87.0)
- **Python:** 3.14.2

## R Package Dependencies

Core (installed): bench, data.table, splines, mgcv, tsModel, arrow, testthat, pkgload
Mission-specific (installed via init): rextendr, mixmeta

## Python Dependencies

Installed via init: matplotlib, numpy, pandas

## Benchmark Data

Pre-generated in `benchmarks/data/`:
- `scale_10mb.rds` — 122,736 rows (24 cities)
- `scale_100mb.rds` — 1,201,790 rows (235 cities)
- `scale_1gb.rds` — 12,017,900 rows (2,350 cities)
- `scale_10gb/` — 47 parquet chunks (23,500 cities, ~120M rows)

## ARM-Specific Notes

- Rust SIMD should target ARM NEON (`std::arch::aarch64`), NOT x86 AVX2
- Auto-vectorization with `-C target-cpu=native` is often sufficient on Apple Silicon
- R's BLAS uses Apple Accelerate framework (already optimized for ARM)
