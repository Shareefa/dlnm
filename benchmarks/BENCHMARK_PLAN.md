# DLNM Benchmarking & Rust Optimization Plan

## Context

The dlnm R package fits Distributed Lag Non-Linear Models -- bi-dimensional models capturing both exposure-response and lag-response relationships via tensor products of basis functions. The core pipeline is: `crossbasis()` -> `glm()`/`gam()` -> `crosspred()` -> `crossreduce()`. We need to (1) benchmark this pipeline at 10MB, 100MB, 1GB, and 10GB data scales, and (2) identify specific functions to rewrite in Rust with SIMD for performance gains.

---

## Part 1: Dataset Generation

**File to create:** `benchmarks/generate_data.R`

### Approach: Multi-City Panel Replication
Scale by replicating `chicagoNMMAPS` (5,114 rows x 14 cols, ~85 bytes/row) across virtual "cities" with realistic perturbations:

| Scale | Cities | Rows | Strategy |
|-------|--------|------|----------|
| 10 MB | 24 | ~123K | In-memory, save as `.rds` |
| 100 MB | 235 | ~1.2M | In-memory, save as `.rds` |
| 1 GB | 2,350 | ~12M | In-memory (needs ~4GB RAM for generation), save as `.rds` |
| 10 GB | 23,500 | ~120M | Chunked generation in batches of 500 cities, save as partitioned `.parquet` files |

### Per-City Perturbation Logic
1. Compute correlation matrix of `(temp, dptp, rhum, pm10, o3)` from original data
2. For each city: add correlated noise via Cholesky decomposition of the correlation matrix
3. Shift temperature mean by `rnorm(1, 0, 5)` to simulate different climates
4. Scale pollution by `runif(1, 0.5, 2.0)` to simulate different pollution levels
5. Regenerate mortality counts: `rpois(n, original_death * runif(1, 0.5, 3.0))`
6. Add `city` column as factor (used by `crossbasis()` `group` argument)
7. Preserve original seasonality structure (`doy`, `dow`, `month`)

### Output Structure
```
benchmarks/
  generate_data.R
  data/                     # .gitignore this directory
    scale_10mb.rds
    scale_100mb.rds
    scale_1gb.rds
    scale_10gb/             # partitioned parquet for chunked access
      chunk_001.parquet
      ...
```

---

## Part 2: Benchmarking Plan

**File to create:** `benchmarks/benchmark_dlnm.R`

### Model Configurations to Test

| Config | Var Basis | Lag Basis | Lag Range | CB Columns | Purpose |
|--------|-----------|-----------|-----------|------------|---------|
| C1 | `lin` (df=1) | `poly(degree=4)` (df=5) | 0-15 | 5 | Minimal/DLM |
| C2 | `ns(df=5)` | `ns(df=4)` | 0-21 | 20 | Typical epidemiology |
| C3 | `bs(df=6)` | `ns(df=4)` | 0-40 | 24 | Extended lag |
| C4 | `ps(df=10)` | `ps(df=5)` | 0-30 | 50 | Penalized high-df |
| C5 | `ps(df=15)` | `ps(df=8)` | 0-60 | 120 | Stress test |

### Stages to Benchmark Separately
For each (scale x config) combination, time these independently:

1. **`crossbasis()`** -- cross-basis matrix construction
2. **`glm()`** -- model fitting (`death ~ cb + ns(time, 7*n_years) + dow`, family=quasipoisson)
3. **`crosspred()`** -- predictions at 50 exposure values, with `cumul=TRUE`
4. **`crossreduce()`** -- overall cumulative reduction

### Sub-function Profiling
Create instrumented wrappers to time internal operations within `crossbasis()` and `crosspred()`:

**Inside `crossbasis()`:**
- `onebasis()` for variable basis (`R/crossbasis.R` line 31)
- `onebasis()` for lag basis (`R/crossbasis.R` line 49)
- `Lag()` matrix construction (`R/crossbasis.R` line 64, from `tsModel` package)
- `mat %*% basislag[,l]` multiply loop (`R/crossbasis.R` lines 66-68)

**Inside `crosspred()`:**
- `mkXpred()` call including `tensor.prod.model.matrix()` (`R/crosspred.R` line 101)
- `Xpred %*% coef` prediction (`R/crosspred.R` line 104)
- `rowSums((Xpred %*% vcov) * Xpred)` SE quadratic form (`R/crosspred.R` line 105)
- Cumulative accumulation loop (`R/crosspred.R` lines 126-133)

### Tooling
- Use `bench::mark()` for timing (provides GC info, memory, iterations)
- Use `profmem::profmem()` for allocation tracking
- Use `gc()` deltas for peak memory measurement
- 3-5 iterations for large stages, 10+ for fast stages

### 10 GB Special Handling
The 10 GB dataset cannot use the standard pipeline (cross-basis alone would be 48GB for 50 columns). Options:
- **Benchmark only crossbasis construction** using chunked reads via `arrow::open_dataset()`
- **Time extrapolation**: fit scaling curves from 10MB/100MB/1GB results to predict 10GB behavior
- **This becomes the primary motivation for the Rust streaming implementation** (Part 4)

### Output
```
benchmarks/results/
  timing_results.csv        # scale, config, stage, median_time_sec, mem_alloc_mb
  substage_timings.csv      # scale, config, substage, time_seconds
  scaling_analysis.csv      # for regression: stage, config, n_rows, time_sec
```

---

## Part 3: Rust/SIMD Optimization Opportunities

### P1: Cross-Basis Kernel (Highest Priority)

**File:** `R/crossbasis.R` lines 61-69

**Current code:**
```r
crossbasis <- matrix(0, nrow=dim[1], ncol=ncol(basisvar)*ncol(basislag))
for(v in seq(length=ncol(basisvar))) {
    if(dim[2]==1L) {
      mat <- as.matrix(Lag(basisvar[, v], seqlag(lag), group=group))
    } else mat <- matrix(basisvar[,v], ncol=diff(lag)+1)
    for(l in seq(length=ncol(basislag))) {
      crossbasis[,ncol(basislag)*(v-1)+l] <- mat %*% (basislag[,l])
    }
}
```

**Current:** Nested R loop -- for each of `nv` variable basis columns, calls `Lag()` to create (n x lag_range) matrix, then for each of `nl` lag basis columns, does `mat %*% basislag[,l]`. Total: `nv * nl` BLAS calls + `nv` large matrix allocations.

**Rust approach:** Fused sliding-window dot product. Instead of materializing the lag matrix, compute `sum_j(basisvar[i-j, v] * basislag[j, l])` directly for each row. This is a convolution that maps perfectly to SIMD horizontal sums (AVX2: 4 f64/cycle). Batch all (v,l) output columns in a single pass. Parallelize rows with `rayon`.

**Expected speedup:** 10-50x
**Memory reduction:** Eliminates nv lag matrices (each n x lag_range)

---

### P2: Quadratic Form for Standard Errors (High Priority)

**File:** `R/crosspred.R` lines 105, 131, 135

**Current code:**
```r
# Line 105: lag-specific SE
matse <- matrix(sqrt(pmax(0,rowSums((Xpred%*%vcov)*Xpred))), ...)

# Line 131: cumulative SE (inside loop over lags)
cumse[, i] <- sqrt(pmax(0,rowSums((Xpredall%*%vcov)*Xpredall)))

# Line 135: overall SE
allse <- sqrt(pmax(0,rowSums((Xpredall%*%vcov)*Xpredall)))
```

**Current:** `rowSums((X %*% V) * X)` computes `diag(X V X')` by forming full (m x p) temporary. Called once for lag-specific SEs, once per lag in cumulative loop, and once for overall. The cumulative loop (lines 126-133) recomputes from scratch at each lag.

**Rust approach:** Row-wise `x_i' V x_i` with SIMD inner products, no temporary matrix. For cumulative variant, use incremental update: `(x_old + x_delta)' V (x_old + x_delta) = x_old'Vx_old + 2*x_old'Vx_delta + x_delta'Vx_delta` avoids full recomputation per lag.

**Expected speedup:** 5-20x

---

### P3: Tensor Product Model Matrix (Medium Priority)

**File:** `R/mkXpred.R` line 24

**Current code:**
```r
Xpred <- tensor.prod.model.matrix(list(basisvar, basislag))
```

**Current:** `tensor.prod.model.matrix()` from mgcv computes row-wise Kronecker product of two matrices.

**Rust approach:** Direct row-wise outer product with SIMD vectorization. Trivially parallelizable.

**Expected speedup:** 3-10x

---

### P4: B-Spline Basis Evaluation (Medium Priority)

**File:** `R/ps.R` line 33

**Current code:**
```r
basis <- splineDesign(knots, x, degree+1, x*0, TRUE)
```

**Current:** de Boor's algorithm via R's `splines` C code.

**Rust approach:** Vectorized de Boor processing multiple evaluation points simultaneously via SIMD. Batch evaluation eliminates per-point R overhead.

**Expected speedup:** 3-8x

---

### P5: Lag Matrix Construction (Fused with P1)

**File:** `R/crossbasis.R` line 64, `tsModel::Lag()`

**Current:** Creates (n x lag_range) matrix per variable basis column by copying+shifting.

**Rust approach:** Eliminated entirely by the fused kernel in P1 -- no lag matrix materialization needed.

---

### P6: Penalty Matrix Operations (Lower Priority)

**File:** `R/cbPen.R` lines 28-33

**Current code:**
```r
Slist <- c(Slist,list(Svar=attr$argvar$S %x% diag(attr$df[2])))
Slist <- c(Slist,list(Slag=diag(attr$df[1]) %x% attr$arglag$S))
Slist <- lapply(Slist, function(X)
    X/eigen(X, symmetric=TRUE, only.values=TRUE)$values[1])
```

**Current:** Kronecker products `S %x% diag(d)` creating dense matrices, eigenvalue decomposition for rescaling.

**Rust approach:** Block-diagonal storage for structured Kronecker products. Power iteration for leading eigenvalue only (O(p^2) vs O(p^3)).

**Expected speedup:** 2-5x (one-time setup cost only)

---

### Priority Summary

| # | Target | Speedup | Impact at Scale | Complexity |
|---|--------|---------|-----------------|------------|
| P1 | Cross-basis kernel | 10-50x | Very High | Medium |
| P2 | Quadratic form SE | 5-20x | High | Low-Medium |
| P3 | Tensor product | 3-10x | Medium | Low |
| P4 | B-spline basis | 3-8x | Medium | Medium |
| P5 | Lag matrix | Fused w/ P1 | Medium | Low |
| P6 | Penalty matrices | 2-5x | Low | Medium-High |

---

## Part 4: 10 GB Strategy -- Streaming Rust Architecture

The 10 GB case cannot fit in R's memory model. The cross-basis matrix alone at 120M rows x 50 cols = 48 GB.

### Solution: Streaming Sufficient Statistics

For Poisson GLM, IRLS needs only the accumulated `X'WX` (p x p) and `X'Wy` (p x 1) matrices. Rust processes the data in ~1M-row chunks:

1. Memory-map the data file (Rust `memmap2` crate)
2. For each chunk: compute basis functions, fused cross-basis (no full matrix), accumulate `X'WX` and `X'Wy`
3. Return small sufficient statistics to R for final solve

### Memory Budget Per Chunk (~1M rows)

| Component | Size | Notes |
|-----------|------|-------|
| Memory-mapped data | 10 GB virtual | OS pages in/out as needed |
| One chunk raw data | ~85 MB | 1M rows x 14 cols |
| Basisvar for chunk | 8-120 MB | 1M x (1-15 df) x 8 bytes |
| Basislag (constant) | <50 KB | lag_range x nl x 8 bytes |
| Cross-basis chunk | 40-960 MB | 1M x (5-120) x 8 bytes |
| X'WX accumulator | <120 KB | p x p x 8 bytes |
| **Total working memory** | **~200 MB - 1.1 GB** | vs 48 GB+ for full materialization |

---

## Implementation Order

1. `benchmarks/generate_data.R` -- data generation script
2. `benchmarks/benchmark_dlnm.R` -- benchmarking script with all configs
3. Run baseline R benchmarks at 10MB, 100MB, 1GB
4. Analyze scaling behavior, produce timing/memory reports
5. Set up Rust crate with `extendr` (or `savvy`) bindings
6. Implement P1 (cross-basis kernel) and P2 (quadratic form) in Rust
7. Re-benchmark with Rust kernels
8. Implement streaming pipeline for 10GB case

## Verification

- **Data generation:** Compare summary statistics (mean, sd, correlation structure) of generated data vs original `chicagoNMMAPS`
- **Benchmarks:** Validate that model results (coefficients, SEs) match across all scales
- **Rust correctness:** Compare Rust kernel outputs to R outputs with `all.equal(r_result, rust_result, tolerance=1e-10)`
- **Run:** `Rscript benchmarks/generate_data.R` then `Rscript benchmarks/benchmark_dlnm.R`

## Critical Files

| File | Role |
|------|------|
| `R/crossbasis.R` | Primary bottleneck: nested loop lines 61-69 |
| `R/crosspred.R` | Quadratic form SE: lines 104-106, 126-135 |
| `R/mkXpred.R` | Tensor product: line 24 |
| `R/ps.R` | B-spline basis: line 33 |
| `R/cbPen.R` | Penalty matrices: lines 28-33 |
| `R/crossreduce.R` | Dimension reduction with matrix ops |
| `R/onebasis.R` | Basis function dispatcher |
| `data/chicagoNMMAPS.rda` | Base dataset for replication |
