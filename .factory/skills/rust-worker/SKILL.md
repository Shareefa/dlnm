---
name: rust-worker
description: Implements Rust optimizations for dlnm R package via extendr and wires them into R
---

# Rust Optimization Worker

NOTE: Startup and cleanup are handled by `worker-base`. This skill defines the WORK PROCEDURE.

## When to Use This Skill

Use for features involving:
- Setting up the extendr/rextendr R-Rust bridge
- Implementing Rust kernels (P1 cross-basis, P2 quadratic form SE)
- Wiring Rust functions into existing R functions with fallback
- Correctness testing of Rust vs R implementations

## Work Procedure

### Step 1: Understand the Feature

1. Read the feature description, preconditions, and expectedBehavior carefully.
2. Read `.factory/library/extendr-reference.md` for the extendr API reference.
3. Read `.factory/library/architecture.md` for the code structure and optimization targets.
4. Read the relevant R source files (`R/crossbasis.R`, `R/crosspred.R`) to understand exactly what code is being optimized.
5. Read `benchmarks/MODEL_EXPLAINER.md` for mathematical context of the cross-basis computation.

### Step 2: Write Tests First (TDD)

1. Create test files in `tests/testthat/` BEFORE implementing Rust code.
2. For Rust function correctness:
   - Compute expected results using the existing R implementation
   - Write `test_that()` blocks comparing Rust output to R output via `all.equal(tolerance=1e-10)`
   - Cover all model configs (C1-C5) where applicable
   - Include edge cases: NA values, single-group (no grouping), group boundaries, empty groups
3. For R integration:
   - Test that the dispatch mechanism selects Rust when available
   - Test fallback to R when Rust is unavailable
   - Test full pipeline equivalence (crossbasis → glm → crosspred → crossreduce)
4. Run tests — they should FAIL at this point (red phase).

### Step 3: Implement Rust Code

1. All Rust code goes in `src/rust/src/` directory.
2. Use `extendr-api` for R interop — annotate functions with `#[extendr]`.
3. Register all functions in `extendr_module!` macro.
4. Key implementation constraints:
   - **Column-major layout**: R matrices are column-major. `data[row + col * nrows]`.
   - **Zero-copy input**: Use `mat.data()` for read-only access to R matrix data.
   - **ARM NEON**: Use `cargo build --release` with `-C target-cpu=native` for auto-vectorization.
   - **Rayon**: Use `rayon::prelude::*` for data parallelism across rows.
   - **No lag matrix materialization** (P1): Compute sliding-window dot products directly.
   - **NA propagation**: Match R's NA behavior exactly — NAs in lag window produce NA output.
   - **Group boundaries**: Lag windows must NOT cross city/group boundaries.
5. After implementation: `rextendr::document()` to compile and generate R wrappers.

### Step 4: Wire into R Functions

1. Modify R functions to call Rust with fallback:
   ```r
   result <- tryCatch(
     .Call("wrap__rust_function", arg1, arg2),
     error = function(e) {
       # R fallback implementation
       original_r_code(arg1, arg2)
     }
   )
   ```
2. Do NOT delete the original R code — keep it as the fallback path.
3. Ensure the function signature and return type are unchanged.
4. Run `pkgload::load_all(".", quiet=TRUE)` to reload the package.

### Step 5: Verify Correctness

1. Run all tests: `Rscript -e "testthat::test_dir('tests/testthat')"`
2. All tests must pass (green phase).
3. Run a quick manual benchmark to verify Rust is faster:
   ```r
   pkgload::load_all(".", quiet=TRUE)
   dt <- readRDS("benchmarks/data/scale_10mb.rds")
   system.time(crossbasis(dt$temp, lag=c(0,21), argvar=list(fun="ns",df=5), arglag=list(fun="ns",df=4), group=dt$city))
   ```
4. Compare timing with baseline (10MB C2 R baseline: ~0.22 seconds).

### Step 6: Run Validators

1. `Rscript -e "testthat::test_dir('tests/testthat')"` — all tests pass
2. `cd src/rust && cargo test` — Rust unit tests pass
3. `cd src/rust && cargo clippy` — no warnings

## Example Handoff

```json
{
  "salientSummary": "Implemented P1 fused cross-basis kernel in Rust via extendr. The kernel computes sliding-window dot products without materializing lag matrices, using rayon for row parallelism. Verified correctness via all.equal(tolerance=1e-10) for all 5 configs (C1-C5) on 10MB data. Quick benchmark shows 8x speedup at 10MB C2 (0.027s vs 0.22s R).",
  "whatWasImplemented": "src/rust/src/crossbasis_kernel.rs: fused_crossbasis() function taking basisvar matrix, basislag matrix, lag range, and group boundaries. Computes sum_j(basisvar[i-j, v] * basislag[j, l]) directly for each (row, v, l) combination. Uses rayon par_chunks for row parallelism. Handles NA propagation and group boundaries. R wrapper in R/crossbasis.R dispatches to Rust via tryCatch(.Call()) with R fallback.",
  "whatWasLeftUndone": "",
  "verification": {
    "commandsRun": [
      {"command": "cd src/rust && cargo build --release", "exitCode": 0, "observation": "Compiled successfully, 0 warnings"},
      {"command": "cd src/rust && cargo test", "exitCode": 0, "observation": "6 tests passed"},
      {"command": "cd src/rust && cargo clippy", "exitCode": 0, "observation": "No warnings"},
      {"command": "Rscript -e \"testthat::test_dir('tests/testthat')\"", "exitCode": 0, "observation": "12 tests passed, 0 failures"},
      {"command": "Rscript -e \"pkgload::load_all('.', quiet=TRUE); dt <- readRDS('benchmarks/data/scale_10mb.rds'); system.time(crossbasis(dt$temp, lag=c(0,21), argvar=list(fun='ns',df=5), arglag=list(fun='ns',df=4), group=dt$city))\"", "exitCode": 0, "observation": "elapsed: 0.027 (vs 0.22s baseline = 8.1x speedup)"}
    ],
    "interactiveChecks": []
  },
  "tests": {
    "added": [
      {"file": "tests/testthat/test-rust-crossbasis.R", "cases": [
        {"name": "P1 matches R for C1", "verifies": "Rust fused kernel output equals R Lag+multiply for C1 config"},
        {"name": "P1 matches R for C2", "verifies": "Same for C2"},
        {"name": "P1 handles NAs", "verifies": "NA propagation matches R"},
        {"name": "P1 handles group boundaries", "verifies": "No cross-city lag leakage"},
        {"name": "P1 single group (no grouping)", "verifies": "Works with group=NULL"}
      ]}
    ]
  },
  "discoveredIssues": []
}
```

## When to Return to Orchestrator

- The Rust crate fails to compile due to missing system dependencies (e.g., missing C toolchain)
- `rextendr::use_extendr()` or `rextendr::document()` fails with unclear errors
- R package loading fails after Rust integration
- Numerical discrepancies exceed tolerance=1e-10 and root cause is unclear
- Memory issues during testing at 100MB+ scale
- Feature requires modifying R function signatures (API changes not allowed)
