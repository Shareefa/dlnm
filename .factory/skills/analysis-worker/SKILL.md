---
name: analysis-worker
description: Implements R/Python analysis scripts for benchmarking, mixmeta pipeline, convergence analysis, and report generation
---

# Analysis Worker

NOTE: Startup and cleanup are handled by `worker-base`. This skill defines the WORK PROCEDURE.

## When to Use This Skill

Use for features involving:
- Running benchmarks with the optimized dlnm package
- Setting up the mixmeta two-stage pipeline
- Generating synthetic truth data with known exposure-lag-response surfaces
- Running comparison analyses (single-stage vs two-stage)
- Computing convergence and statistical accuracy metrics
- Generating PNG plots and markdown reports

## Work Procedure

### Step 1: Understand the Feature

1. Read the feature description, preconditions, and expectedBehavior carefully.
2. Read `.factory/library/architecture.md` for code structure.
3. Read `benchmarks/MODEL_EXPLAINER.md` for mathematical context.
4. Read existing benchmark scripts (`benchmarks/benchmark_dlnm.R`, `benchmarks/run_1gb_c1c2c3.R`) to follow established patterns.
5. Check existing results in `benchmarks/results/timing_results.csv` for baseline reference.

### Step 2: Write Tests First (TDD)

1. For benchmark scripts: Write a small test that runs the script on 10MB with one config and verifies CSV output structure.
2. For analysis scripts: Write tests that verify metric computation on small known datasets.
3. For data generation: Write tests that verify generated data has expected properties (dimensions, column names, value ranges).
4. For plots: Write tests that verify PNG file creation.
5. For reports: Write tests that verify markdown file existence and section headers.
6. Run tests — they should FAIL initially (red phase).

### Step 3: Implement

#### For Benchmark Scripts:
1. Load the package with `pkgload::load_all(".", quiet=TRUE)`.
2. Follow the existing `benchmark_dlnm.R` pattern for timing and output format.
3. Save results to `benchmarks/results/` with `optimized_` prefix for Rust-optimized runs.
4. Include both `crossbasis` and full pipeline stages where feasible.
5. For 10GB: use `arrow::open_dataset()` for parquet reading; crossbasis-only if GLM won't fit in memory.

#### For Mixmeta Pipeline:
1. Install mixmeta if not present: `install.packages("mixmeta", repos="https://cloud.r-project.org")`.
2. Implement the two-stage workflow:
   ```r
   # Stage 1: Per-city DLNM
   for each city:
     cb <- crossbasis(city_temp, lag=..., argvar=..., arglag=...) 
     model <- glm(death ~ cb + ns(time, df=...) + dow, family=quasipoisson(), data=city_data)
     reduced <- crossreduce(cb, model, type="overall")
     store coef(reduced) and vcov(reduced)
   # Stage 2: Pool with mixmeta
   pooled <- mixmeta(coef_matrix ~ 1, S=vcov_list)
   ```
3. Record timing for each stage separately.

#### For Truth Data Generation:
1. Define a TRUE exposure-lag-response surface mathematically:
   - Choose a known bi-linear or spline-based surface: `f(temp, lag) = sum_vl(beta_true[v,l] * basis_v(temp) * basis_l(lag))`
   - Set explicit `beta_true` coefficient values
   - Generate outcomes: `death ~ Poisson(exp(intercept + crossbasis %*% beta_true + confounders))`
2. Include confounders (time trend, day-of-week) with known coefficients.
3. Generate multi-city data with between-city heterogeneity in the coefficients.
4. Save true coefficients alongside the data for convergence comparison.
5. Document the DGP specification in script comments.

#### For Convergence Analysis:
1. Run both approaches on truth data:
   - Single-stage: `crossbasis → glm → crosspred` on full multi-city data
   - Two-stage: per-city DLNM → crossreduce → mixmeta pool
2. Extract estimated coefficients/RR curves from each approach.
3. Compare against known true surface:
   - Bias = estimated - true (at each temperature-lag point)
   - MSE = mean((estimated - true)^2)
   - CI coverage = proportion of true values within 95% CI
   - RMSE = sqrt(MSE)
   - Relative bias = bias / true_value (at key temperature points)
4. Record numerical convergence: GLM iteration count, convergence status.
5. Save all results to CSV.

#### For Plots (Python):
1. Read CSV results with pandas.
2. Create plots with matplotlib:
   - Speedup comparison bar chart (R vs Rust across scales/configs)
   - Timing scaling curves (log-log with scaling exponent)
   - Convergence/bias plots (bias vs temperature for both approaches)
   - CI coverage comparison
   - Memory usage comparison (if data available)
3. Save as PNG to `benchmarks/results/report/`.
4. Use clear labels, titles, and legends.

#### For Report (Markdown):
1. Create `benchmarks/results/report/report.md`.
2. Required sections: Introduction, Optimization Methodology, Benchmark Results, Mixmeta Comparison, Convergence Analysis, Memory and Scaling, Conclusions.
3. Embed PNG plots using relative paths: `![Caption](filename.png)`.
4. Include formatted tables with timing comparisons.
5. Reference actual CSV data — do NOT fabricate numbers.

### Step 4: Verify

1. Run tests: `Rscript -e "testthat::test_dir('tests/testthat')"`
2. Verify output files exist and contain expected content.
3. For benchmarks: spot-check timing values are reasonable (not 0, not negative, not absurdly large).
4. For convergence: verify metrics are within reasonable ranges (e.g., CI coverage between 0 and 1).
5. For plots: verify PNG files are non-empty.
6. For report: verify all section headers present and tables formatted correctly.

## Example Handoff

```json
{
  "salientSummary": "Built two-stage mixmeta pipeline script and ran it at 10MB and 100MB scales. Per-city DLNM fits (24 cities at 10MB) complete in 3.2s total, pooling with mixmeta in 0.1s. At 100MB (235 cities), per-city fitting takes 31s. Heterogeneity I²=42%, τ²=0.015. Pooled cumulative RR curve matches expected shape.",
  "whatWasImplemented": "benchmarks/run_mixmeta_pipeline.R: Full two-stage pipeline (per-city crossbasis → glm → crossreduce → mixmeta pool). Records per-city timing, convergence status, and pooled estimates. Outputs timing to benchmarks/results/mixmeta_timing.csv and estimates to benchmarks/results/mixmeta_estimates.csv. tests/testthat/test-mixmeta-pipeline.R with 5 test cases.",
  "whatWasLeftUndone": "",
  "verification": {
    "commandsRun": [
      {"command": "Rscript benchmarks/run_mixmeta_pipeline.R 10mb C2", "exitCode": 0, "observation": "24 cities fitted, pooling complete, CSV saved"},
      {"command": "Rscript benchmarks/run_mixmeta_pipeline.R 100mb C2", "exitCode": 0, "observation": "235 cities fitted in 31.2s, pooling in 0.08s"},
      {"command": "Rscript -e \"testthat::test_dir('tests/testthat')\"", "exitCode": 0, "observation": "18 tests passed, 0 failures"}
    ],
    "interactiveChecks": [
      {"action": "Verified mixmeta_timing.csv has expected columns", "observed": "scale,config,stage,n_cities,time_sec columns present, 4 rows for 10mb+100mb"},
      {"action": "Checked pooled RR at temperature 30°C", "observed": "RR=1.23 (95% CI: 1.15-1.32), consistent with expected temperature-mortality relationship"}
    ]
  },
  "tests": {
    "added": [
      {"file": "tests/testthat/test-mixmeta-pipeline.R", "cases": [
        {"name": "pipeline produces pooled estimates", "verifies": "Two-stage pipeline completes and returns mixmeta object"},
        {"name": "per-city models converge", "verifies": "All city GLMs converge (convergence=TRUE)"},
        {"name": "heterogeneity stats present", "verifies": "I² and tau² are computed and in plausible ranges"},
        {"name": "timing CSV has expected columns", "verifies": "Output CSV schema matches specification"},
        {"name": "pipeline runs at 10MB and 100MB", "verifies": "Multi-scale execution completes without error"}
      ]}
    ]
  },
  "discoveredIssues": []
}
```

## When to Return to Orchestrator

- Rust-optimized package fails to load (needs rust-worker to fix)
- Benchmark data files are missing or corrupted
- 10GB benchmarks hit OOM — need orchestrator guidance on what to skip
- mixmeta package installation fails
- Truth data generation produces degenerate data (all zeros, all NAs)
- Model fitting fails to converge at large scale
- Requirements are ambiguous about what configs or scales to run
