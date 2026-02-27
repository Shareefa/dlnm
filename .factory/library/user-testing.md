# User Testing

Testing surface: tools, URLs, setup steps, isolation notes, known quirks.

**What belongs here:** How to test the mission outputs, what tools to use, setup for validation.

---

## Testing Surface

This is a CLI-based computational research project. There is no web UI or TUI. All validation is done via:

1. **Rscript commands** — Run R scripts and check output/exit codes
2. **File existence and content checks** — Verify CSVs, PNGs, markdown files exist with expected content
3. **Numerical comparison** — Use `all.equal()` in R to compare Rust vs R outputs

## How to Run Tests

### Rust Correctness Tests
```bash
# After Rust is wired in:
cd /Users/abdullahshareef/Documents/Projects/dlnm
Rscript -e "testthat::test_dir('tests/testthat')"
```

### Benchmarks
```bash
cd /Users/abdullahshareef/Documents/Projects/dlnm
Rscript benchmarks/benchmark_dlnm.R 10mb C1,C2
```

### Package Load Test
```bash
Rscript -e 'pkgload::load_all(".", quiet=TRUE); cat("OK\n")'
```

## Validation Approach

- **Correctness:** Run `all.equal(rust_result, r_result, tolerance=1e-10)` for all configs
- **Performance:** Compare timing CSVs (baseline in `benchmarks/results/timing_results.csv`)
- **Report:** Check file existence and section headers in markdown

## Known Quirks

- `pkgload::load_all()` may show dlnm version message — this is normal
- At 1GB+ scale, benchmarks take several minutes per config
- C4/C5 configs may not be feasible at 10GB scale due to memory constraints
- The 10GB dataset is parquet-partitioned, requiring `arrow` package for reading
