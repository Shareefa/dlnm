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

### Python (use project venv)
```bash
.venv/bin/python3 benchmarks/generate_plots.py
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

---

## Flow Validator Guidance: CLI

This project is entirely CLI-based. All testing is done via `Rscript` commands. There is no web UI, API server, or TUI.

### Testing Approach
- Use `Execute` tool to run `Rscript -e '...'` commands
- Each assertion is tested by running R code that produces terminal output (TRUE/FALSE, numeric values, etc.)
- Capture the terminal output as evidence for pass/fail determination

### Isolation Rules
- All subagents test against the same installed package — this is safe because tests are read-only (no mutations)
- Each subagent can run `pkgload::load_all(".", quiet=TRUE)` independently
- Benchmark data in `benchmarks/data/` is read-only shared data
- No test accounts or namespaces needed — CLI tests don't have state conflicts

### Key Commands
- Load package: `pkgload::load_all(".", quiet=TRUE)`
- The 10MB benchmark data: `readRDS("benchmarks/data/scale_10mb.rds")`
- Rust is already compiled; no need to rebuild
- Working directory: `/Users/abdullahshareef/Documents/Projects/dlnm`

### Config Definitions (C1-C5)
These are the basis function configurations used in DLNM testing:
- **C1**: lin x poly(4), lag 0-15
- **C2**: ns(5) x ns(4), lag 0-21
- **C3**: bs(6) x ns(4), lag 0-40
- **C4**: ps(10) x ps(5), lag 0-30
- **C5**: ps(15) x ps(8), lag 0-60

### Report Format
Write flow report as JSON to the specified path with this structure:
```json
{
  "flowId": "<group-id>",
  "assertions": {
    "<assertion-id>": {
      "status": "pass|fail|blocked",
      "evidence": "terminal output or observation",
      "reason": "why it passed/failed/blocked (if not pass)"
    }
  },
  "frictions": [],
  "blockers": [],
  "toolsUsed": ["Execute"]
}
```
