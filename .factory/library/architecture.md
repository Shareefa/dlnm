# Architecture

Architectural decisions, patterns discovered, and design notes.

**What belongs here:** Key design decisions, code organization patterns, module boundaries.

---

## DLNM Package Structure

```
dlnm/
├── R/                    # R source functions
│   ├── crossbasis.R      # PRIMARY BOTTLENECK — builds cross-basis matrix
│   ├── crosspred.R       # Predictions + SE computation (P2 target)
│   ├── crossreduce.R     # Dimension reduction
│   ├── onebasis.R        # Basis function dispatcher
│   ├── mkXpred.R         # Tensor product model matrix
│   ├── ps.R              # P-spline basis (B-spline wrapper)
│   ├── cbPen.R           # Penalty matrix operations
│   └── ...               # Other R files
├── src/                  # Compiled code (Rust via extendr)
│   └── rust/             # Rust crate (created by rextendr::use_extendr())
│       ├── Cargo.toml
│       └── src/
│           └── lib.rs    # Rust implementations of P1, P2
├── benchmarks/           # Benchmark scripts and data
│   ├── benchmark_dlnm.R  # Main benchmark script
│   ├── generate_data.R   # Data generation
│   ├── data/             # Generated datasets (.gitignored)
│   └── results/          # Benchmark results (.gitignored)
├── tests/                # testthat tests
│   └── testthat/
└── DESCRIPTION           # Package metadata
```

## Optimization Targets

### P1: Cross-Basis Kernel (crossbasis.R lines 61-69)
- **Current:** Nested R loop — for each variable basis column, calls `Lag()` to create (n × lag_range) matrix, then multiplies by each lag basis column
- **Rust:** Fused sliding-window dot product — compute `sum_j(basisvar[i-j, v] * basislag[j, l])` directly for each row without materializing lag matrices
- **Key constraint:** Must handle `group` boundaries (cities) — lag windows must NOT cross city boundaries

### P2: Quadratic Form SE (crosspred.R lines 105, 131, 135)
- **Current:** `rowSums((X %*% V) * X)` forms full (m × p) temporary matrix
- **Rust:** Row-wise `x_i' V x_i` with no temporary. Cumulative variant uses incremental update.

## R-Rust Integration Pattern

- Use `rextendr` for scaffolding
- Rust functions exposed via `#[extendr]` macro
- R wrappers call Rust via `.Call()` with automatic fallback to R if Rust unavailable
- Detection: `tryCatch(.Call("wrap__rust_function"), error = function(e) use_r_fallback())`

## Model Configurations

| Config | Var Basis | Lag Basis | Lag Range | CB Columns |
|--------|-----------|-----------|-----------|------------|
| C1 | lin (1) | poly(4) (5) | 0-15 | 5 |
| C2 | ns(5) (5) | ns(4) (4) | 0-21 | 20 |
| C3 | bs(6) (6) | ns(4) (4) | 0-40 | 24 |
| C4 | ps(10) (10) | ps(5) (5) | 0-30 | 50 |
| C5 | ps(15) (15) | ps(8) (8) | 0-60 | 120 |
