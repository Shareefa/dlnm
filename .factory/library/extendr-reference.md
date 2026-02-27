# Extendr/Rextendr Reference

Reference for integrating Rust with the dlnm R package via extendr.

---

## Setup Process

1. `rextendr::use_extendr()` creates scaffolding:
   - `R/extendr-wrappers.R` (auto-generated, DO NOT EDIT)
   - `src/Makevars`, `src/Makevars.win`
   - `src/entrypoint.c`
   - `src/rust/Cargo.toml`
   - `src/rust/src/lib.rs`

2. `rextendr::document()` compiles Rust, generates R wrappers, runs roxygen2.

## Cargo.toml Template

```toml
[package]
name = 'dlnm'
version = '0.1.0'
edition = '2021'

[lib]
crate-type = ['staticlib']

[dependencies]
extendr-api = '0.8'
rayon = '1.10'
```

## Function Annotation

```rust
use extendr_api::prelude::*;

/// Brief description for R docs
/// @export
#[extendr]
fn my_function(mat: RMatrix<f64>) -> RMatrix<f64> {
    // Implementation
}

extendr_module! {
    mod dlnm;
    fn my_function;
}
```

## Matrix Operations (CRITICAL)

### Zero-Copy Access to R Matrix Data
```rust
#[extendr]
fn process_matrix(mat: RMatrix<f64>) -> Doubles {
    let nrows = mat.nrows();
    let ncols = mat.ncols();
    let data = mat.data();  // Zero-copy: direct pointer to R memory
    // R matrices are COLUMN-MAJOR: index = row + col * nrows
    // ...
}
```

### Creating Result Matrices
```rust
RMatrix::new_matrix(nrows, ncols, |row, col| {
    // compute value at (row, col)
    0.0_f64
})
```

## Key Types

| Rust Type | R Type |
|-----------|--------|
| `f64` | numeric scalar |
| `i32` | integer scalar |
| `Doubles` | numeric vector |
| `Integers` | integer vector |
| `RMatrix<f64>` | numeric matrix |
| `&[f64]` | numeric slice (borrowed, zero-copy) |
| `Rfloat` | NA-aware double |

## Critical Gotchas

1. **Column-major order**: R matrices store columns contiguously. `data[row + col * nrows]`
2. **NA handling**: Use `Rfloat`/`Rint` for NA-aware; raw `f64` ignores NA
3. **Every `#[extendr]` fn must be in `extendr_module!`**
4. **Always run `rextendr::document()` after Rust changes**
5. **ARM NEON**: Use `-C target-cpu=native` in `.cargo/config.toml` for auto-vectorization
6. **Rayon**: Use `rayon` crate for data parallelism across rows

## Workflow

```r
# After editing src/rust/src/lib.rs:
rextendr::document()       # Compile Rust + generate R wrappers
pkgload::load_all(".")     # Reload package
# Test functions
```
