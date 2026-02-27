use extendr_api::prelude::*;

/// Compute the dot product of two numeric vectors.
/// This is a trivial test function to verify the R-Rust bridge works.
/// @param x A numeric vector.
/// @param y A numeric vector of the same length as x.
/// @return A scalar numeric value representing the dot product.
#[extendr]
fn rust_dot_product(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "Vectors must have the same length");
    x.iter().zip(y.iter()).map(|(a, b)| a * b).sum()
}

/// Sum elements of a numeric vector using Rayon parallel reduction.
/// This verifies that the rayon dependency is correctly linked.
/// @param x A numeric vector.
/// @return A scalar numeric value representing the sum.
#[extendr]
fn rust_parallel_sum(x: &[f64]) -> f64 {
    use rayon::prelude::*;
    x.par_iter().sum()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod dlnm;
    fn rust_dot_product;
    fn rust_parallel_sum;
}
