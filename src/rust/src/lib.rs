use extendr_api::prelude::*;
use rayon::prelude::*;

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
    x.par_iter().sum()
}

/// Fused cross-basis kernel: computes sliding-window dot products without
/// materializing lag matrices.
///
/// For each output element cb[i, nl*(v-1)+l], computes:
///   sum_j( basisvar[i-j, v] * basislag[j, l] )
/// where j ranges from 0 to (lag_max - lag_min), respecting group boundaries.
///
/// Group boundaries: lag windows must not cross city/group boundaries.
/// NA handling: if any basisvar[i-j, v] is NA (NaN), the output is NA.
///
/// @param basisvar Numeric matrix (n x nv) - variable basis.
/// @param basislag Numeric matrix (lag_range x nl) - lag basis.
/// @param lag_min Integer - minimum lag (typically 0).
/// @param lag_max Integer - maximum lag.
/// @param group_starts Integer vector - 1-based start indices of each group.
/// @param group_ends Integer vector - 1-based end indices of each group.
/// @return Numeric matrix (n x nv*nl) - the cross-basis matrix.
/// @export
#[extendr]
fn fused_crossbasis(
    basisvar: RMatrix<f64>,
    basislag: RMatrix<f64>,
    lag_min: i32,
    lag_max: i32,
    group_starts: Integers,
    group_ends: Integers,
) -> Robj {
    let n = basisvar.nrows();
    let nv = basisvar.ncols();
    let lag_len = basislag.nrows(); // number of lag points = lag_max - lag_min + 1
    let nl = basislag.ncols();

    // Validate lag parameters match basislag dimensions
    let expected_lag_len = (lag_max - lag_min + 1) as usize;
    assert_eq!(
        lag_len, expected_lag_len,
        "basislag nrows ({}) must equal lag_max - lag_min + 1 ({})",
        lag_len, expected_lag_len
    );
    let ncb = nv * nl; // total output columns

    // Access raw data (column-major layout)
    let bv_data = basisvar.data();
    let bl_data = basislag.data();

    // Convert 1-based R indices to 0-based
    let n_groups = group_starts.len();
    let g_starts: Vec<usize> = (0..n_groups)
        .map(|i| {
            let val = group_starts.elt(i);
            if val.is_na() {
                panic!("group_starts contains NA");
            }
            (val.inner() - 1) as usize
        })
        .collect();
    let g_ends: Vec<usize> = (0..n_groups)
        .map(|i| {
            let val = group_ends.elt(i);
            if val.is_na() {
                panic!("group_ends contains NA");
            }
            (val.inner() - 1) as usize
        })
        .collect();

    // Build a lookup: for each row i, store (group_start_0based, group_end_0based)
    let mut row_group: Vec<(usize, usize)> = vec![(0, 0); n];
    for g in 0..n_groups {
        let gs = g_starts[g];
        let ge = g_ends[g];
        for item in row_group.iter_mut().take(ge + 1).skip(gs) {
            *item = (gs, ge);
        }
    }

    // Pre-compute per-row results in a row-major temporary buffer,
    // then transpose to column-major for R.
    // Row-major: row_buf[i * ncb + col]
    let row_buf: Vec<f64> = (0..n)
        .into_par_iter()
        .flat_map_iter(|i| {
            let (grp_start, grp_end) = row_group[i];
            let mut row = vec![0.0_f64; ncb];

            for v in 0..nv {
                // Check if any lag position for this (row, v) is NA or out-of-bounds
                // This determines NA status for all lag basis columns of this v
                let mut var_vals: Vec<f64> = Vec::with_capacity(lag_len);
                let mut has_na = false;

                for j in 0..lag_len {
                    let lag_offset = lag_min as isize + j as isize;
                    let src_row = i as isize - lag_offset;

                    if src_row < grp_start as isize || src_row > grp_end as isize {
                        has_na = true;
                        break;
                    }

                    let src_row_usize = src_row as usize;
                    let bv_val = bv_data[src_row_usize + v * n];

                    if bv_val.is_nan() {
                        has_na = true;
                        break;
                    }

                    var_vals.push(bv_val);
                }

                if has_na {
                    // All lag basis columns for this v are NA
                    for l in 0..nl {
                        row[nl * v + l] = f64::NAN;
                    }
                } else {
                    // Compute dot products for each lag basis column
                    for l in 0..nl {
                        let mut sum = 0.0_f64;
                        for j in 0..lag_len {
                            sum += var_vals[j] * bl_data[j + l * lag_len];
                        }
                        row[nl * v + l] = sum;
                    }
                }
            }

            row.into_iter()
        })
        .collect();

    // Create R matrix (column-major) from row-major buffer
    let result = RMatrix::new_matrix(n, ncb, |row, col| row_buf[row * ncb + col]);

    result.into_robj()
}

/// Compute row-wise quadratic form SE: sqrt(max(0, x_i' V x_i)) for each row i.
///
/// This replaces the R pattern: sqrt(pmax(0, rowSums((X %*% V) * X)))
/// without materializing the full (m x p) temporary matrix X %*% V.
///
/// For each row i of xpred (length p), computes:
///   se_i = sqrt(max(0, sum_j sum_k x_ij * V_jk * x_ik))
///
/// @param xpred Numeric matrix (m x p) - prediction matrix.
/// @param vcov Numeric matrix (p x p) - variance-covariance matrix.
/// @return Numeric vector of length m - the standard errors.
/// @export
#[extendr]
fn quad_form_se(xpred: RMatrix<f64>, vcov: RMatrix<f64>) -> Doubles {
    let m = xpred.nrows();
    let p = xpred.ncols();

    assert_eq!(vcov.nrows(), p, "vcov nrows ({}) must equal xpred ncols ({})", vcov.nrows(), p);
    assert_eq!(vcov.ncols(), p, "vcov must be square ({}x{})", vcov.nrows(), vcov.ncols());

    let x_data = xpred.data();   // column-major: x_data[row + col * m]
    let v_data = vcov.data();    // column-major: v_data[row + col * p]

    // Compute row-wise x_i' V x_i using rayon for parallelism
    let results: Vec<f64> = (0..m)
        .into_par_iter()
        .map(|i| {
            let mut qf = 0.0_f64;
            // For each pair (j, k), accumulate x_ij * V_jk * x_ik
            // Optimize: loop over j, compute (V * x_i) for column j, then dot with x_ij
            for j in 0..p {
                let x_ij = x_data[i + j * m];
                let mut vx_j = 0.0_f64;
                for k in 0..p {
                    vx_j += v_data[j + k * p] * x_data[i + k * m];
                }
                qf += x_ij * vx_j;
            }
            // pmax(0, ...) clamping for near-singular vcov
            if qf < 0.0 { 0.0_f64 } else { qf.sqrt() }
        })
        .collect();

    Doubles::from_values(results)
}

/// Compute incremental cumulative SE for crosspred cumulative predictions.
///
/// In the R code (crosspred.R lines 126-135), the pattern is:
///   Xpredall <- 0
///   for (i in seq(length(predlag))) {
///     ind <- seq(n_at) + n_at*(i-1)
///     Xpredall <- Xpredall + Xpred[ind, , drop=FALSE]
///     cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% vcov) * Xpredall)))
///   }
///   allse <- sqrt(pmax(0, rowSums((Xpredall %*% vcov) * Xpredall)))
///
/// This Rust function computes both cumse (n_at x n_lag) and allse (n_at)
/// incrementally without materializing full temporaries.
///
/// The xpred matrix has (n_at * n_lag) rows and p columns, where rows are
/// organized as n_lag blocks of n_at rows each.
///
/// Returns an R matrix with (n_at) rows and (n_lag + 1) columns:
///   - Columns 1..n_lag contain cumulative SE at each lag step
///   - Column n_lag+1 contains the overall (final) SE (same as last cumulative column)
///
/// @param xpred Numeric matrix ((n_at * n_lag) x p) - prediction matrix.
/// @param vcov Numeric matrix (p x p) - variance-covariance matrix.
/// @param n_at Integer - number of prediction points (at values).
/// @param n_lag Integer - number of lag steps.
/// @return Numeric matrix (n_at x (n_lag + 1)) - cumulative SE + overall SE.
/// @export
#[extendr]
fn cumulative_quad_form_se(
    xpred: RMatrix<f64>,
    vcov: RMatrix<f64>,
    n_at: i32,
    n_lag: i32,
) -> Robj {
    let n_at = n_at as usize;
    let n_lag = n_lag as usize;
    let p = xpred.ncols();
    let total_rows = xpred.nrows();

    assert_eq!(total_rows, n_at * n_lag,
        "xpred nrows ({}) must equal n_at * n_lag ({})", total_rows, n_at * n_lag);
    assert_eq!(vcov.nrows(), p, "vcov nrows ({}) must equal xpred ncols ({})", vcov.nrows(), p);
    assert_eq!(vcov.ncols(), p, "vcov must be square ({}x{})", vcov.nrows(), vcov.ncols());

    let x_data = xpred.data();   // column-major: x_data[row + col * total_rows]
    let v_data = vcov.data();    // column-major: v_data[row + col * p]

    // Output: n_at rows x (n_lag + 1) columns
    // Columns 0..n_lag-1 = cumulative SE at each lag step
    // Column n_lag = overall SE (= same as last cumulative column)
    let out_cols = n_lag + 1;

    // For each "at" row, we maintain an accumulated vector xpredall_i of length p
    // and at each lag step, add the corresponding slice, then compute quadratic form.
    //
    // Parallelize across "at" rows
    let results: Vec<f64> = (0..n_at)
        .into_par_iter()
        .flat_map_iter(|i| {
            let mut xpredall = vec![0.0_f64; p];
            let mut row_results = vec![0.0_f64; out_cols];

            for (lag_step, cum_se) in row_results.iter_mut().enumerate().take(n_lag) {
                // Row index in xpred for this (at_index=i, lag_step)
                let src_row = i + n_at * lag_step;

                // Accumulate: xpredall[k] += xpred[src_row, k]
                for k in 0..p {
                    xpredall[k] += x_data[src_row + k * total_rows];
                }

                // Compute quadratic form: xpredall' V xpredall
                let mut qf = 0.0_f64;
                for j in 0..p {
                    let mut vx_j = 0.0_f64;
                    for k in 0..p {
                        vx_j += v_data[j + k * p] * xpredall[k];
                    }
                    qf += xpredall[j] * vx_j;
                }

                // pmax(0, ...) clamping and sqrt
                *cum_se = if qf < 0.0 { 0.0 } else { qf.sqrt() };
            }

            // Overall SE = final accumulated (same as last lag step)
            row_results[n_lag] = row_results[n_lag - 1];

            row_results.into_iter()
        })
        .collect();

    // results is row-major: results[i * out_cols + col]
    // Convert to column-major R matrix
    let result = RMatrix::new_matrix(n_at, out_cols, |row, col| {
        results[row * out_cols + col]
    });

    result.into_robj()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod dlnm;
    fn rust_dot_product;
    fn rust_parallel_sum;
    fn fused_crossbasis;
    fn quad_form_se;
    fn cumulative_quad_form_se;
}
