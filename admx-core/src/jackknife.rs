//! Weighted block-jackknife estimators.
//!
//! Ports of `wjackvest` and `weightjack` from `qpsubs.c`.
//! `weightjackfourier` is explicitly out of scope (rolloff only).

use crate::error::AdmxResult;

/// Weighted jackknife mean and variance estimator for a scalar quantity.
///
/// Given a global mean and per-block leave-one-out means with block weights,
/// computes the jackknife estimate and standard error.
///
/// # Arguments
/// * `mean` — full-sample estimate
/// * `jmean` — leave-one-out estimates, length `g`
/// * `jwt` — block weights (number of SNPs in each block), length `g`
///
/// # Returns
/// `(estimate, standard_error)`
///
/// Mirrors `wjackest` in `qpsubs.c`.
pub fn wjackest(mean: f64, jmean: &[f64], jwt: &[f64]) -> (f64, f64) {
    let g = jmean.len();
    debug_assert_eq!(jwt.len(), g);

    if g == 0 {
        return (mean, 0.0);
    }

    let wtot: f64 = jwt.iter().sum();
    if wtot <= 0.0 {
        return (mean, 0.0);
    }

    let mut jackest = 0.0;
    let mut sum_jwt_jmean = 0.0;
    for k in 0..g {
        jackest += mean - jmean[k];
        sum_jwt_jmean += jmean[k] * jwt[k];
    }
    jackest += sum_jwt_jmean / wtot;

    let gf = g as f64;
    let mut var = 0.0;
    for k in 0..g {
        let h = wtot / jwt[k];
        let xtau = h * mean - (h - 1.0) * jmean[k] - jackest;
        var += xtau * xtau / (h - 1.0);
    }
    var /= gf;

    (jackest, var.sqrt())
}

/// Weighted jackknife variance estimator for a vector quantity.
///
/// Given per-block leave-one-out mean vectors and block weights,
/// computes the jackknife variance-covariance matrix.
///
/// # Arguments
/// * `d` — dimension of the vector quantity
/// * `mean` — full-sample mean vector, length `d`
/// * `jmean` — leave-one-out mean vectors, shape `[g][d]` (row-major)
/// * `jwt` — block weights, length `g`
///
/// # Returns
/// `(jackest, var)` — the jackknife-corrected mean vector of length `d`,
/// and the variance-covariance matrix as a flat `d × d` array (row-major).
///
/// Mirrors `wjackvest` in `qpsubs.c`.
pub fn wjackvest(d: usize, mean: &[f64], jmean: &[Vec<f64>], jwt: &[f64]) -> AdmxResult<(Vec<f64>, Vec<f64>)> {
    let g = jmean.len();
    debug_assert_eq!(jwt.len(), g);
    debug_assert_eq!(mean.len(), d);

    let mut jackest = vec![0.0; d];
    let mut var = vec![0.0; d * d];

    if g == 0 {
        jackest.copy_from_slice(mean);
        return Ok((jackest, var));
    }

    let wtot: f64 = jwt.iter().sum();
    if wtot <= 0.0 {
        jackest.copy_from_slice(mean);
        return Ok((jackest, var));
    }

    let gf = g as f64;

    // Jackknife estimate: jackest = Σ_k (mean - jmean[k]) + Σ_k jmean[k] * jwt[k] / wtot
    for i in 0..d {
        let mut sum_diff = 0.0;
        let mut sum_jwt_jmean = 0.0;
        for k in 0..g {
            sum_diff += mean[i] - jmean[k][i];
            sum_jwt_jmean += jmean[k][i] * jwt[k];
        }
        jackest[i] = sum_diff + sum_jwt_jmean / wtot;
    }

    // Accumulate covariance
    for k in 0..g {
        let h = wtot / jwt[k];
        let mut xtau = vec![0.0; d];
        for i in 0..d {
            xtau[i] = h * mean[i] - (h - 1.0) * jmean[k][i] - jackest[i];
        }

        for i in 0..d {
            for j in 0..d {
                var[i * d + j] += xtau[i] * xtau[j] / (h - 1.0);
            }
        }
    }

    // Normalize
    for v in var.iter_mut() {
        *v /= gf;
    }

    Ok((jackest, var))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wjackest_trivial() {
        // All blocks identical ⇒ zero variance.
        let mean = 1.0;
        let jmean = vec![1.0, 1.0, 1.0];
        let jwt = vec![1.0, 1.0, 1.0];
        let (est, sig) = wjackest(mean, &jmean, &jwt);
        assert!((est - 1.0).abs() < 1e-15);
        assert!(sig.abs() < 1e-15);
    }

    #[test]
    fn test_wjackest_empty() {
        let (est, sig) = wjackest(42.0, &[], &[]);
        assert!((est - 42.0).abs() < 1e-15);
        assert!(sig.abs() < 1e-15);
    }
}
