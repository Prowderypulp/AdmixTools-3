use admx_rank::ranktest::doranktest;
use admx_core::types::F4Info;
use crate::calcadm::calcadm;
use crate::prng::LegacyLcg;

/// Fill `x` with standard normals, exactly as C `gaussa` (`nicksrc/gauss.c`):
/// `for i: x[i] = gauss()`.  Element order and the per-element draw via the
/// cached Marsaglia-polar `gauss()` must match C so the consumed-uniform stream
/// stays aligned.
fn gaussa(x: &mut [f64], rng: &mut LegacyLcg) {
    for xi in x.iter_mut() {
        *xi = rng.gauss();
    }
}

/// Lower-triangular Cholesky factor `L` of the `dim`×`dim` SPD `a` (row-major),
/// returning `L` row-major (`L[j*dim+k]`, nonzero for `k <= j`).
///
/// Exact port of C `choldc`/`cholesky` (`nicksrc/linsubs.c`), NOT LAPACK
/// `dpotrf` — the Numerical-Recipes algorithm with its **descending** inner
/// accumulation (`k = i-1 .. 0`) is required for bit-identical rounding against
/// the legacy `genmultgauss`.  Returns `Err` if `a` is not positive definite.
fn choldc_lower(a: &[f64], dim: usize) -> Result<Vec<f64>, i32> {
    // Work in C's `choldc` layout: diagonal in `p`, sub-diagonal stored in the
    // lower triangle of `t` (a scratch copy of `a`).
    let mut t = a.to_vec();
    let mut p = vec![0.0_f64; dim];
    for i in 0..dim {
        for j in i..dim {
            let mut sum = t[i * dim + j];
            // Descending accumulation, matching C `for (k=i-1; k>=0; --k)`.
            for k in (0..i).rev() {
                sum -= t[i * dim + k] * t[j * dim + k];
            }
            if i == j {
                if sum <= 0.0 {
                    return Err(-1); // not positive definite
                }
                p[i] = sum.sqrt();
            } else {
                t[j * dim + i] = sum / p[i];
            }
        }
    }
    // Assemble L row-major: L[i][i] = p[i], L[i][j<i] = t[i][j], else 0.
    let mut l = vec![0.0_f64; dim * dim];
    for i in 0..dim {
        l[i * dim + i] = p[i];
        for j in 0..i {
            l[i * dim + j] = t[i * dim + j];
        }
    }
    Ok(l)
}

/// Computes the bootstrap covariance of the admixture weights.
/// This matches `calcevarboot` in `qpAdm.c`.
pub fn calcevarboot(
    jmean: &mut [f64],
    var: &mut [f64],
    ymean: &[f64],
    yvar: &[f64],
    nl: usize,
    nr: usize,
    dim: usize,
    numboot: usize,
    rng: &mut LegacyLcg,
) -> Result<(), i32> {
    // C `genmultgauss` (gds.c): cholesky → transpose → gaussa → mulmat.  The
    // transpose+mulmat composition reduces to `rvec[b][j] = sum_{k<=j} g[b][k]
    // * L[j][k]` with L the *lower* factor (see the noise loop below), so we
    // only need L itself, computed by the legacy `choldc`.
    let cf = choldc_lower(yvar, dim)?;

    let mut x = vec![0.0_f64; numboot * dim];
    gaussa(&mut x, rng);

    // rvec will store both +noise and -noise samples
    let mut rvec = vec![0.0_f64; 2 * numboot * dim];
    for b in 0..numboot {
        for j in 0..dim {
            let mut sum = 0.0;
            for k in 0..=j {
                sum += x[b * dim + k] * cf[j * dim + k];
            }
            let noise = sum;
            rvec[b * dim + j] = ymean[j] + noise;
            rvec[(numboot + b) * dim + j] = ymean[j] - noise;
        }
    }

    // Per-sample rank test + admixture solve. NOTE: this loop is kept SERIAL on
    // purpose. The 2*numboot samples are independent and would parallelize
    // cleanly, but `doranktest` calls into OpenBLAS (dspev/dgemm/pdinv), and the
    // system OpenBLAS build is not safe to *enter* concurrently from multiple
    // application threads — a rayon parallel version deadlocked. Each call here
    // still uses multi-threaded BLAS internally, which is safe and gives good
    // throughput on the ~55×55 matrices.
    let mut tmean = vec![0.0_f64; 2 * numboot * nl];
    for b in 0..(2 * numboot) {
        let sample_mean = &rvec[b * dim .. (b + 1) * dim];

        let mut f4info = F4Info {
            nl,
            nr,
            rank: nl - 1,
            dof_jack: 0.0,
            dof: 0.0,
            dof_diff: 0.0,
            chisq: 0.0,
            chisq_diff: 0.0,
            a: vec![0.0; nl * (nl - 1)],
            b: vec![0.0; nr * (nl - 1)],
            mean: vec![0.0; nl * nr],
            resid: vec![0.0; nl * nr],
        };

        doranktest(sample_mean, yvar, nl, nr, nl - 1, 0.0001, &mut f4info);

        let mut ans = vec![0.0_f64; nl];
        if let Ok(()) = calcadm(&mut ans, &f4info.a, nl) {
            for i in 0..nl {
                tmean[b * nl + i] = ans[i];
            }
        }
    }

    let n_samples = (2 * numboot) as f64;
    for i in 0..nl {
        let mut sum = 0.0;
        for b in 0..(2 * numboot) {
            sum += tmean[b * nl + i];
        }
        jmean[i] = sum / n_samples;
    }

    for i in 0..nl {
        for j in 0..nl {
            let mut sum = 0.0;
            for b in 0..(2 * numboot) {
                let diff_i = tmean[b * nl + i] - jmean[i];
                let diff_j = tmean[b * nl + j] - jmean[j];
                sum += diff_i * diff_j;
            }
            var[i * nl + j] = sum / n_samples;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values produced by the legacy C library (`nicksrc/libnick.a`).
    // Reproduce by linking this against libnick (`gcc t.c -Inicksrc
    // nicksrc/libnick.a -lm -llapack -lblas`):
    //
    //   #include <stdio.h>
    //   #include "ranmath.h"
    //   void gaussa(double*, int); void genmultgauss(double*, int, int, double*);
    //   int main(){
    //     double g[8]; SRAND(1923698036); gaussa(g,8);
    //     for(int i=0;i<8;i++) printf("%.17g\n", g[i]);
    //     double covar[9]={4,1,0,1,3,1,0,1,2}, r[6];
    //     SRAND(1923698036); genmultgauss(r,2,3,covar);
    //     for(int i=0;i<6;i++) printf("%.17g\n", r[i]); return 0; }
    //
    // This pins the gauss() + choldc + multiply chain bit-for-bit given identical
    // inputs.  It does NOT exercise the production `yvar`, which is computed
    // upstream and may itself differ from C at last digits.
    #[test]
    fn gauss_matches_legacy_c() {
        let _guard = crate::prng::RNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let expected: [f64; 8] = [
            -0.048790706840763566,
            1.4638641257289076,
            1.2976350718888496,
            -0.079583705641651098,
            0.6306728841276914,
            0.054550983893108897,
            2.7680219923793556,
            -0.41270728870127615,
        ];
        let mut rng = LegacyLcg::new(1923698036);
        for (i, &e) in expected.iter().enumerate() {
            let got = rng.gauss();
            assert_eq!(got.to_bits(), e.to_bits(), "gauss[{i}]: {got} != {e}");
        }
    }

    #[test]
    fn genmultgauss_matches_legacy_c() {
        let _guard = crate::prng::RNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // C `genmultgauss(rvec, num=2, n=3, covar)` == gaussa(2*3) then
        // rvec[b][j] = sum_{k<=j} g[b*3+k] * L[j][k] with L = choldc_lower(covar).
        let covar = [4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0];
        let expected: [f64; 6] = [
            -0.097581413681527132,
            2.4031486711318326,
            2.5426836430245512,
            -0.1591674112833022,
            1.0060608082305942,
            0.45009191162763784,
        ];
        let dim = 3;
        let num = 2;
        let l = choldc_lower(&covar, dim).unwrap();
        let mut rng = LegacyLcg::new(1923698036);
        let mut g = vec![0.0_f64; num * dim];
        gaussa(&mut g, &mut rng);
        let mut rvec = vec![0.0_f64; num * dim];
        for b in 0..num {
            for j in 0..dim {
                let mut sum = 0.0;
                for k in 0..=j {
                    sum += g[b * dim + k] * l[j * dim + k];
                }
                rvec[b * dim + j] = sum;
            }
        }
        for (i, &e) in expected.iter().enumerate() {
            assert_eq!(rvec[i].to_bits(), e.to_bits(), "rvec[{i}]: {} != {e}", rvec[i]);
        }
    }
}
