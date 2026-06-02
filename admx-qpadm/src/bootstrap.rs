use admx_core::linalg::{dpotrf, dgemm};
use admx_rank::ranktest::doranktest;
use admx_core::types::F4Info;
use crate::calcadm::calcadm;
use crate::prng::LegacyLcg;

/// Generate independent N(0,1) samples using Box-Muller on the legacy PRNG.
fn gaussa(x: &mut [f64], rng: &mut LegacyLcg) {
    let mut i = 0;
    while i < x.len() {
        let u1 = rng.next_f64();
        let u2 = rng.next_f64();
        
        let r = (-2.0 * (1.0 - u1).ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;
        
        x[i] = r * theta.cos();
        if i + 1 < x.len() {
            x[i + 1] = r * theta.sin();
        }
        i += 2;
    }
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
    let mut cf = yvar.to_vec();
    dpotrf(dim, &mut cf)?;

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
