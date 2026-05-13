use admx_rank::ranktest::doranktest;
use admx_core::types::F4Info;
use crate::calcadm::calcadm;
use crate::prng::LegacyLcg;
use rayon::prelude::*;
use std::sync::OnceLock;

fn bootstrap_pool() -> &'static rayon::ThreadPool {
    static POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();
    POOL.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .start_handler(|_| admx_core::linalg::set_blas_single_threaded())
            .build()
            .expect("failed to build bootstrap rayon pool")
    })
}

fn cholesky_lower_rowmajor(covar: &[f64], n: usize) -> Result<Vec<f64>, i32> {
    let mut l = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let mut s = covar[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k];
            }
            if i == j {
                if s <= 0.0 || !s.is_finite() {
                    return Err(-1);
                }
                l[i * n + j] = s.sqrt();
            } else {
                let d = l[j * n + j];
                if d == 0.0 || !d.is_finite() {
                    return Err(-1);
                }
                l[i * n + j] = s / d;
            }
        }
    }
    Ok(l)
}

/// Generate independent N(0,1) samples using Box-Muller on the legacy PRNG.
///
/// This must match C `gauss()` which uses `DRAND2()` (not `DRAND()`).
fn gaussa(x: &mut [f64], rng: &mut LegacyLcg) {
    let mut iset = 0;
    let mut gset = 0.0;
    for i in 0..x.len() {
        if iset == 0 {
            let mut rsq;
            let mut v1;
            let mut v2;
            loop {
                v1 = 2.0 * rng.next_f64() - 1.0;
                v2 = 2.0 * rng.next_f64() - 1.0;
                rsq = v1 * v1 + v2 * v2;
                if rsq < 1.0 && rsq != 0.0 {
                    break;
                }
            }
            let fac = (-2.0 * rsq.ln() / rsq).sqrt();
            gset = v1 * fac;
            iset = 1;
            x[i] = v2 * fac;
        } else {
            iset = 0;
            x[i] = gset;
        }
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
    // C `genmultgauss` uses row-major `cholesky` followed by transpose.
    // If L is lower-triangular (row-major), multiplying standard normals by L^T
    // yields covariance yvar.
    let cf = cholesky_lower_rowmajor(yvar, dim)?;

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

    // Worker: process one bootstrap iteration `k` (both passes), writing
    // its outputs to the slot `tmean[2*k*nl..(2*k+2)*nl]`. `ranktest` is
    // deterministic — the only PRNG in admx-rank lives in `ranktestfix`,
    // which is not on this path — so iterations are pure functions of
    // their `rvec` slices and may run in any order.
    let make_f4info = || F4Info {
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

    let process_k = |k: usize, slot: &mut [f64]| {
        let mut f4info = make_f4info();
        for pass in 0..2 {
            let sample_mean = if pass == 0 {
                &rvec[k * dim..(k + 1) * dim]
            } else {
                &rvec[(numboot + k) * dim..(numboot + k + 1) * dim]
            };
            doranktest(sample_mean, yvar, nl, nr, nl - 1, 0.0001, &mut f4info);
            let mut ans = vec![0.0_f64; nl];
            if calcadm(&mut ans, &f4info.a, nl).is_ok() {
                slot[pass * nl..(pass + 1) * nl].copy_from_slice(&ans);
            }
        }
    };

    // k=0 runs first (serially) to emit the legacy `zzevarboot`/`zzevarboot2`
    // diagnostic lines in deterministic order, then we parallelize the rest.
    if numboot > 0 {
        let (head, tail) = tmean.split_at_mut(2 * nl);
        process_k(0, head);
        println!("zzevarboot");
        for i in 0..nl { print!("{:15.9} ", head[i]); }
        println!();
        println!("zzevarboot2");
        for i in 0..nl { print!("{:15.9} ", head[nl + i]); }
        println!();

        // Use a dedicated pool whose workers force single-threaded BLAS, so
        // each rayon worker that calls into OpenBLAS does not fan out across
        // libgomp again — that nested threading is what made the naive
        // par_iter ~50× slower than the serial loop. The pool is shared via
        // a process-global OnceLock so we don't pay setup cost per qpAdm call.
        bootstrap_pool().install(|| {
            tail.par_chunks_mut(2 * nl)
                .enumerate()
                .for_each(|(idx, slot)| {
                    process_k(idx + 1, slot);
                });
        });
    }

    let n_samples = (2 * numboot) as f64;
    for i in 0..nl {
        let mut sum = 0.0;
        for b in 0..(2 * numboot) {
            sum += tmean[b * nl + i];
        }
        jmean[i] = sum / n_samples;
    }

    let cov_denom = (2 * numboot - 1) as f64;
    for i in 0..nl {
        for j in 0..nl {
            let mut sum = 0.0;
            for b in 0..(2 * numboot) {
                let diff_i = tmean[b * nl + i] - jmean[i];
                let diff_j = tmean[b * nl + j] - jmean[j];
                sum += diff_i * diff_j;
            }
            var[i * nl + j] = sum / cov_denom;
        }
    }

    Ok(())
}
