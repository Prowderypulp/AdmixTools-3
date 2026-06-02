use crate::stats::rtlchsq;
use admx_core::linalg::{dspev_l, pdinv, dgemm};
use admx_core::types::F4Info;
use crate::normab::normab;

pub fn dofrank(m: usize, n: usize, rank: usize) -> usize {
    if rank >= m || rank >= n { return 0; }
    if rank == 0 { return m * n; }
    let mut t = (m + n) * rank;
    t -= rank * rank;
    if m * n > t { m * n - t } else { 0 }
}

pub fn dofrankfix(m: usize, n: usize, rank: usize, nfix: usize) -> usize {
    if nfix > rank || rank >= m || rank >= n { return 0; }
    if rank == 0 { return m * n; }
    let mut t = m * (rank - nfix);
    t += n * rank;
    t -= rank * (rank - nfix);
    if m * n > t { m * n - t } else { 0 }
}

pub fn doranktest(mean: &[f64], var: &[f64], m: usize, n: usize, rank: usize, yscale: f64, f4pt: &mut F4Info) {
    let mut a = vec![0.0; m * rank];
    let mut b = vec![0.0; rank * n];
    
    let y = ranktest(
        mean, var, m, n, rank, yscale, 
        if rank > 0 { Some(&mut a) } else { None }, 
        if rank > 0 { Some(&mut b) } else { None }
    );
    
    let dof = dofrank(m, n, rank);
    
    f4pt.chisq = 0.0;
    if dof > 0 {
        f4pt.chisq = y;
    }
    
    f4pt.dof = dof as f64;
    f4pt.mean.copy_from_slice(mean);
    
    if rank > 0 {
        f4pt.a = a;
        f4pt.b = b;
        
        // resid = A * B
        for i in 0..m {
            for j in 0..n {
                let mut sum = 0.0;
                for r in 0..rank {
                    sum += f4pt.a[i * rank + r] * f4pt.b[r * n + j];
                }
                f4pt.resid[i * n + j] = sum;
            }
        }
    } else {
        f4pt.resid.fill(0.0);
        f4pt.a.clear();
        f4pt.b.clear();
    }
    
    // resid = mean - A * B
    for i in 0..m * n {
        f4pt.resid[i] = f4pt.mean[i] - f4pt.resid[i];
    }
}

pub fn doranktestfix(mean: &[f64], var: &[f64], m: usize, n: usize, xrank: usize, yscale: f64, f4pt: &mut F4Info, vfix: &[i32]) {
    let nfix: usize = vfix.iter().sum::<i32>() as usize;
    let rank = xrank;
    
    let mut y;
    let dof;
    let mut a = vec![0.0; m * rank];
    let mut b = vec![0.0; rank * n];
    
    if nfix == m {
        y = ranktest(mean, var, m, n, 0, yscale, None, None);
        dof = m * n;
    } else {
        y = ranktestfix(
            mean, var, m, n, rank, yscale,
            if rank > 0 { Some(&mut a) } else { None },
            if rank > 0 { Some(&mut b) } else { None },
            vfix
        );
        dof = dofrankfix(m, n, rank, nfix);
    }
    
    f4pt.chisq = 0.0;
    if dof > 0 {
        f4pt.chisq = y;
    }
    f4pt.dof = dof as f64;
    
    f4pt.mean.copy_from_slice(mean);
    if rank > 0 && nfix < m {
        f4pt.a = a;
        f4pt.b = b;
        
        for i in 0..m {
            for j in 0..n {
                let mut sum = 0.0;
                for r in 0..rank {
                    sum += f4pt.a[i * rank + r] * f4pt.b[r * n + j];
                }
                f4pt.resid[i * n + j] = sum;
            }
        }
    } else {
        f4pt.resid.fill(0.0);
        f4pt.a.clear();
        f4pt.b.clear();
    }
    
    for i in 0..m * n {
        f4pt.resid[i] = f4pt.mean[i] - f4pt.resid[i];
    }
}

/// Helper to add a scalar to the diagonal of a matrix
fn addscaldiag(mat: &mut [f64], scale: f64, n: usize) {
    let mut trace = 0.0;
    for i in 0..n {
        trace += mat[i * n + i];
    }
    let y = scale * trace;
    for i in 0..n {
        mat[i * n + i] += y;
    }
}

/// Helper to compute x^T * A * x
fn scx(a: &[f64], mean: &[f64], x: &[f64], n: usize) -> f64 {
    let mut sum = 0.0;
    for i in 0..n {
        for j in 0..n {
            sum += a[i * n + j] * (x[i] - mean[i]) * (x[j] - mean[j]);
        }
    }
    sum
}

/// Helper to solve a symmetric positive-definite system Ax = b using pdinv
fn solvit(a: &mut [f64], b: &[f64], n: usize, x: &mut [f64]) -> Result<(), i32> {
    pdinv(n, a)?;
    for i in 0..n {
        x[i] = 0.0;
        for j in 0..n {
            x[i] += a[i * n + j] * b[j];
        }
    }
    Ok(())
}

pub fn ranktest(
    mean: &[f64],
    var: &[f64],
    m: usize,
    n: usize,
    rank: usize,
    yscale: f64,
    pa: Option<&mut [f64]>,
    pb: Option<&mut [f64]>,
) -> f64 {
    let d = m * n;
    
    if rank == m {
        if let Some(pa) = pa {
            for i in 0..m {
                for j in 0..m {
                    pa[i * m + j] = if i == j { 1.0 } else { 0.0 };
                }
            }
        }
        if let Some(pb) = pb {
            pb.copy_from_slice(mean);
        }
        return 0.0;
    }

    let mut varinv = vec![0.0; d * d];
    varinv.copy_from_slice(var);
    addscaldiag(&mut varinv, yscale, d);
    if pdinv(d, &mut varinv).is_err() {
        varinv.copy_from_slice(var);
        for i in 0..d {
            for j in 0..d {
                varinv[i * d + j] = if i == j { 1.0e12 } else { 0.0 };
            }
        }
    }
    
    if rank == 0 {
        let ww = vec![0.0; d];
        return scx(&varinv, mean, &ww, d);
    }

    let adim = rank * m;
    let bdim = rank * n;
    let tdim = adim.max(bdim).max(m * n);
    
    let mut a = vec![0.0; tdim];
    let mut b = vec![0.0; tdim];
    let mut wright = vec![0.0; n * n];
    let mut evecs = vec![0.0; n * n];
    let mut mt = vec![0.0; m * n];
    let mut xmean = vec![0.0; m * n];
    let mut coeffs = vec![0.0; tdim * tdim];
    let mut rhs = vec![0.0; tdim];
    let mut ans = vec![0.0; tdim];
    let mut ww = vec![0.0; n];

    // transpose mean -> mt
    for i in 0..m {
        for j in 0..n {
            mt[j * m + i] = mean[i * n + j];
        }
    }

    // wright = mt * mean (n x m) * (m x n) -> n x n
    for i in 0..n {
        for j in 0..n {
            let mut sum = 0.0;
            for k in 0..m {
                sum += mt[i * m + k] * mean[k * n + j];
            }
            wright[i * n + j] = sum;
        }
    }

    // eigenvalue decomposition
    let mut wright_packed = vec![0.0; n * (n + 1) / 2];
    let mut idx = 0;
    for j in 0..n {
        for i in j..n {
            wright_packed[idx] = -wright[i * n + j];
            idx += 1;
        }
    }
    dspev_l(n, &mut wright_packed, &mut ww, &mut evecs).expect("dspev failed");
    for i in 0..n {
        ww[i] = -ww[i];
    }
    
    // LAPACK dspev stores eigenvector i contiguously (column i of Z), so
    // eigenvector i is evecs[i*n + j]. C reads it the same way
    // (f4rank.c:555, copyarr(evecs + i*n, ...)). Reading evecs[j*n+i] would
    // transpose the eigenvector matrix and seed the ALS with the wrong B.
    for i in 0..rank {
        for j in 0..n {
            b[i * n + j] = evecs[i * n + j];
        }
    }

    for i in 0..m {
        for r in 0..rank {
            let mut sum = 0.0;
            for j in 0..n {
                sum += mean[i * n + j] * b[r * n + j];
            }
            a[i * rank + r] = sum;
        }
    }

    normab(&mut a, &mut b, m, n, rank);

    for i in 0..m {
        for j in 0..n {
            let mut sum = 0.0;
            for r in 0..rank {
                sum += a[i * rank + r] * b[r * n + j];
            }
            xmean[i * n + j] = sum;
        }
    }

    let mut y1 = scx(&varinv, mean, &xmean, d);
    let numiter = 20;

    for _iter in 1..=numiter {
        coeffs.fill(0.0);
        rhs.fill(0.0);

        for i in 0..m {
            for s in 0..n {
                let a_idx = i * n + s;
                for j in 0..m {
                    for t in 0..n {
                        let b_idx = j * n + t;
                        let y = varinv[a_idx * d + b_idx];
                        for r1 in 0..rank {
                            let k1 = i * rank + r1;
                            let u1 = r1 * n + s;
                            rhs[u1] += y * a[k1] * mean[b_idx];
                            for r2 in 0..rank {
                                let k2 = j * rank + r2;
                                let u2 = r2 * n + t;
                                coeffs[u1 * bdim + u2] += y * a[k1] * a[k2];
                            }
                        }
                    }
                }
            }
        }

        addscaldiag(&mut coeffs, yscale, bdim);
        solvit(&mut coeffs, &rhs, bdim, &mut ans).expect("bad ranktest B");

        b[..bdim].copy_from_slice(&ans[..bdim]);
        normab(&mut a, &mut b, m, n, rank);

        coeffs.fill(0.0);
        rhs.fill(0.0);

        for i in 0..m {
            for s in 0..n {
                let a_idx = i * n + s;
                for j in 0..m {
                    for t in 0..n {
                        let b_idx = j * n + t;
                        let y = varinv[a_idx * d + b_idx];
                        for r1 in 0..rank {
                            let k1 = i * rank + r1;
                            let u1 = r1 * n + s;
                            rhs[k1] += y * b[u1] * mean[b_idx];
                            for r2 in 0..rank {
                                let k2 = j * rank + r2;
                                let u2 = r2 * n + t;
                                coeffs[k1 * adim + k2] += y * b[u1] * b[u2];
                            }
                        }
                    }
                }
            }
        }

        ans.fill(0.0);
        addscaldiag(&mut coeffs, yscale, adim);
        solvit(&mut coeffs, &rhs, adim, &mut ans).expect("bad ranktest A");

        a[..adim].copy_from_slice(&ans[..adim]);
        normab(&mut a, &mut b, m, n, rank);

        for i in 0..m {
            for j in 0..n {
                let mut sum = 0.0;
                for r in 0..rank {
                    sum += a[i * rank + r] * b[r * n + j];
                }
                xmean[i * n + j] = sum;
            }
        }

        y1 = scx(&varinv, mean, &xmean, d);
    }

    if let Some(pa) = pa {
        pa.copy_from_slice(&a[..adim]);
    }
    if let Some(pb) = pb {
        pb.copy_from_slice(&b[..bdim]);
    }

    y1
}

fn solvitforcez(coeffs: &mut [f64], rhs: &mut [f64], dim: usize, ans: &mut [f64], vl: &[usize]) -> Result<(), i32> {
    // `coeffs`/`rhs` are the caller's oversized work buffers (tdim²/tdim); the
    // meaningful system is the leading dim×dim block (row stride dim) and the
    // first dim entries of rhs, matching how solvit reads them.
    let mut tco = vec![0.0; dim * dim];
    let mut trhs = vec![0.0; dim];
    tco.copy_from_slice(&coeffs[..dim * dim]);
    trhs.copy_from_slice(&rhs[..dim]);

    for &v in vl {
        for k in 0..dim {
            trhs[k] -= tco[k * dim + v] * 0.0; // vval is 0.0
            tco[k * dim + v] = 0.0;
            tco[v * dim + k] = 0.0;
        }
        tco[v * dim + v] = 1.0;
        trhs[v] = 0.0;
    }

    solvit(&mut tco, &trhs, dim, ans)
}

pub fn ranktestfix(
    mean: &[f64],
    var: &[f64],
    m: usize,
    n: usize,
    rank: usize,
    yscale: f64,
    pa: Option<&mut [f64]>,
    pb: Option<&mut [f64]>,
    vfix: &[i32],
) -> f64 {
    let d = m * n;
    
    let nfix: i32 = vfix.iter().sum();
    if nfix as usize > rank {
        panic!("(ranktestfix) too many fixed variables");
    }

    let mut vl = Vec::with_capacity(m * rank);
    let mut k_idx = 0;
    for i in 0..m {
        if vfix[i] == 0 { continue; }
        for j in 0..m {
            if i == j { continue; }
            vl.push(j * rank + k_idx);
        }
        k_idx += 1;
    }

    let mut varinv = vec![0.0; d * d];
    varinv.copy_from_slice(var);
    addscaldiag(&mut varinv, yscale, d);
    if pdinv(d, &mut varinv).is_err() {
        varinv.copy_from_slice(var);
        for i in 0..d {
            for j in 0..d {
                varinv[i * d + j] = if i == j { 1.0e12 } else { 0.0 };
            }
        }
    }

    let adim = rank * m;
    let bdim = rank * n;
    let tdim = adim.max(bdim).max(m * n);
    
    let mut a = vec![0.0; tdim];
    let mut b = vec![0.0; tdim];
    let mut xmean = vec![0.0; m * n];
    let mut coeffs = vec![0.0; tdim * tdim];
    let mut rhs = vec![0.0; tdim];
    let mut ans = vec![0.0; tdim];

    // C gaussa initializes with standard normal. 
    // Here we use a pseudo-random approach or ones just to be deterministic.
    // Given iter goes to 50, standard normal is better. Let's use basic deterministic LCG
    let mut seed = 123456789u32;
    let mut next_rand = || -> f64 {
        seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
        (seed as f64 / (std::u32::MAX as f64)) * 2.0 - 1.0
    };

    let mut iset = false;
    let mut gset = 0.0;
    let mut gauss = || -> f64 {
        if iset {
            iset = false;
            return gset;
        }
        loop {
            let v1 = next_rand();
            let v2 = next_rand();
            let rsq = v1 * v1 + v2 * v2;
            if rsq < 1.0 && rsq > 0.0 {
                let fac = (-2.0 * rsq.ln() / rsq).sqrt();
                gset = v1 * fac;
                iset = true;
                return v2 * fac;
            }
        }
    };

    for i in 0..adim { a[i] = gauss(); }
    for i in 0..bdim { b[i] = gauss(); }
    
    normab(&mut a, &mut b, m, n, rank);

    for i in 0..m {
        for j in 0..n {
            let mut sum = 0.0;
            for r in 0..rank {
                sum += a[i * rank + r] * b[r * n + j];
            }
            xmean[i * n + j] = sum;
        }
    }

    let mut y1 = scx(&varinv, mean, &xmean, d);
    let numiter = 50;

    for _iter in 1..=numiter {
        coeffs.fill(0.0);
        rhs.fill(0.0);

        for i in 0..m {
            for s in 0..n {
                let a_idx = i * n + s;
                for j in 0..m {
                    for t in 0..n {
                        let b_idx = j * n + t;
                        let y = varinv[a_idx * d + b_idx];
                        for r1 in 0..rank {
                            let k1 = i * rank + r1;
                            let u1 = r1 * n + s;
                            rhs[u1] += y * a[k1] * mean[b_idx];
                            for r2 in 0..rank {
                                let k2 = j * rank + r2;
                                let u2 = r2 * n + t;
                                coeffs[u1 * bdim + u2] += y * a[k1] * a[k2];
                            }
                        }
                    }
                }
            }
        }

        addscaldiag(&mut coeffs, yscale, bdim);
        solvit(&mut coeffs, &rhs, bdim, &mut ans).expect("bad ranktestfix B");

        b[..bdim].copy_from_slice(&ans[..bdim]);
        normab(&mut a, &mut b, m, n, rank);

        coeffs.fill(0.0);
        rhs.fill(0.0);

        for i in 0..m {
            for s in 0..n {
                let a_idx = i * n + s;
                for j in 0..m {
                    for t in 0..n {
                        let b_idx = j * n + t;
                        let y = varinv[a_idx * d + b_idx];
                        for r1 in 0..rank {
                            let k1 = i * rank + r1;
                            let u1 = r1 * n + s;
                            rhs[k1] += y * b[u1] * mean[b_idx];
                            for r2 in 0..rank {
                                let k2 = j * rank + r2;
                                let u2 = r2 * n + t;
                                coeffs[k1 * adim + k2] += y * b[u1] * b[u2];
                            }
                        }
                    }
                }
            }
        }

        ans.fill(0.0);
        addscaldiag(&mut coeffs, yscale, adim);
        solvitforcez(&mut coeffs, &mut rhs, adim, &mut ans, &vl).expect("bad ranktestfix A");

        a[..adim].copy_from_slice(&ans[..adim]);
        normab(&mut a, &mut b, m, n, rank);

        for i in 0..m {
            for j in 0..n {
                let mut sum = 0.0;
                for r in 0..rank {
                    sum += a[i * rank + r] * b[r * n + j];
                }
                xmean[i * n + j] = sum;
            }
        }

        y1 = scx(&varinv, mean, &xmean, d);
    }

    if let Some(pa) = pa {
        pa.copy_from_slice(&a[..adim]);
    }
    if let Some(pb) = pb {
        pb.copy_from_slice(&b[..bdim]);
    }

    y1
}
