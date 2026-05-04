//! Thin wrappers around LAPACK/BLAS routines used by AdmixTools.
//!
//! The C code calls `dspev`, `dpotrf`, `dgetrf`, `dgetri`, `dgetrs`, `dsygv`
//! via a mix of direct calls and `eigx.c` wrappers.  We wrap the same
//! LAPACK entry points here so the numerical path is identical.

/// Compute eigenvalues and eigenvectors of a real symmetric matrix
/// stored in packed upper-triangular form.
///
/// Wraps LAPACK `dspev`.
///
/// # Arguments
/// * `n` — order of the matrix
/// * `ap` — packed upper triangle, length `n*(n+1)/2` (overwritten)
/// * `w` — eigenvalues output, length `n`
/// * `z` — eigenvectors output, column-major `n × n`
///
/// # Returns
/// `Ok(())` on success, `Err` with LAPACK info code on failure.
pub fn dspev(n: usize, ap: &mut [f64], w: &mut [f64], z: &mut [f64]) -> Result<(), i32> {
    dspev_uplo(b'U', n, ap, w, z)
}

pub fn dspev_l(n: usize, ap: &mut [f64], w: &mut [f64], z: &mut [f64]) -> Result<(), i32> {
    dspev_uplo(b'L', n, ap, w, z)
}

fn dspev_uplo(uplo: u8, n: usize, ap: &mut [f64], w: &mut [f64], z: &mut [f64]) -> Result<(), i32> {
    let n_i32 = n as i32;
    let mut work = vec![0.0_f64; 3 * n];
    let mut info: i32 = 0;

    unsafe {
        lapack::dspev(
            b'V',  // compute eigenvalues and eigenvectors
            uplo,  // upper or lower triangle
            n_i32,
            ap,
            w,
            z,
            n_i32, // ldz
            &mut work,
            &mut info,
        );
    }

    if info == 0 {
        Ok(())
    } else {
        Err(info)
    }
}

/// Cholesky factorization of a symmetric positive-definite matrix.
///
/// Wraps LAPACK `dpotrf`.  On return, `a` contains the upper Cholesky factor.
pub fn dpotrf(n: usize, a: &mut [f64]) -> Result<(), i32> {
    let n_i32 = n as i32;
    let mut info: i32 = 0;

    unsafe {
        lapack::dpotrf(b'U', n_i32, a, n_i32, &mut info);
    }

    if info == 0 {
        Ok(())
    } else {
        Err(info)
    }
}

/// Inverse of a symmetric positive-definite matrix via Cholesky.
///
/// Wraps LAPACK `dpotrf` followed by `dpotri`.  Input/output is stored in
/// column-major layout.  On success the result is the full symmetric inverse
/// (both upper and lower triangles populated).
///
/// Mirrors `pdinv` in `nicksrc/linsubs.c`.
pub fn pdinv(n: usize, a: &mut [f64]) -> Result<(), i32> {
    let n_i32 = n as i32;
    let mut info: i32 = 0;

    unsafe {
        lapack::dpotrf(b'U', n_i32, a, n_i32, &mut info);
        if info != 0 { return Err(info); }
        lapack::dpotri(b'U', n_i32, a, n_i32, &mut info);
        if info != 0 { return Err(info); }
    }

    // LAPACK (column-major) with uplo='U' populates elements A(i,j) for i<=j
    // at `a[j*n + i]`. Mirror those into the unpopulated lower triangle at
    // `a[i*n + j]` so the caller can treat the buffer as a full symmetric
    // matrix in either layout.
    for j in 0..n {
        for i in 0..j {
            a[i * n + j] = a[j * n + i];
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use rand::Rng;

    #[test]
    fn test_pdinv_random_spd() {
        let n = 10;
        let mut rng = rand::thread_rng();
        
        // Generate a random matrix A
        let mut a_data = vec![0.0; n * n];
        for val in a_data.iter_mut() {
            *val = rng.gen_range(-1.0..1.0);
        }
        let a = Array2::from_shape_vec((n, n), a_data).unwrap();
        
        // S = A.T * A + n * I (ensure it's well-conditioned SPD)
        let mut s = a.t().dot(&a);
        for i in 0..n {
            s[[i, i]] += n as f64;
        }
        
        let mut s_buf = s.clone().into_raw_vec();
        pdinv(n, &mut s_buf).expect("pdinv failed");
        
        let sinv = Array2::from_shape_vec((n, n), s_buf).unwrap();
        
        // Check S * Sinv \approx I
        let identity = s.dot(&sinv);
        for j in 0..n {
            for i in 0..n {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((identity[[i, j]] - expected).abs() < 1e-10, 
                        "At [{}, {}]: expected {}, got {}", i, j, expected, identity[[i, j]]);
            }
        }
    }
}

/// Rank-1 update of a symmetric matrix: A := alpha * x * x^T + A.
///
/// Wraps BLAS `dsyr`.  `a` is stored in column-major upper triangle.
#[inline]
pub fn dsyr(n: usize, alpha: f64, x: &[f64], a: &mut [f64]) {
    let n_i32 = n as i32;
    unsafe {
        blas::dsyr(b'U', n_i32, alpha, x, 1, a, n_i32);
    }
}

/// General matrix-matrix multiply: C := alpha * A * B + beta * C.
///
/// Wraps BLAS `dgemm`.
#[allow(clippy::too_many_arguments)]
pub fn dgemm(
    transa: u8,
    transb: u8,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    beta: f64,
    c: &mut [f64],
    ldc: usize,
) {
    unsafe {
        blas::dgemm(
            transa,
            transb,
            m as i32,
            n as i32,
            k as i32,
            alpha,
            a,
            lda as i32,
            b,
            ldb as i32,
            beta,
            c,
            ldc as i32,
        );
    }
}

/// Solves a general system of linear equations A * X = B.
///
/// Wraps LAPACK `dgesv`.
/// Note: LAPACK expects column-major layout. The matrix A will be overwritten with the LU factors.
/// On return, `b` contains the solution vector.
///
/// Returns `Ok(())` on success, or `Err(info)` if the matrix is exactly singular.
pub fn linsolv(n: usize, a: &mut [f64], b: &mut [f64]) -> Result<(), i32> {
    let n_i32 = n as i32;
    let mut ipiv = vec![0_i32; n];
    let mut info = 0_i32;

    unsafe {
        lapack::dgesv(n_i32, 1, a, n_i32, &mut ipiv, b, n_i32, &mut info);
    }

    if info == 0 {
        Ok(())
    } else {
        Err(info)
    }
}
