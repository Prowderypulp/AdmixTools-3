use admx_core::linalg::linsolv;

/// Weight solver for qpAdm.
/// `ans` gets length `n`. `A` is `n x (n-1)` row-major matrix from `doranktest`.
pub fn calcadm(ans: &mut [f64], a: &[f64], n: usize) -> Result<(), i32> {
    if n == 1 {
        ans[0] = 1.0;
        return Ok(());
    }

    // Build coeff matrix of size n x n in column-major layout for LAPACK
    let mut coeff = vec![0.0_f64; n * n];
    for j in 0..n {
        for i in 0..(n - 1) {
            // coeff(i, j) = A(j, i). In column-major, this is index j * n + i.
            // A is row-major n x (n-1). A(j, i) is a[j * (n - 1) + i].
            coeff[j * n + i] = a[j * (n - 1) + i];
        }
        // Last row is 1.0. coeff(n-1, j)
        coeff[j * n + n - 1] = 1.0;
    }

    let mut rhs = vec![0.0_f64; n];
    rhs[n - 1] = 1.0;

    let mut baditer = 0;
    let mut success = false;

    // trace is the sum of diagonal elements of coeff.
    let mut ytrace = 0.0;
    for i in 0..n {
        ytrace += coeff[i * n + i];
    }

    loop {
        if baditer >= 10 {
            break;
        }

        let mut coeff_copy = coeff.clone();
        let mut rhs_copy = rhs.clone();

        match linsolv(n, &mut coeff_copy, &mut rhs_copy) {
            Ok(()) => {
                let sum: f64 = rhs_copy.iter().sum();
                if sum.is_finite() {
                    ans.copy_from_slice(&rhs_copy);
                    success = true;
                    break;
                }
            }
            Err(_) => {}
        }

        baditer += 1;
        for i in 0..n {
            coeff[i * n + i] += ytrace * 0.001;
        }
    }

    if success {
        let sum: f64 = ans.iter().sum();
        if sum > 0.0 {
            for v in ans.iter_mut() {
                *v /= sum;
            }
        }
        Ok(())
    } else {
        Err(-1) // Failed to solve
    }
}

/// Weight solver for qpAdm with fixed parameters.
/// `vf` is slice of size `n`, where `vf[i] == 1` means `ans[i]` is fixed to 0.
pub fn calcadmfix(ans: &mut [f64], a: &[f64], n: usize, vf: &[i32]) -> Result<(), i32> {
    if n == 1 {
        ans[0] = 1.0;
        return Ok(());
    }

    let mut coeff = vec![0.0_f64; n * n];
    for j in 0..n {
        for i in 0..(n - 1) {
            coeff[j * n + i] = a[j * (n - 1) + i];
        }
        coeff[j * n + n - 1] = 1.0;
    }

    let mut rhs = vec![0.0_f64; n];
    rhs[n - 1] = 1.0;

    let mut free_vars = Vec::new();
    for i in 0..n {
        if vf[i] != 1 {
            free_vars.push(i);
        }
    }
    let nfree = free_vars.len();

    let mut baditer = 0;
    let mut success = false;

    let mut ytrace = 0.0;
    for i in 0..n {
        ytrace += coeff[i * n + i];
    }

    loop {
        if baditer >= 10 {
            break;
        }

        // Compute prod = coeff^T * coeff.
        let mut prod = vec![0.0_f64; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut sum = 0.0;
                for k in 0..n {
                    sum += coeff[i * n + k] * coeff[j * n + k];
                }
                prod[j * n + i] = sum;
            }
        }

        // rr = coeff^T * rhs
        let mut rr = vec![0.0_f64; n];
        for i in 0..n {
            let mut sum = 0.0;
            for k in 0..n {
                sum += coeff[i * n + k] * rhs[k];
            }
            rr[i] = sum;
        }

        if nfree == 0 {
            for i in 0..n {
                ans[i] = 0.0;
            }
            success = true;
            break;
        }

        // Build reduced system
        let mut reduced_prod = vec![0.0_f64; nfree * nfree];
        let mut reduced_rr = vec![0.0_f64; nfree];
        for (i, &free_i) in free_vars.iter().enumerate() {
            reduced_rr[i] = rr[free_i];
            for (j, &free_j) in free_vars.iter().enumerate() {
                reduced_prod[j * nfree + i] = prod[free_j * n + free_i];
            }
        }

        match linsolv(nfree, &mut reduced_prod, &mut reduced_rr) {
            Ok(()) => {
                let sum: f64 = reduced_rr.iter().sum();
                if sum.is_finite() {
                    for i in 0..n {
                        ans[i] = 0.0;
                    }
                    for (i, &free_i) in free_vars.iter().enumerate() {
                        ans[free_i] = reduced_rr[i];
                    }
                    success = true;
                    break;
                }
            }
            Err(_) => {}
        }

        baditer += 1;
        for i in 0..n {
            coeff[i * n + i] += ytrace * 0.001;
        }
    }

    if success {
        let sum: f64 = ans.iter().sum();
        if sum > 0.0 {
            for v in ans.iter_mut() {
                *v /= sum;
            }
        }
        Ok(())
    } else {
        Err(-1)
    }
}
