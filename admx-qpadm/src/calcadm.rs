fn linsolv_rowmajor(n: usize, matr: &[f64], rhs: &[f64], sol: &mut [f64]) -> bool {
    let mut a = matr.to_vec();
    let mut b = rhs.to_vec();

    for k in 0..(n.saturating_sub(1)) {
        // Pivot search (max abs in column k)
        let mut max_elem = a[k * n + k].abs();
        let mut m = k;
        for i in (k + 1)..n {
            let v = a[i * n + k].abs();
            if max_elem < v {
                // Match C linsolvx pivot-tracking quirk: it compares with abs(...)
                // but stores the raw (possibly negative) matrix entry.
                max_elem = a[i * n + k];
                m = i;
            }
        }

        // Swap rows k and m
        if m != k {
            for i in k..n {
                a.swap(k * n + i, m * n + i);
            }
            b.swap(k, m);
        }

        if a[k * n + k] == 0.0 {
            return false;
        }

        // Elimination
        for j in (k + 1)..n {
            let f_acc = -a[j * n + k] / a[k * n + k];
            for i in k..n {
                a[j * n + i] += f_acc * a[k * n + i];
            }
            b[j] += f_acc * b[k];
        }
    }

    for k in (0..n).rev() {
        let mut f_acc = 0.0;
        for i in (k + 1)..n {
            f_acc += a[k * n + i] * sol[i];
        }
        if a[k * n + k] == 0.0 {
            return false;
        }
        sol[k] = (b[k] - f_acc) / a[k * n + k];
    }

    true
}

/// Weight solver for qpAdm.
/// `ans` gets length `n`. `A` is `n x (n-1)` row-major matrix from `doranktest`.
pub fn calcadm(ans: &mut [f64], a: &[f64], n: usize) -> Result<(), i32> {
    if n == 1 {
        ans[0] = 1.0;
        return Ok(());
    }

    // C path:
    // transpose(coeff, A, n, n-1); vclear(last_row, 1.0, n);
    // with row-major storage and custom Gaussian elimination.
    let mut coeff = vec![0.0_f64; n * n];
    for i in 0..(n - 1) {
        for j in 0..n {
            coeff[i * n + j] = a[j * (n - 1) + i];
        }
    }

    let mut rhs = vec![0.0_f64; n];
    rhs[n - 1] = 1.0;

    let mut baditer = 0;
    let mut success = false;

    // C computes trace before forcing the last row to 1's.
    let mut ytrace = 0.0;
    for i in 0..n {
        ytrace += coeff[i * n + i];
    }
    for j in 0..n {
        coeff[(n - 1) * n + j] = 1.0;
    }

    loop {
        if baditer >= 10 {
            break;
        }

        let mut candidate = vec![0.0; n];
        if linsolv_rowmajor(n, &coeff, &rhs, &mut candidate) {
            let sum: f64 = candidate.iter().sum();
            if sum.is_finite() {
                ans.copy_from_slice(&candidate);
                success = true;
                break;
            }
        }

        baditer += 1;
        for i in 0..(n - 1) {
            coeff[i * n + i] += ytrace * 0.001;
        }
    }

    if !success {
        ans.fill(0.0);
    }

    let sum: f64 = ans.iter().sum();
    if sum > 0.0 {
        for v in ans.iter_mut() {
            *v /= sum;
        }
    }
    Ok(())
}

fn linsolv_colmajor_via_rowmajor(n: usize, a_colmajor: &[f64], rhs: &[f64], sol: &mut [f64]) -> bool {
    // Convert A from column-major to row-major for the local solver.
    let mut a_row = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            a_row[i * n + j] = a_colmajor[j * n + i];
        }
    }
    linsolv_rowmajor(n, &a_row, rhs, sol)
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

        let mut reduced_ans = vec![0.0_f64; nfree];
        if linsolv_colmajor_via_rowmajor(nfree, &reduced_prod, &reduced_rr, &mut reduced_ans) {
            let sum: f64 = reduced_ans.iter().sum();
                if sum.is_finite() {
                    for i in 0..n {
                        ans[i] = 0.0;
                    }
                    for (i, &free_i) in free_vars.iter().enumerate() {
                        ans[free_i] = reduced_ans[i];
                    }
                    success = true;
                    break;
                }
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
