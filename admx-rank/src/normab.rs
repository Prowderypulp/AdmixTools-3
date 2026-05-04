pub fn normab(a: &mut [f64], b: &mut [f64], m: usize, n: usize, rank: usize) {
    for r in 0..rank {
        let bpt = &mut b[r * n .. (r + 1) * n];
        let mut y: f64 = bpt.iter().map(|&x| x * x).sum();
        y = (y / (n as f64)).sqrt();
        let sum_b: f64 = bpt.iter().sum();
        if sum_b < 0.0 {
            y = -y;
        }
        for x in bpt.iter_mut() {
            *x /= y;
        }
        for i in 0..m {
            a[i * rank + r] *= y;
        }
    }
}
