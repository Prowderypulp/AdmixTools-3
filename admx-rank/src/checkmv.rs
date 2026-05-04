pub fn checkmv(mean: &[f64], var: &[f64], m: usize, n: usize) -> i32 {
    let d = m * n;
    let y1: f64 = mean.iter().map(|&x| x * x).sum();
    let y2: f64 = var.iter().step_by(d + 1).sum();
    if y1.abs() < 1e-12 { return -1; }
    if y2 < 1e-12 { return -2; }
    1
}
