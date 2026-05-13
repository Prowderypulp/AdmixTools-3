//! Tolerance-aware diff between C and Rust outputs.
//!
//! Tolerances are intentionally tight. They reflect what we actually expect
//! once Phase-2 parity is preserved: f-stat means agree to a few ULPs in
//! file-units (×1000), covariances to ~1e-6 relative, qpAdm coefficients to
//! 1e-4, and chi-sq tails to 1e-3 absolute.

use crate::parsers::{AdmLog, FstatsFile, WaveLog};

#[derive(Debug, Clone, Default)]
pub struct ToleranceProfile {
    pub fstats_mean_abs: f64,
    pub fstats_cov_rel: f64,
    pub fstats_cov_abs_floor: f64,
    pub chisq_abs: f64,
    pub tail_abs: f64,
    pub coeff_abs: f64,
    pub stderr_rel: f64,
    pub jack_abs: f64,
    pub worst_z_abs: f64,
}

impl ToleranceProfile {
    pub fn strict() -> Self {
        Self {
            fstats_mean_abs: 1e-6,        // file units (×1000)
            fstats_cov_rel: 1e-6,
            fstats_cov_abs_floor: 1e-6,
            chisq_abs: 1e-4,
            tail_abs: 1e-4,
            coeff_abs: 1e-4,
            stderr_rel: 1e-3,
            jack_abs: 1e-3,
            worst_z_abs: 1e-3,
        }
    }
    /// Looser profile for stochastic outputs (bootstrap means).
    pub fn stochastic() -> Self {
        Self { coeff_abs: 5e-3, jack_abs: 5e-3, worst_z_abs: 5e-3, ..Self::strict() }
    }
}

#[derive(Debug, Default)]
pub struct ParityReport {
    pub failures: Vec<Failure>,
    pub worst: Vec<Worst>,
}

#[derive(Debug, Clone)]
pub struct Failure {
    pub field: String,
    pub key: String,
    pub c_value: f64,
    pub rust_value: f64,
    pub residual: f64,
    pub tolerance: f64,
}

#[derive(Debug, Clone)]
pub struct Worst { pub field: String, pub residual: f64 }

impl ParityReport {
    pub fn ok(&self) -> bool { self.failures.is_empty() }

    pub fn render_markdown(&self) -> String {
        let mut s = String::new();
        s.push_str("# Parity report\n\n");
        if self.ok() { s.push_str("**PASS** — all fields within tolerance.\n\n"); }
        else { s.push_str(&format!("**FAIL** — {} field(s) exceed tolerance.\n\n", self.failures.len())); }
        if !self.worst.is_empty() {
            s.push_str("## Worst residuals (informational)\n\n| field | residual |\n|---|---|\n");
            for w in &self.worst {
                s.push_str(&format!("| {} | {:.3e} |\n", w.field, w.residual));
            }
            s.push('\n');
        }
        if !self.failures.is_empty() {
            s.push_str("## Failures\n\n| field | key | c | rust | residual | tol |\n|---|---|---|---|---|---|\n");
            for f in &self.failures {
                s.push_str(&format!("| {} | {} | {:.6e} | {:.6e} | {:.3e} | {:.1e} |\n",
                    f.field, f.key, f.c_value, f.rust_value, f.residual, f.tolerance));
            }
        }
        s
    }
}

pub fn diff_fstats(c: &FstatsFile, r: &FstatsFile, tol: &ToleranceProfile) -> ParityReport {
    let mut rep = ParityReport::default();
    let mut worst_mean = 0.0_f64;
    let mut worst_cov = 0.0_f64;

    for (k, &cv) in &c.means {
        match r.means.get(k) {
            Some(&rv) => {
                let d = (cv - rv).abs();
                if d > worst_mean { worst_mean = d; }
                if d > tol.fstats_mean_abs {
                    rep.failures.push(Failure {
                        field: "fstats.mean".into(),
                        key: format!("{}|{}", k.0, k.1),
                        c_value: cv, rust_value: rv, residual: d,
                        tolerance: tol.fstats_mean_abs,
                    });
                }
            }
            None => rep.failures.push(Failure {
                field: "fstats.mean".into(), key: format!("{}|{}", k.0, k.1),
                c_value: cv, rust_value: f64::NAN, residual: f64::INFINITY,
                tolerance: tol.fstats_mean_abs,
            }),
        }
    }
    for (k, &cv) in &c.covar {
        if let Some(&rv) = r.covar.get(k) {
            let scale = cv.abs().max(rv.abs()).max(tol.fstats_cov_abs_floor);
            let rel = (cv - rv).abs() / scale;
            if rel > worst_cov { worst_cov = rel; }
            if rel > tol.fstats_cov_rel {
                rep.failures.push(Failure {
                    field: "fstats.cov".into(),
                    key: format!("({},{})|({},{})", (k.0).0, (k.0).1, (k.1).0, (k.1).1),
                    c_value: cv, rust_value: rv, residual: rel,
                    tolerance: tol.fstats_cov_rel,
                });
            }
        }
    }
    rep.worst.push(Worst { field: "fstats.mean".into(), residual: worst_mean });
    rep.worst.push(Worst { field: "fstats.cov".into(), residual: worst_cov });
    rep
}

pub fn diff_wave(c: &WaveLog, r: &WaveLog, tol: &ToleranceProfile) -> ParityReport {
    let mut rep = ParityReport::default();
    let n = c.rows.len().min(r.rows.len());
    if c.rows.len() != r.rows.len() {
        rep.failures.push(Failure {
            field: "wave.rows".into(), key: "len".into(),
            c_value: c.rows.len() as f64, rust_value: r.rows.len() as f64,
            residual: f64::INFINITY, tolerance: 0.0,
        });
    }
    let mut worst_chisq = 0.0_f64;
    let mut worst_tail = 0.0_f64;
    for i in 0..n {
        let (cr, rr) = (&c.rows[i], &r.rows[i]);
        let dc = (cr.chisq - rr.chisq).abs();
        let dt = (cr.tail - rr.tail).abs();
        if dc > worst_chisq { worst_chisq = dc; }
        if dt > worst_tail { worst_tail = dt; }
        if dc > tol.chisq_abs {
            rep.failures.push(Failure {
                field: "wave.chisq".into(), key: format!("rank={}", cr.rank),
                c_value: cr.chisq, rust_value: rr.chisq, residual: dc,
                tolerance: tol.chisq_abs,
            });
        }
        if dt > tol.tail_abs {
            rep.failures.push(Failure {
                field: "wave.tail".into(), key: format!("rank={}", cr.rank),
                c_value: cr.tail, rust_value: rr.tail, residual: dt,
                tolerance: tol.tail_abs,
            });
        }
    }
    rep.worst.push(Worst { field: "wave.chisq".into(), residual: worst_chisq });
    rep.worst.push(Worst { field: "wave.tail".into(), residual: worst_tail });
    rep
}

pub fn diff_adm(c: &AdmLog, r: &AdmLog, tol: &ToleranceProfile) -> ParityReport {
    let mut rep = ParityReport::default();

    fn cmp_vec(rep: &mut ParityReport, field: &str, c: &[f64], r: &[f64], tol_abs: f64) -> f64 {
        let mut worst = 0.0_f64;
        let n = c.len().min(r.len());
        if c.len() != r.len() {
            rep.failures.push(Failure {
                field: format!("{}.len", field), key: "".into(),
                c_value: c.len() as f64, rust_value: r.len() as f64,
                residual: f64::INFINITY, tolerance: 0.0,
            });
        }
        for i in 0..n {
            let d = (c[i] - r[i]).abs();
            if d > worst { worst = d; }
            if d > tol_abs {
                rep.failures.push(Failure {
                    field: field.to_string(), key: format!("[{}]", i),
                    c_value: c[i], rust_value: r[i], residual: d, tolerance: tol_abs,
                });
            }
        }
        worst
    }

    let w_coef = cmp_vec(&mut rep, "adm.coeff", &c.coefficients, &r.coefficients, tol.coeff_abs);
    let w_jack = cmp_vec(&mut rep, "adm.zzjmean", &c.jackknife_mean, &r.jackknife_mean, tol.jack_abs);
    let w_boot = cmp_vec(&mut rep, "adm.bootmean", &c.bootstrap_mean, &r.bootstrap_mean, tol.jack_abs);

    // std. errors: relative tol
    let mut worst_se = 0.0_f64;
    let n_se = c.std_errors.len().min(r.std_errors.len());
    for i in 0..n_se {
        let denom = c.std_errors[i].abs().max(1e-9);
        let rel = (c.std_errors[i] - r.std_errors[i]).abs() / denom;
        if rel > worst_se { worst_se = rel; }
        if rel > tol.stderr_rel {
            rep.failures.push(Failure {
                field: "adm.stderr".into(), key: format!("[{}]", i),
                c_value: c.std_errors[i], rust_value: r.std_errors[i],
                residual: rel, tolerance: tol.stderr_rel,
            });
        }
    }

    // fixed pat: match by pattern string
    let mut worst_pat_chisq = 0.0_f64;
    let mut worst_pat_tail = 0.0_f64;
    for cr in &c.fixed_pat {
        if let Some(rr) = r.fixed_pat.iter().find(|p| p.pat == cr.pat) {
            let dc = (cr.chisq - rr.chisq).abs();
            let dt = (cr.tail - rr.tail).abs();
            if dc > worst_pat_chisq { worst_pat_chisq = dc; }
            if dt > worst_pat_tail { worst_pat_tail = dt; }
            if dc > tol.chisq_abs {
                rep.failures.push(Failure {
                    field: "adm.fixed_pat.chisq".into(), key: cr.pat.clone(),
                    c_value: cr.chisq, rust_value: rr.chisq, residual: dc,
                    tolerance: tol.chisq_abs,
                });
            }
            if dt > tol.tail_abs {
                rep.failures.push(Failure {
                    field: "adm.fixed_pat.tail".into(), key: cr.pat.clone(),
                    c_value: cr.tail, rust_value: rr.tail, residual: dt,
                    tolerance: tol.tail_abs,
                });
            }
        } else {
            rep.failures.push(Failure {
                field: "adm.fixed_pat".into(), key: cr.pat.clone(),
                c_value: f64::NAN, rust_value: f64::NAN, residual: f64::INFINITY,
                tolerance: 0.0,
            });
        }
    }

    if let (Some(cz), Some(rz)) = (c.worst_z, r.worst_z) {
        let d = (cz - rz).abs();
        if d > tol.worst_z_abs {
            rep.failures.push(Failure {
                field: "adm.worst_z".into(), key: "".into(),
                c_value: cz, rust_value: rz, residual: d, tolerance: tol.worst_z_abs,
            });
        }
        rep.worst.push(Worst { field: "adm.worst_z".into(), residual: d });
    }

    rep.worst.push(Worst { field: "adm.coeff".into(), residual: w_coef });
    rep.worst.push(Worst { field: "adm.stderr_rel".into(), residual: worst_se });
    rep.worst.push(Worst { field: "adm.zzjmean".into(), residual: w_jack });
    rep.worst.push(Worst { field: "adm.bootmean".into(), residual: w_boot });
    rep.worst.push(Worst { field: "adm.fixed_pat.chisq".into(), residual: worst_pat_chisq });
    rep.worst.push(Worst { field: "adm.fixed_pat.tail".into(), residual: worst_pat_tail });
    rep
}

pub fn merge(reports: Vec<ParityReport>) -> ParityReport {
    let mut out = ParityReport::default();
    for r in reports { out.failures.extend(r.failures); out.worst.extend(r.worst); }
    out
}
