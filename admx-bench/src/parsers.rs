//! Typed parsers for `.fstats` files and qpAdm/qpWave logs.
//!
//! These are intentionally permissive: they extract the load-bearing numerical
//! fields and ignore prose. The C and Rust binaries are both expected to emit
//! the same canonical layout (the changelog's Phase-2 work normalized this).

use std::collections::BTreeMap;
use std::path::Path;

use anyhow::{anyhow, Context, Result};

#[derive(Debug, Clone, Default)]
pub struct FstatsFile {
    pub basepop: Option<String>,
    /// f3 basis means, keyed by sorted (pop_a, pop_b). Stored *1000 (raw file units).
    pub means: BTreeMap<(String, String), f64>,
    /// covariance entries keyed by sorted (a,b),(c,d). Stored *1e6 (raw file units).
    pub covar: BTreeMap<((String, String), (String, String)), f64>,
}

pub fn parse_fstats(path: &Path) -> Result<FstatsFile> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("reading {}", path.display()))?;
    let mut f = FstatsFile::default();
    for line in text.lines() {
        let l = line.trim();
        if l.is_empty() { continue; }
        if l.starts_with("##") {
            if let Some(idx) = l.find("basepop:") {
                let rest = &l[idx + "basepop:".len()..];
                let bp = rest.split("::").next().unwrap_or("").trim();
                if !bp.is_empty() { f.basepop = Some(bp.to_string()); }
            }
            continue;
        }
        let parts: Vec<&str> = l.split_whitespace().collect();
        match parts.len() {
            3 => {
                let v: f64 = parts[2].parse()?;
                f.means.insert(sorted_pair(parts[0], parts[1]), v);
            }
            5 => {
                let v: f64 = parts[4].parse()?;
                let k1 = sorted_pair(parts[0], parts[1]);
                let k2 = sorted_pair(parts[2], parts[3]);
                let key = if k1 <= k2 { (k1, k2) } else { (k2, k1) };
                f.covar.insert(key, v);
            }
            _ => {} // ignore stray lines
        }
    }
    Ok(f)
}

fn sorted_pair(a: &str, b: &str) -> (String, String) {
    if a <= b { (a.to_string(), b.to_string()) } else { (b.to_string(), a.to_string()) }
}

// ---------------------------------------------------------------------------
// qpWave / qpAdm log parsing
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Default)]
pub struct F4RankRow {
    pub rank: i32,
    pub dof: i32,
    pub chisq: f64,
    pub tail: f64,
    pub dofdiff: i32,
    pub chisqdiff: f64,
    pub taildiff: f64,
}

#[derive(Debug, Clone, Default)]
pub struct WaveLog {
    pub rows: Vec<F4RankRow>,
}

pub fn parse_wave_log(path: &Path) -> Result<WaveLog> {
    let text = std::fs::read_to_string(path)?;
    let mut w = WaveLog::default();
    for line in text.lines() {
        let l = line.trim();
        if let Some(row) = parse_f4rank(l)? { w.rows.push(row); }
    }
    Ok(w)
}

fn parse_f4rank(l: &str) -> Result<Option<F4RankRow>> {
    if !l.starts_with("f4rank:") { return Ok(None); }
    // Tokens after 'f4rank:': rank dof: D chisq: C tail: T dofdiff: DD chisqdiff: CD taildiff: TD
    let toks: Vec<&str> = l.split_whitespace().collect();
    let get = |key: &str| -> Result<f64> {
        let i = toks.iter().position(|t| *t == key)
            .ok_or_else(|| anyhow!("missing key {} in: {}", key, l))?;
        toks.get(i + 1).ok_or_else(|| anyhow!("trailing key {}", key))
            .and_then(|s| s.parse::<f64>().map_err(|e| anyhow!("parse {}: {}", s, e)))
    };
    Ok(Some(F4RankRow {
        rank: toks.get(1).ok_or_else(|| anyhow!("no rank"))?.parse()?,
        dof: get("dof:")? as i32,
        chisq: get("chisq:")?,
        tail: get("tail:")?,
        dofdiff: get("dofdiff:")? as i32,
        chisqdiff: get("chisqdiff:")?,
        taildiff: get("taildiff:")?,
    }))
}

#[derive(Debug, Clone, Default)]
pub struct AdmLog {
    pub rows: Vec<F4RankRow>,
    pub coefficients: Vec<f64>,
    pub std_errors: Vec<f64>,
    pub error_cov: Vec<Vec<f64>>,        // K x K
    pub fixed_pat: Vec<FixedPatRow>,
    pub jackknife_mean: Vec<f64>,        // zzjmean
    pub bootstrap_mean: Vec<f64>,        // boot mean
    pub worst_z: Option<f64>,            // worst Z-score with right hand mix
}

#[derive(Debug, Clone, Default)]
pub struct FixedPatRow {
    pub pat: String,
    pub wt: i32,
    pub dof: i32,
    pub chisq: f64,
    pub tail: f64,
    pub coeffs: Vec<f64>,
}

pub fn parse_adm_log(path: &Path) -> Result<AdmLog> {
    let text = std::fs::read_to_string(path)?;
    let mut a = AdmLog::default();
    let mut in_error_cov = false;
    let mut in_fixed_pat = false;
    let mut in_worst_z_block = false;

    for line in text.lines() {
        let l = line.trim();
        if let Some(row) = parse_f4rank(l)? { a.rows.push(row); continue; }

        if let Some(rest) = l.strip_prefix("best coefficients:") {
            a.coefficients = rest.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            continue;
        }
        if let Some(rest) = l.strip_prefix("std. errors:") {
            a.std_errors = rest.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            continue;
        }
        if let Some(rest) = l.strip_prefix("zzjmean") {
            a.jackknife_mean = rest.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            continue;
        }
        if let Some(rest) = l.strip_prefix("boot mean:") {
            a.bootstrap_mean = rest.split_whitespace().filter_map(|s| s.parse().ok()).collect();
            continue;
        }
        if l.starts_with("error covariance") {
            in_error_cov = true;
            a.error_cov.clear();
            continue;
        }
        if in_error_cov {
            if l.is_empty() || !l.chars().next().map_or(false, |c| c == '-' || c.is_ascii_digit()) {
                in_error_cov = false;
            } else {
                let row: Vec<f64> = l.split_whitespace().filter_map(|s| s.parse().ok()).collect();
                if !row.is_empty() { a.error_cov.push(row); }
                continue;
            }
        }
        if l.starts_with("fixed pat") {
            in_fixed_pat = true;
            continue;
        }
        if in_fixed_pat {
            // Format: "<pat> <wt> <dof> <chisq> <tail> <coeffs...>"
            let toks: Vec<&str> = l.split_whitespace().collect();
            if toks.len() >= 5 && toks[0].chars().all(|c| c == '0' || c == '1') {
                let coeffs: Vec<f64> = toks[5..].iter().filter_map(|s| s.parse().ok()).collect();
                a.fixed_pat.push(FixedPatRow {
                    pat: toks[0].to_string(),
                    wt: toks[1].parse().unwrap_or(0),
                    dof: toks[2].parse().unwrap_or(0),
                    chisq: toks[3].parse().unwrap_or(f64::NAN),
                    tail: toks[4].parse().unwrap_or(f64::NAN),
                    coeffs,
                });
            } else if !l.is_empty() && !toks[0].starts_with("best") {
                in_fixed_pat = false;
            }
        }
        if l.starts_with("worst Z-score with right hand mix") {
            in_worst_z_block = true;
            continue;
        }
        if in_worst_z_block && l.contains("Z:") {
            let after = l.split("Z:").nth(1).unwrap_or("");
            if let Some(tok) = after.split_whitespace().next() {
                a.worst_z = tok.parse().ok();
            }
            in_worst_z_block = false;
        }
    }
    Ok(a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn write_tmp(name: &str, body: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        let p = std::env::temp_dir().join(format!("admx-bench-{}-{}.tmp", name, nonce));
        fs::write(&p, body).expect("write temp file");
        p
    }

    #[test]
    fn parse_fstats_canonicalizes_pair_keys() {
        let path = write_tmp(
            "fstats",
            "## basepop: AG ::\nzz BG 1.5\nBG AG QQ zz 2.25\n",
        );
        let parsed = parse_fstats(&path).expect("parse");
        fs::remove_file(&path).ok();

        assert_eq!(parsed.basepop.as_deref(), Some("AG"));
        assert_eq!(parsed.means.get(&(String::from("BG"), String::from("zz"))), Some(&1.5));

        let k = (
            (String::from("AG"), String::from("BG")),
            (String::from("QQ"), String::from("zz")),
        );
        assert_eq!(parsed.covar.get(&k), Some(&2.25));
    }

    #[test]
    fn parse_adm_log_extracts_load_bearing_fields() {
        let path = write_tmp(
            "adm",
            "\
best coefficients: 0.447 0.553
zzjmean 0.114 0.886
boot mean: 0.114 0.886
std. errors: 10.957 10.957
error covariance (* 1,000,000)
120046089 -120046089
-120046089 120046089

fixed pat wt dof chisq tail prob
00 0 2 0.859 0.650754 0.447 0.553
01 1 3 5.909 0.11613 1.000 0.000
10 1 3 6.743 0.08057 0.000 1.000
worst Z-score with right hand mix
f4(Target, Fit, Base, mix of Right pops;  Z: 0.927 sum: 1.000
",
        );
        let parsed = parse_adm_log(&path).expect("parse");
        fs::remove_file(&path).ok();

        assert_eq!(parsed.coefficients, vec![0.447, 0.553]);
        assert_eq!(parsed.jackknife_mean, vec![0.114, 0.886]);
        assert_eq!(parsed.bootstrap_mean, vec![0.114, 0.886]);
        assert_eq!(parsed.std_errors, vec![10.957, 10.957]);
        assert_eq!(parsed.error_cov.len(), 2);
        assert_eq!(parsed.fixed_pat.len(), 3);
        assert_eq!(parsed.fixed_pat[1].pat, "01");
        assert_eq!(parsed.worst_z, Some(0.927));
    }
}
