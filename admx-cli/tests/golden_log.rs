//! Golden-log regression tests.
//!
//! These tests run each binary on the fixtures and diff the
//! load-bearing output fields against reference logs from the C binaries.

use std::process::Command;
use std::path::{Path, PathBuf};
use std::fs;

fn get_bin_path(name: &str) -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.pop(); // to admixtools-rs
    path.push("target");
    path.push("debug");
    path.push(name);
    path
}

fn normalize_log(log: &str) -> Vec<String> {
    log.lines()
        .filter_map(|line| {
            let l = line.trim();
            if l.is_empty() 
                || l.starts_with("##") 
                || l.starts_with("###")
                || l.starts_with("parameter file:")
                || l.starts_with("time in")
                || l.contains("seconds cpu")
                || l.contains("Mbytes in use")
                || l.starts_with("before setwt numsnps:")
                || l.starts_with("number of blocks for moving block jackknife:")
                || l.starts_with("setwt numsnps:")
                || l.starts_with("snps:")
                // Ignore lines with 'set' (diagnostics that vary)
                || l.contains(" set ") || l.ends_with(" set")
                // Ignore parameter lines
                || l.contains(":") && !l.starts_with("basis:") && !l.starts_with("pop:") && !l.starts_with("qscore:")
            {
                None
            } else {
                // Normalize multiple spaces to single space
                let norm = l.split_whitespace().collect::<Vec<_>>().join(" ");
                // Normalize floats to 4 decimal places for robust comparison
                let mut parts = Vec::new();
                for part in norm.split_whitespace() {
                    if let Ok(val) = part.parse::<f64>() {
                        parts.push(format!("{:.4}", val));
                    } else {
                        parts.push(part.to_string());
                    }
                }
                Some(parts.join(" "))
            }
        })
        .collect()
}

fn run_golden_test(bin_name: &str, par_file: &Path, expected_log_path: &Path) {
    let bin = get_bin_path(bin_name);
    assert!(bin.exists(), "Binary not found: {:?}", bin);

    let output = Command::new(bin)
        .arg("-p")
        .arg(par_file)
        .current_dir(par_file.parent().unwrap())
        .output()
        .expect("failed to execute process");

    let actual_log = String::from_utf8_lossy(&output.stdout);
    let expected_log = fs::read_to_string(expected_log_path).expect("failed to read expected log");

    let actual_norm = normalize_log(&actual_log);
    let expected_norm = normalize_log(&expected_log);

    if actual_norm != expected_norm {
        // Find first mismatch for better error message
        let mut mismatch_idx = actual_norm.len().min(expected_norm.len());
        for (i, (a, e)) in actual_norm.iter().zip(expected_norm.iter()).enumerate() {
            if a != e {
                mismatch_idx = i;
                break;
            }
        }

        println!("--- ACTUAL (normalized) ---");
        for line in actual_norm.iter().skip(mismatch_idx.saturating_sub(2)).take(5) { println!("{}", line); }
        println!("--- EXPECTED (normalized) ---");
        for line in expected_norm.iter().skip(mismatch_idx.saturating_sub(2)).take(5) { println!("{}", line); }
        
        panic!("Golden log mismatch for {} at line {}", bin_name, mismatch_idx);
    }
}

#[test]
fn test_qpfstats_golden() {
    let mut fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    fixture_path.pop();
    fixture_path.push("tests");
    fixture_path.push("fixtures");
    fixture_path.push("test1");

    let par_file = fixture_path.join("par.qpfstats");
    let expected_log = fixture_path.join("test1.c.no.log");

    run_golden_test("qpfstats", &par_file, &expected_log);
}

#[test]
fn test_qpwave_golden() {
    let mut fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    fixture_path.pop();
    fixture_path.push("tests");
    fixture_path.push("fixtures");
    fixture_path.push("test_wave");

    let par_file = fixture_path.join("par.qpwave");
    let expected_log = fixture_path.join("test123.c.wave.log");

    run_golden_test("qpWave", &par_file, &expected_log);
}

/// qpAdm has RNG-dependent bootstrap lines (std-errors, error covariance), so a
/// full-log diff is not stable. Instead we lock in the DETERMINISTIC outputs that
/// were validated bit-for-bit against the C binary (`/home/drtex/AdmixTools/bin/qpAdm`)
/// on this fixture and on a real 1.23M-SNP dataset: the admixture coefficients, the
/// fixed-pattern chi-square table, and the nested-model p-value.
#[test]
fn test_qpadm_golden() {
    let mut fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    fixture_path.pop();
    fixture_path.push("tests");
    fixture_path.push("fixtures");
    fixture_path.push("test_wave");
    let par_file = fixture_path.join("par.qpadm");

    let bin = get_bin_path("qpAdm");
    assert!(bin.exists(), "Binary not found: {:?}", bin);
    let output = Command::new(&bin)
        .arg("-p")
        .arg(&par_file)
        .current_dir(&fixture_path)
        .output()
        .expect("failed to execute qpAdm");
    let log = String::from_utf8_lossy(&output.stdout);

    // Helper: nth whitespace token of the first line starting with `prefix`, as f64.
    let field = |prefix: &str, idx: usize| -> f64 {
        let line = log
            .lines()
            .find(|l| l.trim_start().starts_with(prefix))
            .unwrap_or_else(|| panic!("no line starting with {:?}\n--- log ---\n{}", prefix, log));
        let tok = line.split_whitespace().nth(idx)
            .unwrap_or_else(|| panic!("line {:?} has no field {}", line, idx));
        tok.parse::<f64>()
            .unwrap_or_else(|_| panic!("field {} of {:?} not a number", idx, line))
    };
    let approx = |got: f64, want: f64, what: &str| {
        assert!((got - want).abs() < 5e-4, "{}: got {}, want {}", what, got, want);
    };

    // Fixed-pattern table rows: "<pat> <wt> <dof> <chisq> <tail> <c0> <c1>".
    // pat 00 = full model.
    approx(field("00 ", 3), 0.859, "pat00 chisq");
    approx(field("00 ", 4), 0.650754, "pat00 tail");
    approx(field("00 ", 5), 0.447, "pat00 coeff0");
    approx(field("00 ", 6), 0.553, "pat00 coeff1");
    // pat 10 / 01 = one source dropped.
    approx(field("10 ", 3), 6.743, "pat10 chisq");
    approx(field("10 ", 4), 0.080570, "pat10 tail");
    approx(field("01 ", 3), 5.909, "pat01 chisq");
    approx(field("01 ", 4), 0.116130, "pat01 tail");

    // Nested-model p-value (best pat line carries it).
    let nested = log
        .lines()
        .find(|l| l.contains("p-value for nested model"))
        .expect("no nested-model line");
    let pval: f64 = nested.split_whitespace().last().unwrap().parse().unwrap();
    assert!((pval - 0.024632).abs() < 5e-4, "nested p-value: got {}", pval);
}
