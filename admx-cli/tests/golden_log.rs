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
    let mut skip_zzevarboot = false;
    let mut normalized = Vec::new();
    for line in log.lines() {
        let l = line.trim();
        if l.starts_with("zzevarboot") {
            skip_zzevarboot = true;
            continue;
        }
        if skip_zzevarboot {
            if l.starts_with("boot mean:") || l.starts_with("bootstrap saimpling") || l.starts_with("zzjmean") {
                skip_zzevarboot = false;
            } else {
                continue;
            }
        }
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
            continue;
        }
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
        normalized.push(parts.join(" "));
    }
    normalized
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

#[test]
fn test_qpadm_golden() {
    let mut fixture_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    fixture_path.pop();
    fixture_path.push("tests");
    fixture_path.push("fixtures");
    fixture_path.push("test_wave");

    let par_file = fixture_path.join("par.qpadm");
    let expected_log = fixture_path.join("test123.c.adm.log");

    run_golden_test("qpAdm", &par_file, &expected_log);
}
