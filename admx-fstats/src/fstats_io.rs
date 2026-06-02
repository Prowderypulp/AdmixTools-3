//! Ports `dumpfstats`, `dumpfstatshr`, and `loadfstats` from `qpsubs.c`.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use ndarray::{Array1, Array2};
use admx_core::error::{AdmxError, AdmxResult};

/// Saves basis means and covariance to an `.fstats` file.
pub fn dump_fstats<P: AsRef<Path>>(
    path: P,
    means: &Array1<f64>,
    covar: &Array2<f64>,
    pop_labels: &[String],
    basis_indices: &[(usize, usize)], // (a, b) such that stat is f3(anchor, a, b)
    anchor_label: &str,
    hires: bool,
) -> AdmxResult<()> {
    let mut file = File::create(path).map_err(AdmxError::Io)?;
    
    writeln!(file, "##fbasis.  basepop: {} ::  f3*1000 covar*1000000", anchor_label)
        .map_err(AdmxError::Io)?;

    let n = basis_indices.len();

    // Write means
    for i in 0..n {
        let (a, b) = basis_indices[i];
        let val = means[i] * 1000.0;
        if hires {
            writeln!(file, "{:>15} {:>15}  {:12.6}", pop_labels[a], pop_labels[b], val)
                .map_err(AdmxError::Io)?;
        } else {
            writeln!(file, "{:>15} {:>15}  {:9.3}", pop_labels[a], pop_labels[b], val)
                .map_err(AdmxError::Io)?;
        }
    }

    // Write covariance
    for i in 0..n {
        let (a, b) = basis_indices[i];
        for j in i..n {
            let (c, d) = basis_indices[j];
            let val = covar[[i, j]] * 1_000_000.0;
            if hires {
                writeln!(file, "{:>15} {:>15}   {:>15} {:>15}   {:12.6}", 
                    pop_labels[a], pop_labels[b], pop_labels[c], pop_labels[d], val)
                    .map_err(AdmxError::Io)?;
            } else {
                writeln!(file, "{:>15} {:>15}   {:>15} {:>15}   {:9.3}", 
                    pop_labels[a], pop_labels[b], pop_labels[c], pop_labels[d], val)
                    .map_err(AdmxError::Io)?;
            }
        }
    }

    Ok(())
}

/// Loads basis means and covariance from an `.fstats` file.
pub fn load_fstats<P: AsRef<Path>>(
    path: P,
) -> AdmxResult<(Array1<f64>, Array2<f64>, Vec<String>, Vec<(usize, usize)>, String)> {
    let file = File::open(path).map_err(AdmxError::Io)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header = lines.next()
        .ok_or_else(|| AdmxError::Fatal("Empty fstats file".to_string()))?
        .map_err(AdmxError::Io)?;

    if !header.starts_with("##fbasis.") {
        return Err(AdmxError::Fatal("Invalid fstats header".to_string()));
    }

    let anchor_pop = header.split("basepop:")
        .nth(1)
        .and_then(|s| s.split("::").next())
        .map(|s| s.trim().to_string())
        .ok_or_else(|| AdmxError::Fatal("Could not find basepop in fstats header".to_string()))?;

    let mut pop_list = vec![anchor_pop.clone()];
    let mut pop_map = std::collections::HashMap::new();
    pop_map.insert(anchor_pop.clone(), 0);

    let mut get_or_insert_pop = |name: &str| -> usize {
        *pop_map.entry(name.to_string()).or_insert_with(|| {
            let idx = pop_list.len();
            pop_list.push(name.to_string());
            idx
        })
    };

    let mut parsed_means = Vec::new();
    let mut parsed_covars = Vec::new();

    for line in lines {
        let line = line.map_err(AdmxError::Io)?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() == 3 {
            let a = get_or_insert_pop(parts[0]);
            let b = get_or_insert_pop(parts[1]);
            let val: f64 = parts[2].parse().map_err(|_| AdmxError::Fatal("Invalid float".into()))?;
            parsed_means.push(((a, b), val / 1000.0));
        } else if parts.len() == 5 {
            let a = get_or_insert_pop(parts[0]);
            let b = get_or_insert_pop(parts[1]);
            let c = get_or_insert_pop(parts[2]);
            let d = get_or_insert_pop(parts[3]);
            let val: f64 = parts[4].parse().map_err(|_| AdmxError::Fatal("Invalid float".into()))?;
            parsed_covars.push(((a, b, c, d), val / 1_000_000.0));
        }
    }

    let np = pop_list.len();
    let nbasis = (np - 1) * np / 2;

    let mut basis_indices = Vec::with_capacity(nbasis);
    for i in 1..np {
        for j in i..np {
            basis_indices.push((i, j));
        }
    }

    let mut basis_map = std::collections::HashMap::new();
    for (i, &(a, b)) in basis_indices.iter().enumerate() {
        basis_map.insert((a, b), i);
        basis_map.insert((b, a), i);
    }

    let mut means = Array1::zeros(nbasis);
    for ((a, b), val) in parsed_means {
        if let Some(&idx) = basis_map.get(&(a, b)) {
            means[idx] = val;
        }
    }

    let mut covar = Array2::zeros((nbasis, nbasis));
    for ((a, b, c, d), val) in parsed_covars {
        if let (Some(&i), Some(&j)) = (basis_map.get(&(a, b)), basis_map.get(&(c, d))) {
            covar[[i, j]] = val;
            covar[[j, i]] = val;
        }
    }

    Ok((means, covar, pop_list, basis_indices, anchor_pop))
}
