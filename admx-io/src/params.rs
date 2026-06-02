//! Parameter file parser.
//!
//! Ports the `getpars` logic from `nicksrc/getpars.c`.
//! Only the keys actually read by qpfstats, qpWave, and qpAdm are supported.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use admx_core::error::{AdmxError, AdmxResult};

/// Parameter file handle.
#[derive(Debug, Default, Clone)]
pub struct ParFile {
    params: HashMap<String, String>,
}

impl ParFile {
    /// Load parameters from a file. Supports `+++` includes and string substitution.
    pub fn load<P: AsRef<Path>>(path: P) -> AdmxResult<Self> {
        let mut par_file = Self::default();
        par_file.read_recursive(path.as_ref(), 0)?;
        par_file.do_str_sub();
        Ok(par_file)
    }

    fn read_recursive(&mut self, path: &Path, depth: usize) -> AdmxResult<()> {
        if depth >= 100 {
            return Err(AdmxError::Fatal("include files looping".to_string()));
        }

        let file = File::open(path).map_err(AdmxError::Io)?;
        let reader = BufReader::new(file);

        for line in reader.lines() {
            let line = line.map_err(AdmxError::Io)?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            // Strip inline comments
            let line = if let Some(pos) = line.find('#') {
                &line[..pos]
            } else {
                line
            };

            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }

            let key = parts[0];
            let rest = &line[key.len()..].trim();

            if key == "+++" {
                let include_path = PathBuf::from(rest);
                // If it's relative, it should be relative to the current parfile?
                // AdmixTools C code just does fopen(rest, "r").
                // Usually this means relative to CWD.
                self.read_recursive(&include_path, depth + 1)?;
                continue;
            }

            if !key.ends_with(':') {
                eprintln!("**warning: dubious parameter, please check the parameter name {}", key);
            }

            if self.params.contains_key(key) {
                log::warn!("duplicate parameter: {} overriding", key);
            }

            self.params.insert(key.to_string(), rest.to_string());
        }

        Ok(())
    }

    /// String substitution for UPPER case parameters.
    fn do_str_sub(&mut self) {
        // Ports dostrsub from getpars.c
        // It's a recursive process in C, but we can do it iteratively or until convergence.
        loop {
            let mut changed = false;
            
            // Identify uppercase keys (excluding the trailing colon if present)
            let mut upper_params = Vec::new();
            for (key, val) in &self.params {
                let stem = key.strip_suffix(':').unwrap_or(key);
                if !stem.is_empty() && stem.chars().all(|c| !c.is_lowercase() && (c.is_alphanumeric() || c == '_')) {
                    upper_params.push((stem.to_string(), val.clone()));
                }
            }

            // Sort longest first to match C's stable sort behavior (sort of)
            upper_params.sort_by(|a, b| b.0.len().cmp(&a.0.len()));

            let mut new_params = self.params.clone();
            for (upper_key, upper_val) in upper_params {
                for (_key, val) in new_params.iter_mut() {
                    // C code only substitutes in "lowercase" (non-fully-uppercase) parameters?
                    // Actually, it substitutes in ANY parameter's DATA.
                    if val.contains(&upper_key) {
                        let replaced = val.replace(&upper_key, &upper_val);
                        if replaced != *val {
                            *val = replaced;
                            changed = true;
                        }
                    }
                }
            }

            self.params = new_params;

            if !changed {
                break;
            }
        }
    }

    pub fn get_string(&self, key: &str) -> Option<String> {
        self.params.get(key).map(|s| {
            // AdmixTools getstring returns only the FIRST word of the value
            s.split_whitespace().next().unwrap_or("").to_string()
        })
    }

    pub fn get_long_string(&self, key: &str) -> Option<String> {
        self.params.get(key).cloned()
    }

    pub fn get_int(&self, key: &str) -> Option<i32> {
        self.get_string(key).and_then(|s| {
            if s == "YES" {
                Some(1)
            } else if s == "NO" {
                Some(0)
            } else {
                s.parse().ok()
            }
        })
    }

    pub fn get_dbl(&self, key: &str) -> Option<f64> {
        self.get_string(key).and_then(|s| s.parse().ok())
    }

    pub fn get_bool(&self, key: &str, default: bool) -> bool {
        self.get_int(key).map(|v| v != 0).unwrap_or(default)
    }

    pub fn get_list(&self, key: &str) -> Option<Vec<String>> {
        self.get_string(key).map(|path| {
            if let Ok(content) = std::fs::read_to_string(&path) {
                content.lines()
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty() && !s.starts_with('#'))
                    .collect()
            } else {
                // If it's not a file, it could be a comma-separated list or space-separated list
                path.split(&[',', ' '][..])
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect()
            }
        })
    }

    /// Write parameters to stdout (useful for echoing the parfile).
    pub fn write_pars(&self) {
        let mut keys: Vec<_> = self.params.keys().collect();
        keys.sort();
        for key in keys {
            println!("{} {}", key, self.params[key]);
        }
    }
}
