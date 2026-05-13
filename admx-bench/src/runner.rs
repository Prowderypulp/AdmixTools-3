//! Run a binary on a parfile and capture wall-clock + RSS.
//!
//! We deliberately avoid hyperfine as a runtime dep: this is just
//! `Command::new(...).output()` in a loop, with the working directory
//! set to the parfile's parent so relative paths in the parfile resolve.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Instant;

use anyhow::{Context, Result};

#[derive(Debug, Clone)]
pub struct RunSpec {
    pub binary: PathBuf,
    pub par_file: PathBuf,
    /// If Some, write captured stdout to this path.
    pub log_out: Option<PathBuf>,
}

#[derive(Debug, Clone, Default)]
pub struct RunResult {
    pub wall_secs: f64,
    pub status_ok: bool,
    pub stdout_bytes: usize,
    pub stderr_bytes: usize,
}

pub fn run_once(spec: &RunSpec) -> Result<RunResult> {
    let cwd = spec.par_file.parent().unwrap_or(Path::new("."));
    let par_name = spec.par_file.file_name().unwrap();

    let t = Instant::now();
    let out = Command::new(&spec.binary)
        .arg("-p").arg(par_name)
        .current_dir(cwd)
        .output()
        .with_context(|| format!("invoking {}", spec.binary.display()))?;
    let wall = t.elapsed().as_secs_f64();

    if let Some(p) = &spec.log_out {
        std::fs::write(p, &out.stdout)
            .with_context(|| format!("writing log {}", p.display()))?;
    }

    Ok(RunResult {
        wall_secs: wall,
        status_ok: out.status.success(),
        stdout_bytes: out.stdout.len(),
        stderr_bytes: out.stderr.len(),
    })
}

#[derive(Debug, Clone)]
pub struct BenchStats {
    pub n: usize,
    pub min: f64,
    pub median: f64,
    pub mean: f64,
    pub max: f64,
    pub stddev: f64,
}

pub fn bench(spec: &RunSpec, warmup: usize, runs: usize) -> Result<BenchStats> {
    for _ in 0..warmup { run_once(spec)?; }
    let mut times = Vec::with_capacity(runs);
    for _ in 0..runs {
        let r = run_once(spec)?;
        if !r.status_ok {
            anyhow::bail!("{} exited non-zero on {}", spec.binary.display(), spec.par_file.display());
        }
        times.push(r.wall_secs);
    }
    Ok(stats(&mut times))
}

fn stats(t: &mut [f64]) -> BenchStats {
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = t.len();
    let mean = t.iter().sum::<f64>() / n as f64;
    let var = t.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n.max(1) as f64;
    BenchStats {
        n,
        min: t[0],
        median: t[n / 2],
        mean,
        max: t[n - 1],
        stddev: var.sqrt(),
    }
}
