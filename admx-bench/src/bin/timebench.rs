//! Wallclock benchmark: pair C vs Rust on identical parfiles, repeated.
//!
//! Usage:
//!   admx-timebench --dir fixtures/small --tool qpfstats --warmup 1 --runs 5
//!
//! Emits a markdown table + JSON to stdout / --json-out.

use std::path::PathBuf;

use admx_bench::runner::{bench, BenchStats, RunSpec};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueEnum};
use serde::Serialize;

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Tool { Qpfstats, Qpwave, Qpadm, All }

#[derive(Parser, Debug)]
#[command(version, about = "Wallclock benchmark: C vs Rust admixtools")]
struct Args {
    #[arg(long)] dir: PathBuf,
    #[arg(long, default_value = "all")] tool: Tool,
    #[arg(long, default_value_t = 1)] warmup: usize,
    #[arg(long, default_value_t = 5)] runs: usize,
    #[arg(long)] c_bin_dir: Option<PathBuf>,
    #[arg(long)] rust_bin_dir: Option<PathBuf>,
    #[arg(long)] json_out: Option<PathBuf>,
}

#[derive(Debug, Serialize)]
struct ToolResult {
    tool: String,
    c: BenchStatsSer,
    rust: BenchStatsSer,
    speedup_median: f64,
    speedup_mean: f64,
}

#[derive(Debug, Serialize)]
struct BenchStatsSer { n: usize, min: f64, median: f64, mean: f64, max: f64, stddev: f64 }
impl From<BenchStats> for BenchStatsSer {
    fn from(s: BenchStats) -> Self {
        BenchStatsSer { n: s.n, min: s.min, median: s.median, mean: s.mean, max: s.max, stddev: s.stddev }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let c_dir = args.c_bin_dir.unwrap_or_else(|| PathBuf::from("/home/drtex/AdmixTools/bin"));
    let rust_dir = args.rust_bin_dir.unwrap_or_else(|| {
        let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        p.pop(); p.push("target"); p.push("release"); p
    });

    let tools: Vec<(&str, &str)> = match args.tool {
        Tool::All      => vec![("qpfstats", "par.qpfstats"), ("qpWave", "par.qpwave"), ("qpAdm", "par.qpadm")],
        Tool::Qpfstats => vec![("qpfstats", "par.qpfstats")],
        Tool::Qpwave   => vec![("qpWave",   "par.qpwave")],
        Tool::Qpadm    => vec![("qpAdm",    "par.qpadm")],
    };

    let mut results: Vec<ToolResult> = Vec::new();
    println!("| tool | C median (s) | Rust median (s) | speedup (median) | C±sd | Rust±sd |");
    println!("|---|---|---|---|---|---|");
    for (bin, par) in tools {
        let c_bin = c_dir.join(bin);
        let r_bin = rust_dir.join(bin);
        let par_path = args.dir.join(par);
        if !c_bin.exists()   { return Err(anyhow!("missing C binary {}",    c_bin.display())); }
        if !r_bin.exists()   { return Err(anyhow!("missing Rust binary {}", r_bin.display())); }
        if !par_path.exists(){ return Err(anyhow!("missing parfile {}", par_path.display())); }

        let c = bench(&RunSpec { binary: c_bin, par_file: par_path.clone(), log_out: None },
                      args.warmup, args.runs)?;
        let r = bench(&RunSpec { binary: r_bin, par_file: par_path.clone(), log_out: None },
                      args.warmup, args.runs)?;

        let sp_med = c.median / r.median;
        let sp_mean = c.mean / r.mean;
        println!("| {} | {:.3} | {:.3} | {:.2}× | {:.3}±{:.3} | {:.3}±{:.3} |",
                 bin, c.median, r.median, sp_med, c.mean, c.stddev, r.mean, r.stddev);
        results.push(ToolResult {
            tool: bin.to_string(),
            c: c.into(), rust: r.into(),
            speedup_median: sp_med, speedup_mean: sp_mean,
        });
    }
    if let Some(p) = args.json_out {
        std::fs::write(&p, serde_json::to_string_pretty(&results)?)?;
        eprintln!("json: {}", p.display());
    }
    Ok(())
}
