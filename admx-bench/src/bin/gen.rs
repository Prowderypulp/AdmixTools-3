//! Generate a synthetic EIGENSTRAT dataset tier into a directory.
//!
//! Usage:
//!   admx-gen --tier small --out fixtures/small --seed 42

use std::path::PathBuf;

use admx_bench::gen::{self, DatasetSpec};
use anyhow::{anyhow, Result};
use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about = "Synthetic EIGENSTRAT generator for AdmixTools benchmarks")]
struct Args {
    /// Tier name: tiny | small | med | large
    #[arg(long)]
    tier: String,
    /// Output directory (created if missing).
    #[arg(long)]
    out: PathBuf,
    /// Reproducibility seed.
    #[arg(long, default_value_t = 0xA5A5_A5A5_A5A5_A5A5)]
    seed: u64,
    /// Override stem (default = tier name).
    #[arg(long)]
    stem: Option<String>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let mut spec = DatasetSpec::tier(&args.tier, args.seed)
        .ok_or_else(|| anyhow!("unknown tier '{}': use tiny|small|med|large", args.tier))?;
    if let Some(s) = args.stem { spec.stem = s; }
    let paths = gen::write(&spec, &args.out)?;
    println!("wrote tier={} stem={} dir={}", args.tier, paths.stem, paths.dir.display());
    println!("  n_snps={} n_indivs={}", spec.n_snps, spec.n_indivs());
    Ok(())
}
