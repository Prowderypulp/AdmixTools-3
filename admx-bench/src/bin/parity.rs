//! Run C vs Rust on the same parfile and emit a tolerance-aware parity report.
//!
//! Usage:
//!   admx-parity --dir fixtures/small --tool qpfstats --c-bin /path/qpfstats --rust-bin target/debug/qpfstats
//!
//! Exits non-zero if any failures exceed tolerance.

use std::path::PathBuf;

use admx_bench::parity::{self, ToleranceProfile};
use admx_bench::parsers;
use admx_bench::runner::{run_once, RunSpec};
use anyhow::{anyhow, Result};
use clap::{Parser, ValueEnum};

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum Tool { Qpfstats, Qpwave, Qpadm }

#[derive(Parser, Debug)]
#[command(version, about = "Parity diff: C vs Rust admixtools outputs")]
struct Args {
    /// Fixture directory (must contain par.qpfstats / par.qpwave / par.qpadm).
    #[arg(long)]
    dir: PathBuf,
    #[arg(long)]
    tool: Tool,
    /// Path to C binary (default: discover under /home/drtex/AdmixTools/bin or PATH).
    #[arg(long)]
    c_bin: Option<PathBuf>,
    /// Path to Rust binary (default: target/debug/<tool>).
    #[arg(long)]
    rust_bin: Option<PathBuf>,
    /// Use stochastic-tolerance profile (qpAdm bootstrap means).
    #[arg(long, default_value_t = false)]
    stochastic: bool,
    /// Where to drop the markdown report (default: <dir>/parity_<tool>.md).
    #[arg(long)]
    report: Option<PathBuf>,
}

fn main() -> Result<()> {
    let args = Args::parse();
    let (par_name, log_suffix) = match args.tool {
        Tool::Qpfstats => ("par.qpfstats", "qpfstats"),
        Tool::Qpwave   => ("par.qpwave",   "qpWave"),
        Tool::Qpadm    => ("par.qpadm",    "qpAdm"),
    };
    let par = args.dir.join(par_name);
    if !par.exists() { return Err(anyhow!("missing {}", par.display())); }

    let c_bin = args.c_bin.unwrap_or_else(|| default_c_bin(log_suffix));
    let rust_bin = args.rust_bin.unwrap_or_else(|| default_rust_bin(log_suffix));
    if !c_bin.exists()    { return Err(anyhow!("C binary not found: {}",    c_bin.display())); }
    if !rust_bin.exists() { return Err(anyhow!("Rust binary not found: {}", rust_bin.display())); }

    let c_log = args.dir.join(format!("{}.c.log", log_suffix));
    let r_log = args.dir.join(format!("{}.rust.log", log_suffix));

    let tol = if args.stochastic { ToleranceProfile::stochastic() } else { ToleranceProfile::strict() };

    let report = match args.tool {
        Tool::Qpfstats => {
            // qpfstats writes a .fstats sibling. Capture C and Rust snapshots
            // after each run so we don't rerun tools or rely on rename ordering.
            let stem = read_param(&par, "fstatsoutname").unwrap_or_else(|| "out.fstats".into());
            let fst = args.dir.join(&stem);
            let c_fst = args.dir.join(format!("{}.c", stem));
            let r_fst = args.dir.join(format!("{}.rust", stem));

            println!("running C   : {}", c_bin.display());
            run_once(&RunSpec { binary: c_bin, par_file: par.clone(), log_out: Some(c_log.clone()) })?;
            std::fs::copy(&fst, &c_fst)
                .map_err(|e| anyhow!("copying C fstats {} -> {}: {}", fst.display(), c_fst.display(), e))?;

            println!("running rust: {}", rust_bin.display());
            run_once(&RunSpec { binary: rust_bin, par_file: par.clone(), log_out: Some(r_log.clone()) })?;
            std::fs::copy(&fst, &r_fst)
                .map_err(|e| anyhow!("copying Rust fstats {} -> {}: {}", fst.display(), r_fst.display(), e))?;

            let c = parsers::parse_fstats(&c_fst)?;
            let r = parsers::parse_fstats(&r_fst)?;
            parity::diff_fstats(&c, &r, &tol)
        }
        Tool::Qpwave => {
            println!("running C   : {}", c_bin.display());
            run_once(&RunSpec { binary: c_bin, par_file: par.clone(), log_out: Some(c_log.clone()) })?;
            println!("running rust: {}", rust_bin.display());
            run_once(&RunSpec { binary: rust_bin, par_file: par.clone(), log_out: Some(r_log.clone()) })?;
            let c = parsers::parse_wave_log(&c_log)?;
            let r = parsers::parse_wave_log(&r_log)?;
            parity::diff_wave(&c, &r, &tol)
        }
        Tool::Qpadm => {
            println!("running C   : {}", c_bin.display());
            run_once(&RunSpec { binary: c_bin, par_file: par.clone(), log_out: Some(c_log.clone()) })?;
            println!("running rust: {}", rust_bin.display());
            run_once(&RunSpec { binary: rust_bin, par_file: par.clone(), log_out: Some(r_log.clone()) })?;
            let c = parsers::parse_adm_log(&c_log)?;
            let r = parsers::parse_adm_log(&r_log)?;
            parity::diff_adm(&c, &r, &tol)
        }
    };

    let report_path = args.report.unwrap_or_else(|| args.dir.join(format!("parity_{}.md", log_suffix)));
    std::fs::write(&report_path, report.render_markdown())?;
    println!("report: {}", report_path.display());
    for w in &report.worst {
        println!("  worst {:<24} {:.3e}", w.field, w.residual);
    }
    if report.ok() { println!("PARITY: PASS"); }
    else {
        println!("PARITY: FAIL ({} failures)", report.failures.len());
        std::process::exit(1);
    }
    Ok(())
}

fn default_c_bin(name: &str) -> PathBuf {
    let standard = PathBuf::from("/home/drtex/AdmixTools/bin").join(name);
    if standard.exists() { return standard; }
    PathBuf::from(name)
}

fn default_rust_bin(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.pop();
    p.push("target");
    p.push("debug");
    p.push(name);
    p
}

fn read_param(par: &std::path::Path, key: &str) -> Option<String> {
    let text = std::fs::read_to_string(par).ok()?;
    for line in text.lines() {
        let l = line.trim();
        if let Some(rest) = l.strip_prefix(&format!("{}:", key)) {
            return Some(rest.trim().to_string());
        }
    }
    None
}
