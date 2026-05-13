//! Benchmark + parity utilities for admixtools-rs vs C admixtools.
//!
//! Three concerns, three modules:
//!   - `gen`     — deterministic synthetic EIGENSTRAT datasets at multiple tiers.
//!   - `parsers` — typed records for `.fstats` files and qpAdm/qpWave logs.
//!   - `parity`  — tolerance-aware diff between C and Rust outputs.
//!   - `runner`  — invoke a binary on a parfile, repeated, capture wall/RSS.

pub mod gen;
pub mod parsers;
pub mod parity;
pub mod runner;
