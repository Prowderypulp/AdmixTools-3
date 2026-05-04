//! `admx-fstats` — Basis-driven f-statistics engine.
//!
//! Provides:
//! - Optimized per-SNP accumulation (`accumulator`)
//! - Canonical basis construction and recovery (`basis`)
//! - `fstats` file reader/writer (`fstats_io`)
//! - High-level `qpfstats` driver (`driver`)

pub mod accumulator;
pub mod basis;
pub mod fstats_io;
pub mod driver;
