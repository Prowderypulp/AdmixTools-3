//! `admx-qpwave` — qpWave driver.
//!
//! Ports `qpWave.c`:
//! - `doq4vecb` with `allsnps` short-circuit semantics
//! - `loadymv` path when `fstatsname:` is set
//! - `checkmv` error handling with exact exit codes
//! - `addscaldiag` call with the exact `yscale` constant

pub mod driver;
