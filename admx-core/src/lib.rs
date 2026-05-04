//! `admx-core` — Core types, constants, jackknife, and LAPACK wrappers.
//!
//! This crate defines the fundamental data structures used across the entire
//! AdmixTools Rust port. All types mirror their C counterparts in `admutils.h`
//! and `qpsubs.h` for fidelity.

// Re-export openblas-src so downstream crates get the linker symbols.
extern crate openblas_src;

pub mod types;
pub mod constants;
pub mod jackknife;
pub mod linalg;
pub mod error;
pub mod blocks;
