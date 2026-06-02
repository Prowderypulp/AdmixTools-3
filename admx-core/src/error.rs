//! Core error types for the AdmixTools pipeline.
//!
//! Legacy C code uses `fatalx()` with format strings and `exit(1)`.
//! We mirror those with typed errors that preserve the original message
//! strings and exit codes so downstream scripts see the same behavior.

use thiserror::Error;

/// Top-level error type for AdmixTools operations.
#[derive(Debug, Error)]
pub enum AdmxError {
    #[error("FATAL: {0}")]
    Fatal(String),

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),

    #[error("parse error in parameter file: {0}")]
    ParamParse(String),

    #[error("genotype file format error: {0}")]
    GenoFormat(String),

    #[error("SNP/indiv file error: {0}")]
    SnpIndiv(String),

    #[error("linear algebra error: {0}")]
    Linalg(String),

    #[error("rank test did not converge after {0} iterations")]
    RankConvergence(usize),

    #[error("checkmv: bad mean/variance matrix — {0}")]
    BadMeanVar(String),
}

/// Result alias used throughout the crate.
pub type AdmxResult<T> = Result<T, AdmxError>;
