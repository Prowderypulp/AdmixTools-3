//! Numeric constants that must be preserved exactly to maintain output fidelity.
//!
//! These values appear in the C source and surface in printed output and
//! downstream script parsing.  **Do not change them** (see §5 of the spec).

/// Default block size in Morgans for jackknife blocking.
/// C: `blgsize` default in `setblocksz`.
pub const DEFAULT_BLGSIZE: f64 = 0.05;

/// Diagonal regularizer added in `dofstats` to stabilize covariance inversion.
/// C: `diag = 1e-8` in `qpsubs.c`.
pub const DOFSTATS_DIAG: f64 = 1e-8;

/// Scale factor for qpWave variance regularization.
/// C: `yscale = 0.0001` in `qpWave.c`.
pub const QPWAVE_YSCALE: f64 = 0.0001;

/// Default number of bootstrap replicates in qpAdm.
/// C: `numboot = 1000`.
pub const DEFAULT_NUMBOOT: usize = 1000;

/// Maximum number of populations (C: `MAXG = 200`).
pub const MAXG: usize = 200;

/// Maximum number of populations in wave/adm left/right lists (C: `MAXW = 100`).
pub const MAXW: usize = 100;

/// Sentinel value returned by `fstatx` / `getf4` for missing data.
pub const FSTAT_MISSING: f64 = -9999.0;

/// Maximum chromosome number for standard autosomes (C: `MTCHROM = 90`).
pub const MTCHROM: i32 = 90;

/// Pseudo-chromosome number for XY (C: `XYCHROM = 91`).
pub const XYCHROM: i32 = 91;

/// Pseudo-chromosome number for bad/ignored chromosomes (C: `BADCHROM = 99`).
pub const BADCHROM: i32 = 99;
