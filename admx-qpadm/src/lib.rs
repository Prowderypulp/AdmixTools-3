//! `admx-qpadm` — qpAdm driver.
//!
//! Ports `qpAdm.c`:
//! - `calcadm` / `calcadmfix` — admixture weight solver
//! - Fix-pattern enumeration (`2^nl` masks, parallelized)
//! - `calcevar` (delete-block) and `calcevarboot` (bootstrap)
//! - Legacy LCG PRNG for reproducible bootstrap
//! - Nested-model p-values, `summ:` line, hires output
//! - Direct call into `admx-fstats` instead of `system("qpfstats ...")`

pub mod driver;
pub mod calcadm;
pub mod bootstrap;
pub mod prng;
