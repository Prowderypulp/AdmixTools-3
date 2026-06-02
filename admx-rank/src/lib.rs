//! `admx-rank` — Rank test for f4 matrices.
//!
//! Ports `f4rank.c`:
//! - `ranktest` / `ranktestfix` (iterative least-squares, default path)
//! - Closed-form SVD shortcut (gated behind `newmode: YES`)
//! - `normab` normalization
//! - `checkmv` validation

pub mod ranktest;
pub mod normab;
pub mod checkmv;
pub mod stats;
