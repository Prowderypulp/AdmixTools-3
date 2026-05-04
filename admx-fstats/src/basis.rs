//! F-statistics basis construction and recovery.
//!
//! Handles the transformation from per-population sum-of-products $M_{ab}$
//! to canonical f-statistics ($f2, f3, f4$) and basis propagation.

use ndarray::Array2;
use admx_core::types::FStatKind;

/// Computes an f-statistic from an unbiased sum-of-products matrix $M$.
/// 
/// $M_{ab}$ is expected to be the mean across SNPs of $(p_a p_b - \delta_{ab} aax_a)$.
pub fn compute_fstat(kind: FStatKind, m: &Array2<f64>) -> f64 {
    match kind {
        FStatKind::F2(a, b) => {
            m[[a, a]] + m[[b, b]] - 2.0 * m[[a, b]]
        }
        FStatKind::F3(anchor, a, b) => {
            m[[anchor, anchor]] - m[[anchor, a]] - m[[anchor, b]] + m[[a, b]]
        }
        FStatKind::F4(a, b, c, d) => {
            m[[a, c]] - m[[a, d]] - m[[b, c]] + m[[b, d]]
        }
    }
}

/// Canonical basis for `qpfstats`.
/// Usually all $f3(base; a, b)$ for $1 \le a \le b < np$.
pub struct FBasis {
    pub anchor: usize,
    pub pops: Vec<usize>,
}

impl FBasis {
    pub fn new(anchor: usize, num_pops: usize) -> Self {
        let mut pops = Vec::new();
        for i in 0..num_pops {
            if i != anchor {
                pops.push(i);
            }
        }
        Self { anchor, pops }
    }

    pub fn stats(&self) -> Vec<FStatKind> {
        let mut stats = Vec::new();
        for i in 0..self.pops.len() {
            let a = self.pops[i];
            for j in i..self.pops.len() {
                let b = self.pops[j];
                stats.push(FStatKind::F3(self.anchor, a, b));
            }
        }
        stats
    }

    pub fn coefficients(&self, kind: FStatKind) -> Vec<(usize, f64)> {
        // Map fstat to indices in self.stats()
        // stats() returns F3(anchor, a, b) for 1 <= a <= b < np
        let mut co = Vec::new();
        let n = self.pops.len();

        let get_idx = |p1: usize, p2: usize| {
            if p1 == self.anchor || p2 == self.anchor { return None; }
            let (a, b) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
            // Find (a, b) in self.stats()
            // pops are 1..np relative to anchor 0? No, they are absolute indices.
            let i = self.pops.iter().position(|&x| x == a)?;
            let j = self.pops.iter().position(|&x| x == b)?;
            // index in stats() is sum(k=0..i) (n-k) + (j-i)
            let mut idx = 0;
            for k in 0..i { idx += n - k; }
            idx += j - i;
            Some(idx)
        };

        match kind {
            FStatKind::F2(a, b) => {
                // f2(a, b) = f3(0, a, a) + f3(0, b, b) - 2 f3(0, a, b)
                if let Some(idx) = get_idx(a, a) { co.push((idx, 1.0)); }
                if let Some(idx) = get_idx(b, b) { co.push((idx, 1.0)); }
                if let Some(idx) = get_idx(a, b) { co.push((idx, -2.0)); }
            }
            FStatKind::F3(c, a, b) => {
                // f3(c, a, b) = f3(0, a, b) - f3(0, c, a) - f3(0, c, b) + f3(0, c, c)
                if let Some(idx) = get_idx(a, b) { co.push((idx, 1.0)); }
                if let Some(idx) = get_idx(c, a) { co.push((idx, -1.0)); }
                if let Some(idx) = get_idx(c, b) { co.push((idx, -1.0)); }
                if let Some(idx) = get_idx(c, c) { co.push((idx, 1.0)); }
            }
            FStatKind::F4(a, b, c, d) => {
                // f4(a, b, c, d) = f3(0, a, c) + f3(0, b, d) - f3(0, a, d) - f3(0, b, c)
                if let Some(idx) = get_idx(a, c) { co.push((idx, 1.0)); }
                if let Some(idx) = get_idx(b, d) { co.push((idx, 1.0)); }
                if let Some(idx) = get_idx(a, d) { co.push((idx, -1.0)); }
                if let Some(idx) = get_idx(b, c) { co.push((idx, -1.0)); }
            }
        }
        co
    }

    pub fn full_stats(num_pops: usize) -> Vec<FStatKind> {
        let mut stats = Vec::new();
        // all f2
        for a in 0..num_pops {
            for b in a+1..num_pops {
                stats.push(FStatKind::F2(a, b));
            }
        }
        // f3 anchored at each pop
        for c in 0..num_pops {
            for a in 0..num_pops {
                if a == c { continue; }
                for b in a+1..num_pops {
                    if b == c { continue; }
                    stats.push(FStatKind::F3(c, a, b));
                }
            }
        }
        // all canonical f4
        for a in 0..num_pops {
            for b in a+1..num_pops {
                for c in a+1..num_pops {
                    if c == b { continue; }
                    for d in c+1..num_pops {
                        if d == b { continue; }
                        stats.push(FStatKind::F4(a, b, c, d));
                    }
                }
            }
        }
        stats
    }
}
