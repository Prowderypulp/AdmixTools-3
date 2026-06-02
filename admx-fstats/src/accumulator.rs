//! Per-SNP accumulation engine.
//!
//! Implements the §3.1 optimization: instead of evaluating O(numeg⁴) statistics
//! per SNP, we accumulate the numeg×numeg M₂ matrix and H vector via BLAS
//! rank-1 updates, then recover all statistics through a single sparse
//! contraction at the end.

use ndarray::{Array1, Array2};
use admx_core::types::DetailedPopCounts;

/// Sentinel `aax` value for a population whose frequency is estimable but whose
/// heterozygosity (and thus the unbiased p² correction) is not. Mirrors C's
/// `aaxadd = -999`; any f-statistic that adds it on a diagonal term becomes
/// `< -99` and is dropped (see `eval_fstat` and C dofstats `yy < -99`).
pub const AAX_INVALID: f64 = -999.0;

/// Accumulates sums for a single jackknife block.
#[derive(Debug, Clone)]
pub struct BlockAccumulator {
    /// Sum of $p_a \times p_b$ across SNPs.
    pub sum_prod: Array2<f64>,
    /// Sum of frequencies $p_a$ across SNPs.
    pub sum_p: Array1<f64>,
    /// Sum of heterozygosity corrections $aax_a$ across SNPs.
    pub sum_corr: Array1<f64>,
    /// Count of SNPs where each pair of populations was present.
    pub counts: Array2<f64>,
}

impl BlockAccumulator {
    pub fn new(num_pops: usize) -> Self {
        Self {
            sum_prod: Array2::zeros((num_pops, num_pops)),
            sum_p: Array1::zeros(num_pops),
            sum_corr: Array1::zeros(num_pops),
            counts: Array2::zeros((num_pops, num_pops)),
        }
    }

    /// Update with a single SNP's frequencies and mask.
    /// 
    /// `p` and `aax` should be 0.0 for missing populations.
    /// `mask` should be 1.0 for present, 0.0 for missing.
    pub fn add_snp(&mut self, p: &Array1<f64>, aax: &Array1<f64>, mask: &Array1<f64>) {
        let n = p.len();
        for i in 0..n {
            if mask[i] == 0.0 { continue; }
            self.sum_p[i] += p[i];
            self.sum_corr[i] += aax[i];
            for j in i..n {
                if mask[j] == 0.0 { continue; }
                let val_prod = p[i] * p[j];
                self.sum_prod[[i, j]] += val_prod;
                if i != j {
                    self.sum_prod[[j, i]] += val_prod;
                }
                
                self.counts[[i, j]] += 1.0;
                if i != j {
                    self.counts[[j, i]] += 1.0;
                }
            }
        }
    }
}

/// Estimates $p_a$ and $aax_a$ from detailed population counts.
/// 
/// Returns `(p, aax, is_valid)`.
pub fn estimate_pop_stats(
    counts: &DetailedPopCounts,
    inbreed: bool,
    is_haploid_snp: bool,
) -> (f64, f64, bool) {
    let pc = counts.to_pop_counts();
    let nref = pc.nref as f64;
    let nalt = pc.nalt as f64;
    let s = (nref + nalt) as f64;

    if s < 0.5 {
        return (0.0, 0.0, false);
    }

    let p = nalt / s;

    if is_haploid_snp {
        // Chromosome 24 (Y) or MT, or other haploid markers
        (p, 0.0, true)
    } else if inbreed {
        // C loadaa (qpsubs.c:4251) sets a valid frequency as soon as s >= 0.5
        // (≥1 sample), then leaves aaxadd = -999 if the heterozygosity estimate
        // can't be formed (s < 1.5 or yt < 0.001). The pop is still usable in
        // any f-statistic where it appears OFF the diagonal (only its frequency
        // is needed); fstatx/dofstats drop the stat (yy < -99) only when the
        // -999 aax lands on a diagonal term. So we return frequency-valid with
        // aax = AAX_INVALID rather than failing the whole pop.
        if s < 1.5 {
            return (p, AAX_INVALID, true);
        }

        let yf = (counts.f0 + counts.f1 + counts.f2) as f64;
        let ym = (counts.m0 + counts.m1) as f64;

        let yt = 4.0 * yf * (yf - 1.0) + 4.0 * yf * ym + ym * (ym - 1.0);
        if yt < 0.001 {
            return (p, AAX_INVALID, true);
        }

        let jhet = 4 * counts.f0 * counts.f2 
                 + 2 * counts.f2 * counts.f1 
                 + 2 * counts.f0 * counts.f1 
                 + counts.f1.saturating_sub(1) * counts.f1
                 + (2 * counts.f0 + counts.f1) * counts.m1
                 + (2 * counts.f2 + counts.f1) * counts.m0
                 + counts.m0 * counts.m1;
        
        let h1 = jhet as f64 / yt;
        let hest = h1;
        let aax = (p - hest) - p * p;
        (p, aax, true)
    } else {
        // Diploid, no inbreed
        if s < 1.5 {
            return (0.0, 0.0, false);
        }
        let hest = (nref * nalt) / (s * (s - 1.0));
        let aax = (p - hest) - p * p;
        (p, aax, true)
    }
}
