//! End-to-end `qpfstats` driver.
//!
//! Implements the C `dofstats` flow from `qpsubs.c:4665-5031`: per-SNP
//! evaluation of every f-statistic in the canonical enumeration
//! (all f2 + anchored f3 + canonical f4), per-fstat block-jackknife for
//! per-fstat sigmas, then a weighted least-squares consensus fit that
//! compresses the nfstats estimates onto the anchored-f3 basis.

use ndarray::{Array1, Array2};
use admx_core::error::{AdmxError, AdmxResult};
use admx_core::types::{Snp, Indiv, DetailedPopCounts, FStatKind, Sex};
use admx_core::blocks::assign_blocks;
use admx_core::jackknife::{wjackest, wjackvest};
use admx_core::linalg::{dgemm, pdinv};
use admx_core::constants::DOFSTATS_DIAG;
use admx_io::GenoReader;
use admx_io::codec::{G_MISSING, unpack_genotype};
use crate::accumulator::estimate_pop_stats;
use crate::basis::FBasis;
use rayon::prelude::*;

/// Per-block accumulators for the parallel SNP scan. Each jackknife block is
/// processed by exactly one thread over its SNP range in SNP order, so every
/// field is summed in the same order the serial loop would — bit-identical.
struct BlockAccum {
    btop: Vec<f64>,
    bbot: Vec<f64>,
    wjack: f64,
    pair_btop: Vec<f64>,
    pair_bbot: Vec<f64>,
    pair_bcnt: Vec<f64>,
    // Global hetrate partials, merged across blocks in block order on return.
    sum_p: Vec<f64>,
    sum_p2: Vec<f64>,
    n_p: Vec<f64>,
    het_seen: Vec<bool>,
}

impl BlockAccum {
    fn new(nfstats: usize, num_pops: usize) -> Self {
        BlockAccum {
            btop: vec![0.0; nfstats],
            bbot: vec![0.0; nfstats],
            wjack: 0.0,
            pair_btop: vec![0.0; num_pops * num_pops],
            pair_bbot: vec![0.0; num_pops * num_pops],
            pair_bcnt: vec![0.0; num_pops * num_pops],
            sum_p: vec![0.0; num_pops],
            sum_p2: vec![0.0; num_pops],
            n_p: vec![0.0; num_pops],
            het_seen: vec![false; num_pops],
        }
    }
}

/// Polymorphism test (C setwt: `total_nref > 0 && total_nalt > 0`) with early
/// termination. Mirrors the allele contributions of [`process_record`]'s count
/// loop exactly — the boolean it returns must match, or block boundaries (and
/// hence per-block jackknife sums) would diverge from the serial path.
#[inline]
fn record_polymorphic(rec: &[u8], snp: &Snp, indivs: &[Indiv], active_inds: &[(usize, usize)]) -> bool {
    let mut seen_ref = false;
    let mut seen_alt = false;
    for &(ind_i, _) in active_inds {
        let g = unpack_genotype(rec, ind_i);
        if g == G_MISSING { continue; }
        let is_haploid = match snp.chrom {
            23 => indivs[ind_i].sex == Sex::Male,
            24 | 90 => true,
            _ => false,
        };
        if is_haploid {
            if g == 0 { seen_ref = true; } else if g == 2 { seen_alt = true; }
        } else {
            match g {
                0 => seen_ref = true,
                1 => { seen_ref = true; seen_alt = true; }
                2 => seen_alt = true,
                _ => {}
            }
        }
        if seen_ref && seen_alt { return true; }
    }
    false
}

/// Accumulate one SNP record into the (block-local) accumulators. Byte-for-byte
/// the same arithmetic, in the same order, as the serial dofstats loop body.
/// `counts`/`p`/`aax`/`mask`/`het_valid` are caller-owned scratch reused across
/// SNPs (no per-SNP heap allocation). Returns true if the SNP is polymorphic.
#[inline]
#[allow(clippy::too_many_arguments)]
fn process_record(
    rec: &[u8],
    snp: &Snp,
    indivs: &[Indiv],
    active_inds: &[(usize, usize)],
    nfstats_kinds: &[FStatKind],
    num_pops: usize,
    inbreed: bool,
    allsnps: bool,
    acc: &mut BlockAccum,
    counts: &mut [DetailedPopCounts],
    p: &mut [f64],
    aax: &mut [f64],
    mask: &mut [bool],
    het_valid: &mut [bool],
) -> bool {
    for c in counts.iter_mut() { *c = DetailedPopCounts::default(); }
    let mut total_nref = 0.0_f64;
    let mut total_nalt = 0.0_f64;
    for &(ind_i, p_idx) in active_inds {
        let g = unpack_genotype(rec, ind_i);
        if g == G_MISSING { continue; }
        let is_haploid = match snp.chrom {
            23 => indivs[ind_i].sex == Sex::Male,
            24 | 90 => true,
            _ => false,
        };
        if is_haploid {
            if g == 0 { counts[p_idx].m0 += 1; total_nref += 1.0; }
            else if g == 2 { counts[p_idx].m1 += 1; total_nalt += 1.0; }
        } else if g == 0 {
            counts[p_idx].f0 += 1; total_nref += 2.0;
        } else if g == 1 {
            counts[p_idx].f1 += 1; total_nref += 1.0; total_nalt += 1.0;
            acc.het_seen[p_idx] = true;
        } else if g == 2 {
            counts[p_idx].f2 += 1; total_nalt += 2.0;
        }
    }
    if total_nref == 0.0 || total_nalt == 0.0 { return false; }

    let is_haploid_snp = snp.chrom == 24 || snp.chrom == 90;
    for p_idx in 0..num_pops {
        p[p_idx] = 0.0; aax[p_idx] = 0.0; mask[p_idx] = false; het_valid[p_idx] = false;
        let (pv, aaxv, valid) = estimate_pop_stats(&counts[p_idx], inbreed, is_haploid_snp);
        if valid {
            p[p_idx] = pv;
            aax[p_idx] = aaxv;
            mask[p_idx] = true;
            if aaxv > -99.0 { het_valid[p_idx] = true; }
        }
    }

    let all_present = het_valid.iter().all(|&m| m);
    if !allsnps && !all_present { return true; }

    acc.wjack += 1.0;
    for (j, &kind) in nfstats_kinds.iter().enumerate() {
        if let Some(yy) = eval_fstat(kind, p, aax, mask) {
            acc.btop[j] += yy;
            acc.bbot[j] += 1.0;
        }
    }
    for i in 0..num_pops {
        if !het_valid[i] { continue; }
        acc.sum_p[i] += p[i];
        acc.sum_p2[i] += p[i] * p[i] + aax[i];
        acc.n_p[i] += 1.0;
        for j in i + 1..num_pops {
            if !het_valid[j] { continue; }
            let en = p[i] * p[i] + aax[i] + p[j] * p[j] + aax[j] - 2.0 * p[i] * p[j];
            let hest_i = p[i] - (p[i] * p[i] + aax[i]);
            let hest_j = p[j] - (p[j] * p[j] + aax[j]);
            let ed = en + hest_i + hest_j;
            acc.pair_btop[i * num_pops + j] += en;
            acc.pair_bbot[i * num_pops + j] += ed;
            acc.pair_bcnt[i * num_pops + j] += 1.0;
        }
    }
    true
}

pub struct QpfstatsConfig {
    pub blgsize: f64,
    pub inbreed: bool,
    pub hires: bool,
    pub allsnps: bool,
    pub noxdata: bool,
    pub numchrom: i32,
    pub doscale: bool,
    pub anchor_pop: String,
}

impl Default for QpfstatsConfig {
    fn default() -> Self {
        Self {
            blgsize: 0.05,
            inbreed: false,
            hires: true,
            allsnps: false,
            noxdata: true,
            numchrom: 22,
            doscale: true,
            anchor_pop: String::new(),
        }
    }
}

/// Evaluate an f-statistic from per-population `p` and `aax` arrays.
///
/// Returns `None` if any population involved in the stat is missing.
/// Mirrors `fstatx` in `qpsubs.c` — including the `aaxadd` corrections
/// applied when index pairs collide.
#[inline]
fn eval_fstat(kind: FStatKind, p: &[f64], aax: &[f64], mask: &[bool]) -> Option<f64> {
    let (a, b, c, d) = match kind {
        FStatKind::F2(a, b) => (a, b, a, b),
        FStatKind::F3(anchor, x, y) => (anchor, x, anchor, y),
        FStatKind::F4(a, b, c, d) => (a, b, c, d),
    };

    if a == b || c == d { return Some(0.0); }
    // Participation requires only frequency validity (mask). aax may be the
    // AAX_INVALID sentinel for low-sample pops; it only matters on diagonal
    // terms, where it drives yy below -99 and the stat is dropped — exactly
    // C's fstatx + dofstats `yy < -99` behaviour.
    if !(mask[a] && mask[b] && mask[c] && mask[d]) { return None; }

    let mut yy = (p[a] - p[b]) * (p[c] - p[d]);
    if a == c { yy += aax[a]; }
    if b == d { yy += aax[b]; }
    if a == d { yy -= aax[a]; }
    if b == c { yy -= aax[b]; }
    if yy < -99.0 { return None; }
    Some(yy)
}

pub struct QpfstatsResult {
    pub means: Array1<f64>,
    pub covar: Array2<f64>,
    pub w3: Array1<f64>,
    pub wls_stderr: Array1<f64>,
    pub lambdascale: f64,
    pub num_snps: usize,
    pub valid_snps: usize,
    pub before_setwt_snps: usize,
    pub num_blocks: usize,
    pub pop_stats: Vec<(usize, f64, usize)>, // (index, hetrate, valid_snps)
    pub fst: Array2<f64>,
    pub f2: Array2<f64>,
    /// Raw per-fstat jackknife mean for each basis statistic (C `fsmean[basisfn[i]]`).
    pub fs_mean: Array1<f64>,
    /// Raw per-fstat jackknife sigma for each basis statistic (C `fssig[basisfn[i]]`).
    pub fs_sig: Array1<f64>,
    /// Number of per-fstat sigma adjustments applied (C `numadj`, "adjusted sigs").
    pub num_adjusted: usize,
}

pub fn run_qpfstats(
    geno_reader: &mut dyn GenoReader,
    snps: &[Snp],
    indivs: &[Indiv],
    pop_list: &[String],
    config: &QpfstatsConfig,
) -> AdmxResult<QpfstatsResult> {
    let num_pops = pop_list.len();
    if num_pops == 0 {
        return Err(AdmxError::Fatal("No populations provided".to_string()));
    }

    // Env-gated phase timing to stderr (ADMX_PROFILE=1). Does not affect output.
    let _prof = std::env::var("ADMX_PROFILE").is_ok();
    let mut _t = std::time::Instant::now();
    macro_rules! phase {
        ($name:expr) => {
            if _prof {
                eprintln!("[prof] {:<14} {:>8.3}s", $name, _t.elapsed().as_secs_f64());
                _t = std::time::Instant::now();
            }
        };
    }

    let pop_map: std::collections::HashMap<String, usize> = pop_list.iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), i))
        .collect();

    let ind_pop_indices: Vec<Option<usize>> = indivs.iter()
        .map(|ind| pop_map.get(&ind.egroup).cloned())
        .collect();

    // Block jackknife membership must be computed over exactly the SNPs that
    // are actually used (non-ignored AND polymorphic), matching C, which runs
    // setblocks on the post-rmsnps `xsnplist` (qpfstats.c:480-525). Assigning
    // over all SNPs would let monomorphic SNPs shift block anchors and perturb
    // the per-block jackknife variance. We use assign_blocks over the full list
    // only as a safe upper bound for allocation (used ⊆ all ⇒ blocks ≤ this),
    // and recompute the real block index inside the loop below.
    let block_assignments = assign_blocks(snps, config.blgsize);
    let num_blocks = (block_assignments.iter().max().cloned().unwrap_or(0) + 2) as usize;
    if num_blocks <= 1 {
        return Err(AdmxError::Fatal("No jackknife blocks found".to_string()));
    }

    // Anchor must be index 0 after CLI reorder.
    let anchor_idx = 0usize;
    let basis = FBasis::new(anchor_idx, num_pops);
    let basis_stats = basis.stats();
    let nbasis = basis_stats.len();

    // Full canonical enumeration.
    let nfstats_kinds = FBasis::full_stats(num_pops);
    let nfstats = nfstats_kinds.len();

    // Build dense nfstats × nbasis coefficient matrix (row-major).
    let mut fbcoeffs = vec![0.0_f64; nfstats * nbasis];
    for (i, &kind) in nfstats_kinds.iter().enumerate() {
        for (j, c) in basis.coefficients(kind) {
            fbcoeffs[i * nbasis + j] = c;
        }
    }

    // Per-block per-fstat top/bot accumulators and SNP counts.
    let mut btop: Vec<Vec<f64>> = vec![vec![0.0; nfstats]; num_blocks];
    let mut bbot: Vec<Vec<f64>> = vec![vec![0.0; nfstats]; num_blocks];
    let mut wjack: Vec<f64> = vec![0.0; num_blocks];

    // Per-pair accumulators used for lambdascale FST regression (dofstnumx).
    // pair_btop = Σ numerator (en), pair_bbot = Σ denominator (ed) for the
    // classic FST ratio, pair_bcnt = count of valid SNPs for that pair (C's
    // `bot[k] += 1.0` in fstmode=NO — the denominator of `fstest`).
    let mut pair_btop = vec![vec![0.0f64; num_pops * num_pops]; num_blocks];
    let mut pair_bbot = vec![vec![0.0f64; num_pops * num_pops]; num_blocks];
    let mut pair_bcnt = vec![vec![0.0f64; num_pops * num_pops]; num_blocks];

    // Per-pair accumulators for hetrate.
    let mut n_p = vec![0.0_f64; num_pops];
    let mut sum_p = vec![0.0_f64; num_pops];
    let mut sum_p2_plus_aax = vec![0.0_f64; num_pops];

    // C `counthets` (qpsubs.c:4604): has this pop ever shown an actual
    // heterozygous genotype call (g == 1) on a used SNP? This is distinct from
    // the model-based hetrate — pseudo-haploid samples (genotypes only 0/2)
    // have a positive hetrate but no het CALLS, so `het_seen` is false. Drives
    // the `hashets == 0` sigma adjustment (qpsubs.c:4864-4882).
    let mut het_seen = vec![false; num_pops];

    let mut before_setwt_snps = 0;
    let mut valid_snps_total = 0;

    // In-loop jackknife block tracking over USED SNPs (C setblocks semantics:
    // new block when chrom changes or genpos - block_start >= blgsize, anchored
    // at the block's first used SNP).
    let mut cur_block: i32 = -1;
    let mut blk_start_pos = 0.0_f64;
    let mut blk_chrom: i32 = -1;

    // Compact (individual index, population index) list — only individuals that
    // belong to a requested population. The genotype file holds all individuals
    // but typically only a small fraction are in the poplist, so iterating this
    // instead of every individual avoids hundreds of millions of no-op checks.
    let active_inds: Vec<(usize, usize)> = ind_pop_indices.iter().enumerate()
        .filter_map(|(i, p)| p.map(|pp| (i, pp)))
        .collect();

    if let Some((gbytes, off0, stride, reclen)) = geno_reader.random_access() {
        // ---- Parallel SNP-major scan (mmap / in-memory backed readers) ----
        // Output is bit-identical to the serial path: per-block accumulators are
        // summed in SNP order within each block, and blocks are disjoint SNP
        // ranges processed independently, so no cross-block float regrouping
        // occurs for the jackknife inputs. (The hetrate `sum_p`/`n_p` globals are
        // merged across blocks in block order; verified against the baseline.)
        let nsnp = snps.len();
        before_setwt_snps = nsnp;
        valid_snps_total = (0..nsnp)
            .filter(|&i| !snps[i].ignore && !(config.noxdata && snps[i].chrom > config.numchrom))
            .count();

        // Pass 1 (parallel): polymorphism flag per SNP.
        let poly: Vec<bool> = (0..nsnp).into_par_iter().map(|i| {
            let snp = &snps[i];
            if snp.ignore { return false; }
            if config.noxdata && snp.chrom > config.numchrom { return false; }
            let base = off0 + i * stride;
            record_polymorphic(&gbytes[base..base + reclen], snp, indivs, &active_inds)
        }).collect();

        // Block assignment over polymorphic SNPs (metadata only) — replicates the
        // serial state machine exactly. Each block is a contiguous SNP-index
        // range [lo, hi); non-poly SNPs that fall inside a range are re-checked
        // and skipped by `process_record` in pass 2.
        let mut ranges: Vec<(usize, usize)> = Vec::new();
        let mut blk_chrom: i32 = -1;
        let mut blk_start_pos = 0.0_f64;
        let mut lo = 0usize;
        for i in 0..nsnp {
            if !poly[i] { continue; }
            let snp = &snps[i];
            if snp.chrom != blk_chrom || (snp.genpos - blk_start_pos) >= config.blgsize {
                if cur_block >= 0 { ranges.push((lo, i)); }
                cur_block += 1;
                blk_chrom = snp.chrom;
                blk_start_pos = snp.genpos;
                lo = i;
            }
        }
        if cur_block >= 0 { ranges.push((lo, nsnp)); }
        debug_assert!(ranges.len() <= num_blocks, "ranges {} > alloc {}", ranges.len(), num_blocks);

        // Pass 2 (parallel over blocks). rayon's indexed collect preserves order,
        // so block_results[g] corresponds to ranges[g] (jackknife block g).
        let block_results: Vec<BlockAccum> = ranges.par_iter().map(|&(lo, hi)| {
            let mut acc = BlockAccum::new(nfstats, num_pops);
            let mut counts = vec![DetailedPopCounts::default(); num_pops];
            let mut p = vec![0.0_f64; num_pops];
            let mut aax = vec![0.0_f64; num_pops];
            let mut mask = vec![false; num_pops];
            let mut het_valid = vec![false; num_pops];
            for i in lo..hi {
                let snp = &snps[i];
                if snp.ignore { continue; }
                if config.noxdata && snp.chrom > config.numchrom { continue; }
                let base = off0 + i * stride;
                process_record(
                    &gbytes[base..base + reclen], snp, indivs, &active_inds,
                    &nfstats_kinds, num_pops, config.inbreed, config.allsnps,
                    &mut acc, &mut counts, &mut p, &mut aax, &mut mask, &mut het_valid,
                );
            }
            acc
        }).collect();

        // Merge in block order: per-block arrays placed by index (no summing);
        // hetrate globals summed block-by-block; het_seen OR-ed.
        for (g, acc) in block_results.into_iter().enumerate() {
            btop[g] = acc.btop;
            bbot[g] = acc.bbot;
            wjack[g] = acc.wjack;
            pair_btop[g] = acc.pair_btop;
            pair_bbot[g] = acc.pair_bbot;
            pair_bcnt[g] = acc.pair_bcnt;
            for pp in 0..num_pops {
                if acc.het_seen[pp] { het_seen[pp] = true; }
                sum_p[pp] += acc.sum_p[pp];
                sum_p2_plus_aax[pp] += acc.sum_p2[pp];
                n_p[pp] += acc.n_p[pp];
            }
        }
    } else {
        // ---- Serial fallback (streaming/text readers, e.g. EIGENSTRAT) ----
        let mut record = vec![0u8; geno_reader.record_bytes()];

        for (_snp_i, snp) in snps.iter().enumerate() {
            before_setwt_snps += 1;
            let have_row = geno_reader.read_record(&mut record).map_err(AdmxError::Io)?;
            if !have_row { break; }

            if snp.ignore { continue; }

            if config.noxdata && snp.chrom > config.numchrom { continue; }
            valid_snps_total += 1;

            // Per-SNP detailed counts.
            let mut counts = vec![DetailedPopCounts::default(); num_pops];
            let mut total_nref = 0.0_f64;
            let mut total_nalt = 0.0_f64;
            for (ind_i, &pop_idx) in ind_pop_indices.iter().enumerate() {
                if let Some(p_idx) = pop_idx {
                    let g = unpack_genotype(&record, ind_i);
                    if g == G_MISSING { continue; }
                    let is_haploid = match snp.chrom {
                        23 => indivs[ind_i].sex == Sex::Male,
                        24 | 90 => true,
                        _ => false,
                    };
                    if is_haploid {
                        if g == 0 { counts[p_idx].m0 += 1; total_nref += 1.0; }
                        else if g == 2 { counts[p_idx].m1 += 1; total_nalt += 1.0; }
                    } else if g == 0 {
                        counts[p_idx].f0 += 1; total_nref += 2.0;
                    } else if g == 1 {
                        counts[p_idx].f1 += 1; total_nref += 1.0; total_nalt += 1.0;
                        het_seen[p_idx] = true; // observed a real het call (C counthets)
                    } else if g == 2 {
                        counts[p_idx].f2 += 1; total_nalt += 2.0;
                    }
                }
            }

            // Polymorphism filter (match C setwt).
            if total_nref == 0.0 || total_nalt == 0.0 { continue; }

            // Assign this used SNP to a jackknife block (C setblocks over xsnplist).
            if snp.chrom != blk_chrom || (snp.genpos - blk_start_pos) >= config.blgsize {
                cur_block += 1;
                blk_chrom = snp.chrom;
                blk_start_pos = snp.genpos;
            }
            let block_idx = cur_block as usize;
            // Invariant: used SNPs ⊆ all SNPs, so the in-loop block count never
            // exceeds the up-front allocation bound. Fail loudly if that breaks.
            debug_assert!(block_idx < num_blocks, "block_idx {} >= alloc {}", block_idx, num_blocks);

            // Per-population p, aax, masks. `mask` = frequency valid (pop can appear
            // in any f-statistic). `het_valid` = aax (heterozygosity) also valid
            // (pop can appear on a diagonal term and in the FST/hetrate sums).
            let mut p = vec![0.0_f64; num_pops];
            let mut aax = vec![0.0_f64; num_pops];
            let mut mask = vec![false; num_pops];
            let mut het_valid = vec![false; num_pops];
            let is_haploid_snp = snp.chrom == 24 || snp.chrom == 90;
            for p_idx in 0..num_pops {
                let (pv, aaxv, valid) = estimate_pop_stats(&counts[p_idx], config.inbreed, is_haploid_snp);
                if valid {
                    p[p_idx] = pv;
                    aax[p_idx] = aaxv;
                    mask[p_idx] = true;
                    if aaxv > -99.0 { het_valid[p_idx] = true; }
                }
            }

            // C dofstats (qpsubs.c:4740-4775) has no anchor-presence gate: every
            // non-ignored SNP is counted (wjack++) and each f-statistic is skipped
            // individually when its inputs are missing. Gating on the anchor here
            // wrongly restricted the all-pairs FST accumulation (and thus
            // lambdascale) to anchor-present SNPs. The all_present check below
            // still implements allsnps:NO (it subsumes anchor presence).
            let all_present = het_valid.iter().all(|&m| m);
            if !config.allsnps && !all_present { continue; }

            wjack[block_idx] += 1.0;

            // Per-fstat accumulation.
            let top = &mut btop[block_idx];
            let bot = &mut bbot[block_idx];
            for (j, &kind) in nfstats_kinds.iter().enumerate() {
                if let Some(yy) = eval_fstat(kind, &p, &aax, &mask) {
                    top[j] += yy;
                    bot[j] += 1.0;
                }
            }

            // dofstnumx accumulation.
            let ptop = &mut pair_btop[block_idx];
            let pbot = &mut pair_bbot[block_idx];
            let pcnt = &mut pair_bcnt[block_idx];
            for i in 0..num_pops {
                if !het_valid[i] { continue; }
                sum_p[i] += p[i];
                sum_p2_plus_aax[i] += p[i] * p[i] + aax[i];
                n_p[i] += 1.0;
                for j in i + 1..num_pops {
                    if !het_valid[j] { continue; }
                    let en = p[i] * p[i] + aax[i] + p[j] * p[j] + aax[j] - 2.0 * p[i] * p[j];
                    let hest_i = p[i] - (p[i] * p[i] + aax[i]);
                    let hest_j = p[j] - (p[j] * p[j] + aax[j]);
                    let ed = en + hest_i + hest_j;
                    ptop[i * num_pops + j] += en;
                    pbot[i * num_pops + j] += ed;
                    pcnt[i * num_pops + j] += 1.0;
                }
            }
        }
    }

    phase!("snp_loop");

    let active_blocks: Vec<usize> = (0..num_blocks).filter(|&k| wjack[k] > 0.0).collect();
    let active_wjack: Vec<f64> = active_blocks.iter().map(|&k| wjack[k]).collect();

    // Global per-pair sums.
    let mut g_pair_top = vec![0.0f64; num_pops * num_pops];
    let mut g_pair_bot = vec![0.0f64; num_pops * num_pops];
    let mut g_pair_count = vec![0.0f64; num_pops * num_pops];
    for k in 0..num_blocks {
        for idx in 0..num_pops * num_pops {
            g_pair_top[idx] += pair_btop[k][idx];
            g_pair_bot[idx] += pair_bbot[k][idx];
            g_pair_count[idx] += pair_bcnt[k][idx];
        }
    }

    let mut pair_fst = Array2::zeros((num_pops, num_pops));
    let mut pair_fstest = Array2::zeros((num_pops, num_pops));
    let mut pair_fstsig = Array2::zeros((num_pops, num_pops));

    for i in 0..num_pops {
        for j in i + 1..num_pops {
            let idx = i * num_pops + j;
            // fstest (fstmode=NO) is the per-pair numerator mean: denominator
            // is the count of valid SNPs for THIS pair (C dofstnumx), not the
            // global SNP count. fst stays the classic ratio Σen/Σed.
            let mean_en = g_pair_top[idx] / (g_pair_count[idx] + 1e-10);
            let mean_fst = g_pair_top[idx] / (g_pair_bot[idx] + 1e-10);
            pair_fst[[i, j]] = mean_fst;
            pair_fst[[j, i]] = mean_fst;

            let mut jmeans = Vec::with_capacity(active_blocks.len());
            for &k in &active_blocks {
                // delete-block numerator mean, per-pair count denominator
                let tk = g_pair_top[idx] - pair_btop[k][idx];
                let bk = g_pair_count[idx] - pair_bcnt[k][idx] + 1e-10;
                jmeans.push(tk / bk);
            }
            let (jest, jsig) = wjackest(mean_en, &jmeans, &active_wjack);
            pair_fstest[[i, j]] = jest;
            pair_fstest[[j, i]] = jest;
            pair_fstsig[[i, j]] = jsig;
            pair_fstsig[[j, i]] = jsig;
        }
    }

    let mut lambdascale = 1.0_f64;
    if config.doscale {
        let mut y1 = 0.0;
        let mut y2 = 0.0;
        for i in 0..num_pops {
            for j in 0..num_pops {
                let sig = pair_fstsig[[i, j]] + 1e-10;
                let w1 = pair_fst[[i, j]] / sig;
                let w2 = pair_fstest[[i, j]] / sig;
                y1 += w1 * w2;
                y2 += w2 * w2;
            }
        }
        if y2 > 0.0 { lambdascale = y1 / y2; }
    }

    let mut pair_f2 = pair_fstest.clone();

    // Global per-fstat sums.
    let mut gtop = vec![0.0_f64; nfstats];
    let mut gbot = vec![0.0_f64; nfstats];
    for k in 0..num_blocks {
        for j in 0..nfstats {
            gtop[j] += btop[k][j];
            gbot[j] += bbot[k][j];
        }
    }

    // Apply lambdascale uniformly to every accumulated sum.
    if lambdascale != 1.0 {
        for j in 0..nfstats {
            gtop[j] *= lambdascale;
            for k in 0..num_blocks { btop[k][j] *= lambdascale; }
        }
        for i in 0..num_pops {
            for j in 0..num_pops {
                pair_f2[[i, j]] *= lambdascale;
            }
        }
    }

    // Per-block delete-one means (overrides btop[k] — matches C).
    for k in 0..num_blocks {
        for j in 0..nfstats {
            let wt = gtop[j] - btop[k][j];
            let wb = gbot[j] - bbot[k][j] + 1.0e-12;
            btop[k][j] = wt / wb;
        }
    }

    // Global per-fstat mean (overrides gtop — matches C).
    for j in 0..nfstats {
        gtop[j] /= gbot[j] + 1.0e-12;
    }

    // Per-fstat jackknife. Capture both the jackknife-corrected mean (`jest`)
    // and sigma so the basis report can print the raw per-fstat values (C
    // `fsmean[t]` / `fssig[t]` in qpsubs.c:dofstats).
    let mut jsig = vec![0.0_f64; nfstats];
    let mut jest = vec![0.0_f64; nfstats];
    for j in 0..nfstats {
        if gbot[j] < 0.001 {
            jest[j] = 0.0;
            jsig[j] = 1.0e6;
            continue;
        }
        let jmean_j: Vec<f64> = active_blocks.iter().map(|&k| btop[k][j]).collect();
        let (est, sig) = wjackest(gtop[j], &jmean_j, &active_wjack);
        jest[j] = est;
        jsig[j] = sig;
    }
    // C `hashets == 0` sigma adjustment (qpsubs.c:4864-4882). For each fstat,
    // if a NON-inbred population with no observed het calls appears on a
    // diagonal term (i.e. its index occurs >1 time among the stat's 4 indices),
    // inflate that fstat's sigma: jsig = sqrt(jsig^2 + 100). This massively
    // downweights f2(O,a) self-terms and O-anchored f3 for zero-het pops in the
    // WLS consensus fit, so their fitted f2 collapses toward 0. Gated on the
    // global inbreed flag: under inbreed:YES every pop is inbred and C `continue`s
    // for all of them, so the block is a no-op (matching the validated path).
    // The 4 indices per stat mirror `eval_fstat`'s (a,b,c,d) expansion.
    let mut num_adjusted = 0usize;
    if !config.inbreed {
        for j in 0..nfstats {
            // Mirror C: the adjustment lives after the `gbot[j] < .001` zeroing
            // guard, so fstats with no valid SNPs (jsig already 1e6) are skipped
            // and don't count toward numadj.
            if gbot[j] < 0.001 { continue; }
            let (i0, i1, i2, i3) = match nfstats_kinds[j] {
                FStatKind::F2(a, b) => (a, b, a, b),
                FStatKind::F3(anchor, x, y) => (anchor, x, anchor, y),
                FStatKind::F4(a, b, c, d) => (a, b, c, d),
            };
            let mut dd = vec![0u32; num_pops];
            dd[i0] += 1; dd[i1] += 1; dd[i2] += 1; dd[i3] += 1;
            for t in 0..num_pops {
                if dd[t] > 1 && !het_seen[t] {
                    jsig[j] = (jsig[j] * jsig[j] + 100.0).sqrt();
                    num_adjusted += 1;
                }
            }
        }
    }

    // Flatten: pretend we have gbot observations but a N(0,1) prior on f —
    // stabilizes tiny samples. C flattens BOTH the mean and the sigma
    // (qpsubs.c:4899-4906), not just the sigma.
    for j in 0..nfstats {
        let y1 = gbot[j] + 1.0;
        jest[j] = jest[j] * gbot[j] / y1;
        let y2 = jsig[j] * jsig[j] * gbot[j] + 1.0;
        jsig[j] = (y2 / y1).sqrt();
    }
    for j in 0..nfstats { jsig[j] += 1.0e-12; }

    phase!("jackknife");

    // Weighted design matrix wfb[i, j] = fbcoeffs[i, j] / jsig[i]  (row-major).
    let mut wfb = vec![0.0_f64; nfstats * nbasis];
    for i in 0..nfstats {
        let inv = 1.0 / jsig[i];
        for j in 0..nbasis {
            wfb[i * nbasis + j] = fbcoeffs[i * nbasis + j] * inv;
        }
    }

    // wco = wfbᵀ · wfb  (nbasis × nbasis).
    let mut wco = vec![0.0_f64; nbasis * nbasis];
    dgemm(
        b'N', b'T',
        nbasis, nbasis, nfstats,
        1.0,
        &wfb, nbasis,
        &wfb, nbasis,
        0.0,
        &mut wco, nbasis,
    );

    // Regularize: diag += DOFSTATS_DIAG * trace/nbasis, then + DOFSTATS_DIAG.
    let trace: f64 = (0..nbasis).map(|i| wco[i * nbasis + i]).sum();
    let ridge = DOFSTATS_DIAG * trace / (nbasis as f64);
    for i in 0..nbasis {
        wco[i * nbasis + i] += ridge + DOFSTATS_DIAG;
    }

    // Invert via Cholesky (pdinv).
    let mut wco_buf = wco.clone();
    pdinv(nbasis, &mut wco_buf)
        .map_err(|info| AdmxError::Fatal(format!("pdinv failed on WLS normal equations: info={}", info)))?;
    let wcoinv = wco_buf; 

    // Helper: compute nbasis-vector x = wcoinv · (wfbᵀ · (top_vec / jsig)).
    let wls_solve = |top_vec: &[f64]| -> Vec<f64> {
        let mut w2 = vec![0.0_f64; nfstats];
        for i in 0..nfstats { w2[i] = top_vec[i] / jsig[i]; }
        // w1 = w2 · wfb  → length nbasis.
        let mut w1 = vec![0.0_f64; nbasis];
        for i in 0..nfstats {
            let c = w2[i];
            for j in 0..nbasis {
                w1[j] += c * wfb[i * nbasis + j];
            }
        }
        // out = wcoinv · w1  (symmetric matmul).
        let mut out = vec![0.0_f64; nbasis];
        for i in 0..nbasis {
            let mut s = 0.0;
            for j in 0..nbasis {
                s += wcoinv[i * nbasis + j] * w1[j];
            }
            out[i] = s;
        }
        out
    };

    let w3 = wls_solve(&gtop);
    let mut vjmean: Vec<Vec<f64>> = Vec::with_capacity(num_blocks);
    for k in 0..num_blocks {
        vjmean.push(wls_solve(&btop[k]));
    }

    phase!("wls_solve");

    // Drop empty blocks
    let filtered_jmean: Vec<Vec<f64>> = active_blocks.iter().map(|&k| vjmean[k].clone()).collect();

    let (jackest, covar_flat) = wjackvest(nbasis, &w3, &filtered_jmean, &active_wjack)?;
    
    let pop_stats: Vec<(usize, f64, usize)> = (0..num_pops).map(|i| {
        let hetrate = if n_p[i] > 0.0 {
            let p_mean = sum_p[i] / n_p[i];
            let m_ii = sum_p2_plus_aax[i] / n_p[i];
            p_mean - m_ii
        } else {
            0.0
        };
        (i, hetrate, n_p[i] as usize)
    }).collect();

    let wls_stderr = Array1::from_vec((0..nbasis).map(|i| wcoinv[i * nbasis + i].sqrt()).collect());

    // Raw per-fstat mean/sigma for each basis statistic (C `basisfn[i]` →
    // `fsmean[t]` / `fssig[t]`). A diagonal basis f3(anchor; a, a) is the
    // f2(anchor, a) entry of the full enumeration; an off-diagonal one maps
    // straight to its f3 entry.
    let mut fs_mean = vec![0.0_f64; nbasis];
    let mut fs_sig = vec![0.0_f64; nbasis];
    for (i, &kind) in basis_stats.iter().enumerate() {
        if let FStatKind::F3(c, a, b) = kind {
            let target = if a == b {
                let (lo, hi) = if c < a { (c, a) } else { (a, c) };
                FStatKind::F2(lo, hi)
            } else {
                FStatKind::F3(c, a, b)
            };
            let t = nfstats_kinds.iter().position(|&k| k == target)
                .expect("basis statistic missing from full enumeration");
            fs_mean[i] = jest[t];
            fs_sig[i] = jsig[t];
        }
    }

    Ok(QpfstatsResult {
        means: Array1::from_vec(jackest),
        covar: Array2::from_shape_vec((nbasis, nbasis), covar_flat).unwrap(),
        w3: Array1::from_vec(w3),
        wls_stderr,
        lambdascale,
        num_snps: wjack.iter().sum::<f64>() as usize,
        valid_snps: valid_snps_total,
        before_setwt_snps,
        num_blocks: (cur_block + 1).max(0) as usize,
        pop_stats,
        fst: pair_fst,
        f2: pair_f2,
        fs_mean: Array1::from_vec(fs_mean),
        fs_sig: Array1::from_vec(fs_sig),
        num_adjusted,
    })
}

#[cfg(test)]
mod parallel_equiv_tests {
    use super::*;
    use admx_core::types::{Snp, Indiv, Sex};
    use admx_io::{GenoReader, Layout};
    use admx_io::codec::g_to_2bit;

    /// In-memory PACKEDANCESTRYMAP-style reader: contiguous canonical records,
    /// exposes `random_access` so `run_qpfstats` takes the parallel path.
    struct MemReader { bytes: Vec<u8>, nind: usize, nsnp: usize, rec_bytes: usize, next: usize }
    impl GenoReader for MemReader {
        fn nind(&self) -> usize { self.nind }
        fn nsnp(&self) -> usize { self.nsnp }
        fn layout(&self) -> Layout { Layout::SnpMajor }
        fn record_bytes(&self) -> usize { self.rec_bytes }
        fn read_record(&mut self, dst: &mut [u8]) -> std::io::Result<bool> {
            if self.next >= self.nsnp { return Ok(false); }
            let s = self.next * self.rec_bytes;
            dst.copy_from_slice(&self.bytes[s..s + self.rec_bytes]);
            self.next += 1;
            Ok(true)
        }
        fn random_access(&self) -> Option<(&[u8], usize, usize, usize)> {
            Some((&self.bytes, 0, self.rec_bytes, self.rec_bytes))
        }
    }

    /// Wrapper that hides `random_access` (returns `None`) so `run_qpfstats`
    /// takes the serial fallback path — the reference for the equivalence check.
    struct ForceSerial(MemReader);
    impl GenoReader for ForceSerial {
        fn nind(&self) -> usize { self.0.nind() }
        fn nsnp(&self) -> usize { self.0.nsnp() }
        fn layout(&self) -> Layout { Layout::SnpMajor }
        fn record_bytes(&self) -> usize { self.0.record_bytes() }
        fn read_record(&mut self, dst: &mut [u8]) -> std::io::Result<bool> { self.0.read_record(dst) }
        // random_access stays None (trait default) → serial path.
    }

    fn pack_row(gs: &[u8]) -> Vec<u8> {
        let rec = (gs.len() * 2 + 7) / 8;
        let mut bytes = vec![0u8; rec];
        for (i, &g) in gs.iter().enumerate() {
            let shift = 6 - 2 * (i % 4);
            bytes[i / 4] |= g_to_2bit(g) << shift;
        }
        bytes
    }

    fn indiv(id: &str, egroup: &str) -> Indiv {
        Indiv { id: id.into(), egroup: egroup.into(), sex: Sex::Unknown,
                idnum: 0, affstatus: 1, ignore: false, gkode: 0 }
    }

    fn snp(id: &str, chrom: i32, genpos: f64) -> Snp {
        Snp { id: id.into(), chrom, cchrom: chrom.to_string(), genpos, physpos: genpos * 1e8,
              alleles: ['A', 'C'], ignore: false, tagnumber: 0, weight: 1.0 }
    }

    /// Build a multi-block fixture: 8 individuals across 4 pops, ~60 SNPs over
    /// two chromosomes with genpos crossing several `blgsize` boundaries, a mix
    /// of polymorphic / monomorphic / partially-missing rows.
    fn fixture() -> (Vec<Indiv>, Vec<Snp>, Vec<String>, MemReader) {
        let indivs = vec![
            indiv("a0", "P0"), indiv("a1", "P0"),
            indiv("b0", "P1"), indiv("b1", "P1"),
            indiv("c0", "P2"), indiv("c1", "P2"),
            indiv("d0", "P3"), indiv("d1", "P3"),
        ];
        let pop_list = vec!["P0".to_string(), "P1".into(), "P2".into(), "P3".into()];

        let nind = indivs.len();
        let mut snps = Vec::new();
        let mut bytes = Vec::new();
        // Deterministic pseudo-genotype generator (no RNG dependency).
        let mut state = 0x1234_5678u32;
        let mut nextg = || { state = state.wrapping_mul(1103515245).wrapping_add(12345); (state >> 16) & 0xffff };
        let nsnp = 64usize;
        for s in 0..nsnp {
            let chrom = if s < 32 { 1 } else { 2 };
            // genpos rises 0.012 per SNP → with blgsize 0.05 a new block ~every 4-5 SNPs.
            let genpos = 0.012 * (s as f64);
            snps.push(snp(&format!("rs{s}"), chrom, genpos));
            let mut row = vec![0u8; nind];
            for (i, gi) in row.iter_mut().enumerate() {
                let r = nextg();
                *gi = match s % 7 {
                    0 => 0,                              // monomorphic ref (filtered out)
                    1 if i >= 6 => G_MISSING,            // P3 missing → exercises allsnps:NO
                    _ => match r % 4 { 0 => 0, 1 => 1, 2 => 2, _ => if s % 5 == 0 { G_MISSING } else { 2 } },
                };
            }
            bytes.extend_from_slice(&pack_row(&row));
        }
        let rec_bytes = (nind * 2 + 7) / 8;
        let reader = MemReader { bytes, nind, nsnp, rec_bytes, next: 0 };
        (indivs, snps, pop_list, reader)
    }

    fn cfg(allsnps: bool) -> QpfstatsConfig {
        QpfstatsConfig { blgsize: 0.05, inbreed: false, hires: true, allsnps,
                         noxdata: true, numchrom: 22, doscale: true, anchor_pop: "P0".into() }
    }

    fn assert_equiv(allsnps: bool) {
        let (indivs, snps, pops, mem) = fixture();
        let mut par = mem;            // parallel path (random_access = Some)
        let mut ser = ForceSerial(fixture().3); // serial path (random_access = None)

        let rp = run_qpfstats(&mut par, &snps, &indivs, &pops, &cfg(allsnps)).unwrap();
        let rs = run_qpfstats(&mut ser, &snps, &indivs, &pops, &cfg(allsnps)).unwrap();

        // Counts and block structure: exact.
        assert_eq!(rp.num_snps, rs.num_snps, "num_snps (allsnps={allsnps})");
        assert_eq!(rp.valid_snps, rs.valid_snps, "valid_snps");
        assert_eq!(rp.before_setwt_snps, rs.before_setwt_snps, "before_setwt_snps");
        assert_eq!(rp.num_blocks, rs.num_blocks, "num_blocks");
        assert_eq!(rp.num_adjusted, rs.num_adjusted, "num_adjusted");
        assert!(rp.num_blocks > 4, "fixture must span several blocks to exercise parallelism, got {}", rp.num_blocks);

        // btop-derived outputs are summed in identical per-block SNP order → bit-identical.
        let bits = |a: &[f64]| a.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(rp.means.as_slice().unwrap()), bits(rs.means.as_slice().unwrap()), "means bits");
        assert_eq!(bits(rp.covar.as_slice().unwrap()), bits(rs.covar.as_slice().unwrap()), "covar bits");
        assert_eq!(bits(rp.w3.as_slice().unwrap()), bits(rs.w3.as_slice().unwrap()), "w3 bits");
        assert_eq!(bits(rp.fs_mean.as_slice().unwrap()), bits(rs.fs_mean.as_slice().unwrap()), "fs_mean bits");
        assert_eq!(bits(rp.fs_sig.as_slice().unwrap()), bits(rs.fs_sig.as_slice().unwrap()), "fs_sig bits");
        assert_eq!(rp.lambdascale.to_bits(), rs.lambdascale.to_bits(), "lambdascale bits");
        assert_eq!(bits(rp.fst.as_slice().unwrap()), bits(rs.fst.as_slice().unwrap()), "fst bits");
        assert_eq!(bits(rp.f2.as_slice().unwrap()), bits(rs.f2.as_slice().unwrap()), "f2 bits");

        // hetrate (sum_p) is merged per-block, so it is equal only up to fp reorder
        // (this is the one quantity that is NOT bit-identical, and stays far below
        // the printed precision). Valid-SNP counts and pop indices are exact.
        for ((ip, hp, np), (is, hs, ns)) in rp.pop_stats.iter().zip(rs.pop_stats.iter()) {
            assert_eq!(ip, is); assert_eq!(np, ns);
            assert!((hp - hs).abs() < 1e-9, "hetrate reorder too large: {hp} vs {hs}");
        }
    }

    #[test]
    fn parallel_matches_serial_allsnps_yes() { assert_equiv(true); }

    #[test]
    fn parallel_matches_serial_allsnps_no() { assert_equiv(false); }
}
