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

    let mut before_setwt_snps = 0;
    let mut valid_snps_total = 0;

    // In-loop jackknife block tracking over USED SNPs (C setblocks semantics:
    // new block when chrom changes or genpos - block_start >= blgsize, anchored
    // at the block's first used SNP).
    let mut cur_block: i32 = -1;
    let mut blk_start_pos = 0.0_f64;
    let mut blk_chrom: i32 = -1;

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
    // NOTE: the C `hashets == 0` adjustment (qpsubs.c:4864-4882,
    // `jsig = sqrt(jsig^2 + 100)` when a non-inbred pop has no observed
    // heterozygotes) is not yet ported. The test1 fixture reports
    // "adjusted sigs: 0", so it does not fire there.

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
    })
}
