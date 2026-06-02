# Changelog — AdmixTools Rust Port

## [2026-06-02] — Fidelity pass: qpfstats/qpWave/qpAdm vs C binaries

Validated against the installed C AdmixTools (`/home/drtex/AdmixTools/bin`) on both
the small golden fixtures and a real 1.23M-SNP × 893-indiv dataset. All golden tests
pass; qpfstats now matches C to ~5-6 decimals under `inbreed:YES`/`allsnps:YES`.

### Fixed — qpfstats f-stat engine
- **Jackknife block assignment over the wrong SNP set (P0, real-data only).** Blocks were
  assigned over all SNPs up front, but C assigns over only the used (non-ignored +
  polymorphic) SNPs — `setblocks` on the post-`rmsnps` `xsnplist` (qpfstats.c:480-525).
  Monomorphic SNPs shifted block anchors, so per-block membership differed (713 vs 718
  blocks): the jackknife *mean* matched but the *variance* (basis covariance) diverged
  ~1-2%, which qpWave/qpAdm then amplified through the covariance inversion. Fixed by
  computing the block index inside the SNP loop, over used SNPs only
  (`admx-fstats/src/driver.rs`). After this, qpfstats matches C to 6 decimals (means AND
  sigmas) and the qpWave/qpAdm genotype paths match the C binaries.
- **Frequency- vs heterozygosity-validity (P0, real-data only).** A low-sample inbred
  population has an estimable frequency but no aax; C keeps it usable in any f-statistic
  where it appears off-diagonal and drops the statistic (`yy < -99`) only when the
  sentinel aax lands on a diagonal term. Rust had excluded such pops entirely, shrinking
  every multi-population f3/f4 toward zero on real data. Split `estimate_pop_stats` into a
  frequency mask + `het_valid`, returning `AAX_INVALID` for the het, and added the
  `yy < -99` drop in `eval_fstat`. (`admx-fstats/src/accumulator.rs`, `driver.rs`)
- **Anchor-presence gate (P0).** The per-SNP loop skipped SNPs where the anchor was
  absent, wrongly restricting the all-pairs FST accumulation and corrupting `lambdascale`
  and per-pop counts under `allsnps:YES`. Removed; `allsnps:NO` is still handled by the
  all-present check. (`admx-fstats/src/driver.rs`)
- **lambdascale `fstest` denominator (P0).** Must be the per-pair valid-SNP count
  (C `dofstnumx`), not the global SNP count. Added a `pair_bcnt` accumulator;
  `lambdascale` now matches C exactly (2.840 on the benchmark set).
- **qpfstats basis report columns (P1).** Printed the wrong four quantities; now
  `fbmean, sqrt(fbcovar diag) :: fsmean[t], fssig[t]` per C, including flattening the
  per-fstat mean (not just the sigma). Added `fs_mean`/`fs_sig` to the driver result.
- **CLI header `n`/`coefsize`.** Used `np*(np+1)/2`; now the true `full_stats` size.

### Fixed — qpWave / qpAdm
- **Eigenvector layout transpose (P0, ranktest).** Initial B was seeded from a row of
  the LAPACK Z matrix instead of the eigenvector column → wrong ALS fixed point and
  off chisq. (`admx-rank/src/ranktest.rs`)
- **qpWave `yvar` regularizer (P0).** Must be `addscaldiag` = `yscale * trace(yvar)`,
  not `yscale` absolute. (`admx-qpwave/src/driver.rs`)
- **qpAdm `yvar` regularizer (P0).** For qpAdm the default `diagvarplus` is 0.0 (not
  `yscale` as in qpWave), so the driver must add nothing to `yvar`.
  (`admx-qpadm/src/driver.rs`)
- **`solvitforcez` panic (P0).** Copied the full oversized `tdim²` work buffer into a
  `dim²` slice; now copies only the leading `dim×dim`. (`admx-rank/src/ranktest.rs`)

- **Genotype-path `inbreed` ignored (P1).** `QpWaveConfig`/`QpAdmConfig` and both CLIs now
  parse and plumb `inbreed:` through to the internal `run_qpfstats`; it was hardcoded
  `false`. With this plus the block-assignment fix, the genotype path now matches C.

Deterministic qpWave/qpAdm outputs match the C binaries on real 9-pop data via BOTH the
precomputed-`.fstats` path and the direct genotype path (`inbreed:YES`): qpWave chisq
51.110/7.640, qpAdm coefficients 0.757/0.243. Known remaining gaps: the qpAdm bootstrap RNG
(`prng.rs` placeholder LCG, non-deterministic std-errors) and the `hashets==0` variance
adjustment (does not fire under inbreed:YES).

## [2026-04-24] — Phase 1: Core F-Statistics & `qpfstats`

### Added
- **Optimized Accumulation Engine (`admx-fstats`):**
    - Implemented §3.1 matrix-based accumulation. Instead of evaluating $O(N_{stats})$ per SNP, we now accumulate a population-pair matrix $M_{ab}$ and a heterozygosity correction vector $H_a$.
    - This provides an algorithmic speedup from $O(L \cdot N^4)$ to $O(L \cdot N^2)$ for $L$ SNPs and $N$ populations.
- **Numerical Fidelity Core:**
    - **`DetailedPopCounts`:** Added support for the 5-element genotype count used in legacy `loadaa`, distinguishing between female (diploid) and male (haploid/X) genotypes.
    - **Inbreed Logic:** Faithful port of the complex `jhet` calculation and $aax$ (unbiased $p^2$) estimates from `qpsubs.c`.
    - **X-Chromosome Handling:** Integrated automated haploid logic for sex chromosomes and specified haploid markers.
- **I/O & Configuration:**
    - **`ParFile` Parser:** Faithful implementation of the `nicksrc/getpars.c` logic, including `keyword: value` parsing, inline comments, `+++` file includes, and uppercase string substitution (`dostrsub`).
    - **Format Auto-detection:** Added `is_packed_am` and `is_bed` helpers to distinguish between PACKEDANCESTRYMAP and PLINK BED formats.
    - **`.fstats` Writer:** Implemented standard and `hires` output formats for basis means and covariance matrices.
- **CLI Tools:**
    - **`qpfstats`:** End-to-end implementation of the f-statistics basis generator.

### Improved
- **Memory Efficiency:** Replaced massive per-statistic jackknife arrays with a compact $N \times N$ matrix-per-block approach.
- **Type Safety:** Transitioned from raw C pointers and 2D arrays to `ndarray` and strongly-typed Error variants.

### Fixed
Resolved entries previously tracked in `bugs.md`:

- **WLS consensus fit in `qpfstats` (P0).** Replaced the simplified per-pair basis accumulation in `admx-fstats/src/driver.rs` with the full C `dofstats` flow (`qpsubs.c:4665`): per-SNP evaluation of every canonical f-statistic (all f2 + all anchored f3 + all canonical f4), per-fstat block-jackknife for per-fstat sigmas, `gbot`+prior flattening, then a weighted least-squares compression onto the anchored-f3 basis via `wfb = fbcoeffs / jsig`, `wco = wfbᵀ·wfb` (regularized with `DOFSTATS_DIAG = 1e-8 · trace/nbasis + 1e-8`), Cholesky inverse, and a final basis-level `wjackvest`. Produces correct numbers under `allsnps: YES` with missingness, where the prior unweighted path silently returned per-pair-denominator estimates.
- **`pdinv` returned a corrupted inverse (P0, latent).** Added `admx_core::linalg::pdinv` as a `dpotrf → dpotri` pair. The initial triangle-mirror loop copied the uninitialized lower triangle back into the populated upper (the direction was reversed), wiping LAPACK's inverse with whatever junk `dgemm` had written into the lower half. Swapped the loop so populated upper values propagate into the lower triangle. Without this fix the WLS fit above produced basis means inflated by ~10³×.
- **`wjackvest` did not return the jackknife-corrected mean (P1).** The Rust signature only returned the covariance, forcing callers to use the raw full-sample estimate as the mean — whereas the C `wjackvest` writes back a corrected `jackest`. Extended the signature to return `(jackest, var)` and updated the sole caller (`admx-fstats/src/driver.rs`) to use the corrected vector as the output basis mean.
- **Empty jackknife blocks produced NaN (P0).** `wjackest` / `wjackvest` divide by per-block weight, and block assignment can leave holes (max block id > actual occupied blocks) when SNPs cluster. The driver now collects a single `active_blocks` list (`wjack[k] > 0`) and feeds only those into the per-fstat jackknife and the final basis-level jackknife.

- **Reader/SNP desync on ignored or out-of-block SNPs (P0).** `admx-fstats/src/driver.rs` now calls `read_record` at the top of the SNP loop, before the `snp.ignore` / `chrom` / block-index filters, so the `GenoReader` cursor and the `snps: Vec<Snp>` list stay in lockstep when `chrom:`, `nochrom:`, or `badsnpname:` drop SNPs.
- **Per-(individual, chromosome) haploidy (P0).** `Indiv` now carries a `Sex` enum (`admx-core/src/types.rs`). Haploidy is decided at SNP-read time based on chromosome (X/Y/MT) and sex, instead of a sticky per-individual flag. Males on autosomes are no longer treated as haploid.
- **`wjackvest` centering (P0).** Ported the Busing weighted jackknife from `qpsubs.c:wjackvest` — `jackest` accumulation and `xtau = h·mean − (h−1)·jmean[k] − jackest` — replacing the arithmetic-mean centering. `admx-core/src/jackknife.rs`.
- **Anchor population mismatch in `qpfstats` basis labeling (P0).** `admx-cli/src/bin/qpfstats.rs` reorders `pop_list` so `outpop:` sits at index 0 before the driver runs, matching the C reader's expectation and eliminating the means/covar-vs-labels drift when `outpop:` was not first.
- **Golden-log fixture (P1).** Added `tests/fixtures/test1` with sample `.geno` / `.ind` / `.snp` inputs and C reference logs (`test1.c.log`, `test1.c.fstats`, plus the `.no` / `.raw` / `.scaled` variants) alongside the Rust outputs for diffing.
- **Duplicate keys in `params.rs` (P2).** Downgraded from error to `log::warn!` to match C's "last one wins" semantics in `nicksrc/getpars.c`.
- **`dostrsub` macro regex (P2).** Broadened the upper-case key detector to accept digits and underscores, so macros like `DIR_1` expand correctly.
- **`+++` include-path resolution (P2).** Verified (and kept) CWD-relative resolution with no macro expansion on the include path itself — matches C.
- **`.gz` transparent decompression (P2).** `admx-io/src/lib.rs` adds a `Storage` abstraction that is either a `memmap2::Mmap` (uncompressed) or a `Vec<u8>` produced by `flate2::read::GzDecoder`. `is_packed_am` / `is_bed` sniffing works transparently through the compressed path.
- **Orphan modules (P2).** Removed `admx_io::meta` and the doc-only `geno.rs` / `packed.rs` stubs that were not declared in `lib.rs`. `admx_core::linalg` is retained — Phase 2 will need it.

---

## Future Work

### Immediate (Phase 1 Continued)
1. **`admx-rank`:** Implement SVD-based subspace projection and the `ranktest` family of statistics required for `qpWave`.
2. **`qpWave` & `qpAdm`:** Implement the linear model solvers and nested model testing.
3. **Golden-Log Validation:** Execute a comprehensive comparison suite against C AdmixTools outputs to ensure bit-level parity.

### Performance
4. **Rayon Integration:** Parallelize the SNP accumulation loop across jackknife blocks or SNP chunks.
5. **BLAS Optimization:** Replace manual accumulation loops in `accumulator.rs` with `dsyr` (Rank-1 Update) calls to high-performance BLAS backends.

### Robustness
6. **`load_fstats`:** Complete the parser for reading `.fstats` files back into the engine.
7. **Extended `dostrsub`:** Support recursive string substitution to match C's looping behavior for complex parameter nesting.
