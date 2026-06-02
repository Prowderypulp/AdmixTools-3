# Changelog — AdmixTools Rust Port

## [2026-06-02] — Performance: parallel f-stat scan (output bit-identical)

Made the dominant computation ~4× faster on real data with **byte-for-byte
identical output** (verified three ways — see Verification below).

### The hot path
Profiling the maitrus/seq run (`ADMX_PROFILE=1`, env-gated stderr timing added to
`run_qpfstats`/`run_qpadm`) showed the per-SNP f-stat scan dominates: **~35 s** of
a ~50 s run; jackknife/WLS were <1 s. The scan was single-threaded.

### Parallel SNP scan (`admx-fstats/src/driver.rs`)
- The jackknife block accumulators (`btop`/`bbot`/`wjack`/`pair_*`) are already
  *per-block*, and blocks are disjoint contiguous SNP ranges. So parallelizing
  **over whole blocks** keeps every accumulator summed in the exact serial SNP
  order — bit-identical, no cross-block float regrouping. Implementation:
  (1) parallel polymorphism flag per SNP; (2) serial block-boundary assignment
  (metadata only, replicating the serial state machine); (3) rayon over blocks;
  (4) serial merge in block order. The hetrate globals (`sum_p`/`n_p`, qpfstats
  only) are merged per-block; the reorder stays below the printed precision
  (verified — qpfstats output incl. hetrate is byte-identical).
- Needs random SNP access: added `GenoReader::random_access()` returning the
  mmap byte slice + record layout (`admx-io`), implemented for PACKEDANCESTRYMAP.
  Streaming/text readers (EIGENSTRAT) return `None` and use the unchanged serial
  loop. The whole new path is gated on this — the serial loop is preserved.
- Micro-opt: precompute a compact `(ind, pop)` list so the inner loop visits only
  the ~118 in-poplist individuals instead of all 893 per SNP (output-preserving).
- Per-SNP scratch (`counts`/`p`/`aax`/masks) is now reused per thread instead of
  reallocated for every SNP.

### Build profile (`Cargo.toml`)
- Added `[profile.release]` `lto = "fat"`, `codegen-units = 1` — whole-program
  inlining across the io/fstats boundary. Provably output-neutral; ~3% on its own.

### Verification (bit-identical invariant, now permanent)
1. Full qpAdm and qpfstats outputs (incl. hetrate) diffed **byte-identical**
   against pre-optimization baselines on the 1.2M-SNP real dataset.
2. New `cargo test` equivalence tests in `admx-fstats/src/driver.rs`
   (`parallel_matches_serial_allsnps_yes` / `_no`): a multi-block in-memory packed
   fixture is run through the parallel path and through a `ForceSerial` wrapper
   (`random_access()` → `None`), asserting `means`/`covar`/`w3`/`lambdascale`/
   `fst`/`f2`/`fs_*` bit-identical via `f64::to_bits`, with hetrate equal within
   fp-reorder tolerance. Covers **both** `allsnps` modes — the `!all_present`
   branch the real-data run (allsnps:YES) never exercised.
3. The 3 `golden_log` C-fidelity tests and the 2 RNG bit-identity tests stay green.

### Results (maitrus/seq, 16-core, page cache warm)
- SNP scan: **35.1 s → ~8.1 s** (~4.3×). qpfstats end-to-end wall ~10.4 s
  (User CPU ~122 s ≈ 12× parallel utilization).
- Peak RSS unchanged (~2.29 GB): it is the demand-paged mmap of the 2.08 GB
  `final.geno`, not heap — there is little to cut without changing the IO model.

### Not done — bootstrap left serial (deliberate)
- The qpAdm `numboot` bootstrap (2000 independent `doranktest` calls) is the other
  big phase, but each call enters OpenBLAS (`dspev`/`dgemm`/`pdinv`). The system
  OpenBLAS build is **not safe to enter concurrently** from multiple application
  threads — a rayon-parallel version (single-threaded BLAS per worker) deadlocked.
  Reverted; the bootstrap stays serial with multi-threaded BLAS per call. Its
  wall time is BLAS-thread/-load dependent and noisy (~7–23 s observed); it is
  unchanged from before this work. A safe parallelization would require pure-Rust
  replacements for the small-matrix LAPACK calls (future work).

## [2026-06-02] — qpAdm bootstrap RNG made bit-identical to C

The `numboot:` bootstrap std-errors now match the C binary at printed precision
under a fixed `seed:`. Two defects in the noise generator, both fixed:

### Fixed — bootstrap noise generator (`admx-qpadm/src/{prng,bootstrap}.rs`)
- **Wrong gaussian method.** Rust used Box-Muller (sin/cos, 2 uniforms/pair, no
  caching); C `gauss()` (`nicksrc/gauss.c`) is the Marsaglia **polar** method —
  rejection-sampled pair on the unit disc with a static `iset`/`gset` cache that
  returns one value and stores the other. Different algorithm → different sequence
  from identical uniforms. Ported exactly; the cache is now instance state on
  `LegacyLcg` and must persist across all `gauss()` calls in a run.
- **Wrong Cholesky.** `genmultgauss` factored the covariance with LAPACK `dpotrf`;
  C uses the Numerical-Recipes `choldc` (`nicksrc/linsubs.c`) with a **descending**
  inner accumulation (`k = i-1 .. 0`), which rounds differently. Replaced with a
  `choldc_lower` port. The existing multiply scaffold was already correct against
  C's transposed-Cholesky `mulmat` (ascending `k`, terms beyond `j` are exactly 0).

### Verified — bit-for-bit against the legacy library
- New tests `gauss_matches_legacy_c` / `genmultgauss_matches_legacy_c` assert
  `f64::to_bits` equality against values generated by C's own `libnick.a`
  (`SRAND(1923698036)` → `gaussa`, then `genmultgauss` on a known 3×3 covariance).
- Real-data qpAdm (maitrus/seq, `seed: 1923698036`): std-errors
  `0.088 0.127 0.026 0.053` — exact match to C at printed precision (was `0.094 …`).
- The global glibc `random()` state is process-wide, so a `RNG_TEST_LOCK` mutex now
  serializes RNG-touching tests (Rust runs tests in parallel threads).

### Known limitation (RNG ruled out, not yet localized)
- The raw error covariance still differs from C at the 4th significant figure
  (~0.05%, e.g. `7782` vs `7778` ×1e6). The generator is proven bit-identical *given
  identical inputs*, so the RNG is not the cause. The likely dominant driver is the
  `yvar` covariance fed into `genmultgauss`: the upstream f-stat/jackknife pipeline
  already differs from C at ~1e-6 (see `real_benchmark.md`, rank-0 chisq 8699.256 vs
  8699.263), and that propagates through the correct generator. Downstream rank-test /
  ALS / eigensolver / `calcadm` rounding (LAPACK vs nicksrc) may also contribute. Not
  localized further — out of scope for the RNG work. See `bugs.md`.

### Method note
- The original `0.088` could not have been reproduced as-is: the parfile had no
  `seed:`, so C used a time/pid `seednum()`. C *prints* the seed it chose
  (`seed: 1923698036`); set that explicitly in both binaries to get a fixed target.

## [2026-06-02] — `inbreed:NO` fidelity + performance benchmark

Closes the last known real-data divergence. On the maitrus/seq dataset
(`final.geno`, 5 left × 11 right pops, `inbreed:NO allsnps:YES`) qpAdm now
matches the C binary exactly on the headline statistics, and a measured
time/memory benchmark shows the port is competitive-to-better than C.

### Fixed — qpfstats `hashets == 0` per-fstat sigma adjustment (P0, `inbreed:NO` only)
- C `qpsubs.c:4864-4882`: for each f-statistic, if a **non-inbred** population with
  **no observed heterozygous genotype call** (C `counthets`, i.e. any individual
  genotype `g == 1` — distinct from the model-based hetrate) appears on a diagonal
  term (its index occurs >1× among the stat's four `(a,b,c,d)` indices), C inflates
  `jsig = sqrt(jsig^2 + 100)`, heavily downweighting that fstat in the WLS consensus
  fit. Pseudo-haploid aDNA pops (genotypes only `0`/`2`) have a positive hetrate but
  zero het *calls*, so they trigger it. The block is gated on the inbreed flag (under
  `inbreed:YES` every pop is inbred and C `continue`s for all of them), which is
  exactly why the earlier `inbreed:YES` fidelity pass never exercised this path.
- Fix: `admx-fstats/src/driver.rs` now tracks `het_seen` per pop in the genotype loop
  and applies the `sqrt(jsig^2 + 100)` inflation between the per-fstat jackknife and the
  flatten step, after the `gbot[j] < .001` zeroing guard (so no-valid-SNP stats are
  skipped and don't count toward `numadj`). `QpfstatsResult` gains `num_adjusted`;
  `qpfstats` now prints the real `adjusted sigs` count instead of a hardcoded `0`.
- The basis-anchor pop ordering in qpAdm and qpWave drivers was aligned to C's
  `mkfstats` (RIGHT pops first; with no basepop, `poplist[0]` = first right pop). For
  qpWave this is numerically a no-op (anchor-invariant rank test); for qpAdm it changes
  the regularized WLS anchor and is part of the chisq/coefficient correction below.

### Verified — real-data match (maitrus/seq, `inbreed:NO`)
- qpAdm base model: `chisq: 8.544`, coefficients `0.311 0.492 0.125 0.071` — **exact**
  match to C (was `7.692` / `0.255 0.580 0.125 0.040` before the fix).
- qpfstats: `adjusted sigs: 1440` and `lambdascale: 3.160` — both match C.
- All golden fixtures and the `inbreed:YES` path remain unchanged (16/16 tests pass).
- Caveat: bootstrap std-errors still differ at the ~3rd decimal (legacy LCG PRNG stream
  is not bit-identical); headline statistics match.

### Benchmark — time & memory (maitrus/seq, measured via `/usr/bin/time`)
Do **not** use the C `##end of qpAdm: 8.421 seconds / 2.953 Mbytes` footer for
benchmarking — that is the in-binary allocator counter and undercounts wall/RSS by
~25×/~1000×. Against the real measured C run (`c.time`):

| Metric        | C binary | Rust (release) | Verdict          |
|---------------|----------|----------------|------------------|
| Wall (real)   | 210.71 s | **52.06 s**    | ~4× faster       |
| CPU (user)    | 209.41 s | 215.32 s       | ~same (+3%)      |
| Peak RSS      | 3.06 GB  | **2.29 GB**    | ~25% less        |

Same total compute (faithful port), ~4× faster wall via parallelism (C is
single-threaded), and a leaner footprint. Both runs are dominated by the 2.08 GB
`final.geno`.

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
