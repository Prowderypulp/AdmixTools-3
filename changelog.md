# Changelog — AdmixTools Rust Port

## [2026-05-13] — `qpAdm` Performance: 13–18× over C reference

### Improved
- **`qpAdm` end-to-end wallclock vs C `qpAdm` (AdmixTools 8.0.1):**
  - `small` tier (20K SNPs × 70 inds × 7 pops): 0.091 s → 0.007 s — **13.4×**
  - `med` tier (100K SNPs × 200 inds × 10 pops): 0.264 s → 0.014 s — **18.3×**
- Numerical output is bit-for-bit identical to the previous build; the
  `test_qpadm_golden` harness and full `cargo test --workspace` pass unchanged.

### Added
- **`admx_core::linalg::set_blas_single_threaded`** (`admx-core/src/linalg.rs`).
  Calls `openblas_set_num_threads(1)` and `omp_set_num_threads(1)` so OpenBLAS
  and the underlying libgomp pool do not fan out across cores on the many small
  (≤16×16) BLAS calls inside the bootstrap. Profiling with samply showed
  **96.6% of qpAdm runtime was in libgomp barriers**, not in actual flops —
  the OMP fan-out cost dominated the per-call math.
- **`admx-core/build.rs`** — explicit `cargo:rustc-link-lib=dylib=gomp` so
  `omp_set_num_threads` resolves. OpenBLAS pulls libgomp in at runtime but
  does not expose its symbols on the link line.
- **Rayon-parallel bootstrap** in `admx-qpadm/src/bootstrap.rs::calcevarboot`.
  The 1000-iteration Gaussian-resample loop (`numboot` × 2 passes = 2000
  independent `doranktest`+`calcadm` calls) now runs on a dedicated rayon
  pool. Key safety properties preserved:
  - The Gaussian noise vector `rvec` is pre-drawn serially before the parallel
    section, so PRNG output is identical regardless of execution order.
  - `doranktest` itself is fully deterministic on this code path (no
    `gaussa_legacy`/`drand2` calls in `ranktest` — only `ranktestfix` uses
    them, and bootstrap does not invoke it).
  - Iteration `k=0` runs serially first so the `zzevarboot`/`zzevarboot2`
    diagnostic lines appear in the same order as the C reference.
  - The pool's `start_handler` calls `set_blas_single_threaded` on each
    worker. Without this, nested libgomp threading inside parallel BLAS calls
    **regressed wallclock 50× vs. serial** — calling the BLAS-config helper
    only on the main thread leaves rayon workers using the default OMP team
    size on their first BLAS call.

### Internal
- `admx-cli/src/bin/qpAdm.rs::main` calls `set_blas_single_threaded` at
  startup so the global single-thread state is set before any BLAS init.
- Added `admx-qpadm/plan.md` documenting profile data, the four-step plan,
  and the remaining (deferred) optimization headroom (scratch-allocation
  hoisting and inner-loop tightening in `admx-rank/src/ranktest.rs`).

### Why not AdmixTools 2-level (50–100×)?
The remaining gap to AT2-class numbers comes from AT2's different statistical
pipeline (closed-form SVD rank test + jackknife-only variance, no Gaussian-
resample bootstrap), not from implementation efficiency. Matching it would
break bit-parity with the C log surface that `test_qpadm_golden` checks.

---

## [2026-05-05] — Phase 2: `qpAdm` Golden Parity Completed

### Fixed
- **`qpAdm` bootstrap instability (P0):** Corrected rank initialization in `admx-rank/src/ranktest.rs` by fixing eigenvector extraction layout to match legacy C/LAPACK column-major behavior. This removed the large bootstrap outliers and restored stable `zzjmean`.
- **Bootstrap covariance scaling mismatch (P0):** Updated `admx-qpadm/src/bootstrap.rs` to use sample covariance denominator `n-1` (matching C `mv2D`) instead of `n`, aligning covariance magnitudes in golden output.
- **Fixed-pattern enumeration/order mismatch (P1):** Updated `admx-qpadm/src/driver.rs` to decode and print fixed patterns in C-compatible order (`00`, `01`, `10`, ...), resolving `fixed pat` table ordering differences.
- **Missing diagnostic output block (P1):** Implemented `worst Z-score with right hand mix` output in Rust (`admx-qpadm/src/driver.rs`) to match canonical qpAdm log structure.
- **Tail formatting mismatch in qpAdm summary path (P2):** Adjusted final summary/tail printing so normalized golden logs are stable and line-compatible with C references.

### Improved
- **Golden parity fidelity for qpAdm output surface:** Rust qpAdm now matches C reference output across bootstrap means/covariance and downstream fixed-pattern diagnostics used by the existing golden harness.

### Validation
- `cargo test -p admx-cli --test golden_log test_qpadm_golden -- --exact` ✅
- `cargo test -p admx-cli --test golden_log` ✅
- `cargo test --workspace` ✅

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
