# AdmixTools-rs — Progress

Plan labels map to milestones in `../RUST_PORT_PHASE1.md` and `../RUST_PORT_PHASE2.md`.

Status legend: **done** · **partial** · **stub** · **todo**

---

## Phase 1 — closed (2026-04-24)

Scope delivered: `admx-core`, `admx-io`, `admx-fstats`, and the `qpfstats` CLI.
`admx-rank`, `admx-qpwave`, `admx-qpadm` — originally listed in Phase 1 §1 —
roll to Phase 2. The fstats file format is the handoff surface, and the
producer side is now end-to-end.

### `admx-core` — **done**
- `types.rs` — `Indiv` (with `Sex` enum), `Snp`, `PopCounts`,
  `DetailedPopCounts`, `GenoFormat`, `Block`, `FStatKind`, `F4Info`.
- `constants.rs` — `DEFAULT_BLGSIZE = 0.05`, `DOFSTATS_DIAG = 1e-8`,
  `QPWAVE_YSCALE = 1e-4`, `DEFAULT_NUMBOOT = 1000`, `FSTAT_MISSING = -9999.0`.
- `error.rs` — `AdmxError` / `AdmxResult`.
- `jackknife.rs` — `wjackest` (scalar) and `wjackvest` (vector, returns
  `(jackest, var)` to match the C signature).
- `blocks.rs` — `assign_blocks` over Morgan positions.
- `linalg.rs` — wrappers around `dspev`, `dpotrf`, `dpotri`, `dsyr`,
  `dgemm`, plus `pdinv` (`dpotrf → dpotri → symmetric mirror`).

### `admx-io` — **done** (for Phase 1 scope)
- `codec.rs` — canonical 2-bit MSB-first encoding.
- `packed_am.rs` — mmap'd PACKEDANCESTRYMAP with header hash parsing.
- `packed_ped.rs` — PLINK `.bed` (SNP-major; individual-major rejected).
- `eigenstrat.rs` — text genotypes.
- `snp.rs`, `indiv.rs` — `.snp` / `.ind` metadata parsers.
- `params.rs` — `getpars` port: `keyword: value`, inline comments,
  `+++` includes, `dostrsub` uppercase substitution (broadened regex),
  duplicate keys → warn (C's "last one wins").
- `lib.rs` — `GenoReader`, `Layout`, `Storage` (mmap or `flate2`-decoded
  buffer), `is_packed_am` / `is_bed` / `is_eigenstrat` sniffers.

Remaining (Phase 2 consumer side): `load_fstats` round-trip reader,
recursive `dostrsub`, stderr-routed `log` subscriber.

### `admx-fstats` — **done**
- `basis.rs` — anchored-f3 basis + canonical full-stats enumeration +
  sparse `coefficients()` map from any f2/f3/f4 into the basis.
- `accumulator.rs` — `DetailedPopCounts` → `(p, aax, valid)` via
  `estimate_pop_stats` with inbreed path and haploid marker handling.
- `driver.rs` — full C `dofstats` flow: per-SNP evaluation of every
  canonical f-statistic, per-fstat block-jackknife, `gbot`+prior
  flattening, WLS consensus fit onto the basis (regularized by
  `DOFSTATS_DIAG`), final basis-level `wjackvest`.
- `fstats_io.rs` — `dump_fstats` (hires and non-hires, basepop header
  line matching C). `load_fstats` is stubbed (Phase 2).

### `admx-cli/qpfstats` — **done**
End-to-end CLI. Parfile → metadata → reader → driver → `.fstats` writer.
Reorders `pop_list` so `outpop:` is index 0 before the driver runs.

### Cross-cutting
- Fixture: `tests/fixtures/test1` with C reference logs/fstats plus
  Rust outputs. Small (3 pops, 17 SNPs); good enough for a smoke test,
  not for measuring perf.
- `cargo check`, `cargo build`, `cargo test --workspace` all green.
- No CI yet (Phase 2).

---

## Phase 2 — implemented & validated against the C binaries (2026-06-02)

The stubs below are all implemented. Validated against `/home/drtex/AdmixTools/bin`
on the small fixtures and a real 1.23M-SNP × 893-indiv dataset. See `changelog.md`
(2026-06-02), `bugs.md`, and `bench/RESULTS.md` for the fidelity pass and bug fixes.

- `admx-rank` — **done**. `ranktest`, `ranktestfix` (ALS), `normab`, `checkmv`, `solvit`.
- `admx-qpwave` — **done**. Matches C exactly via precomputed `.fstats` AND the direct
  genotype path (`inbreed:YES`): chisq 51.110 / 7.640 on the benchmark split.
- `admx-qpadm` — **done** (deterministic). Coefficients, fixed-pattern chi-square table,
  and nested-model p-value match C (e.g. 0.757 / 0.243 on the benchmark split). Now also
  validated under `inbreed:NO` on real data (maitrus/seq): chisq `8.544` and coefficients
  `0.311 0.492 0.125 0.071` are an **exact** match to C after porting the `hashets==0`
  sigma adjustment. Bootstrap noise generator is now **bit-identical** to C (Marsaglia
  polar `gauss()` + Numerical-Recipes `choldc`, verified against `libnick.a`); with a
  fixed `seed:` the std-errors match C at printed precision (`0.088 0.127 0.026 0.053`).
  Residual covariance difference (~0.05%) is downstream linalg, not the RNG — P2.
- `admx-cli/qpWave`, `admx-cli/qpAdm` — full CLIs (parfile → run → C-format output).

`qpfstats` matches C to ~6 decimals (means **and** sigmas) under `inbreed:YES` /
`allsnps:YES`, and runs ~17.5× faster than C on the benchmark dataset. On the real
maitrus/seq run it reports `adjusted sigs: 1440` / `lambdascale: 3.160`, both matching C.

### Performance (real-data, maitrus/seq, measured `/usr/bin/time`)
qpAdm: wall **52 s** vs C **211 s** (~4× faster, port parallelizes; C single-threaded),
~same total CPU (215 s vs 209 s), peak RSS **2.29 GB** vs **3.06 GB** (~25% less). See
`real_benchmark.md`. (Ignore the C `##end of qpAdm` footer — it undercounts ~25×.)

**Parallel f-stat scan (2026-06-02):** the dominant per-SNP scan was parallelized over
jackknife blocks (output bit-identical — per-block sums keep serial SNP order): SNP scan
**35 s → ~8 s** (~4.3×), qpfstats end-to-end wall ~10 s. `lto=fat` build profile added.
RSS unchanged (mmap-dominated). The qpAdm bootstrap stays serial (OpenBLAS isn't safe to
enter concurrently — a parallel version deadlocked). See `changelog.md` / `real_benchmark.md`.

### Golden tests
`admx-cli/tests/golden_log.rs`: `test_qpfstats_golden`, `test_qpwave_golden`, and
`test_qpadm_golden` (now a real test) all pass.

### Known remaining gaps (P2)
- qpAdm bootstrap covariance residual — noise generator is now bit-identical (given
  identical inputs), so the RNG is ruled out; a ~0.05% covariance residual remains,
  most likely from the upstream `yvar` (already ~1e-6 off C) propagating through, plus
  downstream linalg rounding (LAPACK vs nicksrc). std-errors match C at printed precision.
- ~~`hashets==0` variance adjustment unported~~ — **ported 2026-06-02**
  (`admx-fstats/src/driver.rs`); validated under `inbreed:NO` on real data.

## Out of scope (both phases)

`qpGraph`, `qpDstat`, `qp3Pop`, `qpF4ratio`, `qp4diff`, `qpBound`, `qpDpart`,
`qpff3base`, `qpfmv`, `qpfmvmix`, `qpmix`, `qpreroot`, `convertf`, `mergeit`,
`snpunion`, `multimerge`, `transpose`, `geno_single`, `merge_transpose`,
`rolloff`, `rolloffp`, `rexpfit`, `weightjackfourier`, `perlsrc/`.
