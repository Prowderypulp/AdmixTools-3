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
  and nested-model p-value match C (e.g. 0.757 / 0.243 on the benchmark split). Bootstrap
  std-errors not yet bit-identical (placeholder LCG in `prng.rs` — P2).
- `admx-cli/qpWave`, `admx-cli/qpAdm` — full CLIs (parfile → run → C-format output).

`qpfstats` matches C to ~6 decimals (means **and** sigmas) under `inbreed:YES` /
`allsnps:YES`, and runs ~17.5× faster than C on the benchmark dataset.

### Golden tests
`admx-cli/tests/golden_log.rs`: `test_qpfstats_golden`, `test_qpwave_golden`, and
`test_qpadm_golden` (now a real test) all pass.

### Known remaining gaps (P2)
- qpAdm bootstrap RNG (`admx-qpadm/src/prng.rs`) — non-deterministic std-errors.
- `hashets==0` variance adjustment unported (does not fire under `inbreed:YES`).

## Out of scope (both phases)

`qpGraph`, `qpDstat`, `qp3Pop`, `qpF4ratio`, `qp4diff`, `qpBound`, `qpDpart`,
`qpff3base`, `qpfmv`, `qpfmvmix`, `qpmix`, `qpreroot`, `convertf`, `mergeit`,
`snpunion`, `multimerge`, `transpose`, `geno_single`, `merge_transpose`,
`rolloff`, `rolloffp`, `rexpfit`, `weightjackfourier`, `perlsrc/`.
