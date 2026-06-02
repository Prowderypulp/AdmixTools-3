# AdmixTools → Rust, Phase 2

**Scope:** the rank engine (`f4rank.c`), `qpWave`, `qpAdm`, and the perf pass that
was §M8 of Phase 1. Explicitly out of scope (unchanged from Phase 1 §6):
`qpGraph`, `qpDstat`, `qp3Pop`, `qpF4ratio`, `qp4diff`, `qpBound`, `qpDpart`,
`qpff3base`, `qpfmv`, `qpfmvmix`, `qpmix`, `qpreroot`, `convertf`, `mergeit`,
`snpunion`, `multimerge`, `transpose`, `geno_single`, `merge_transpose`,
`rolloff`, `rolloffp`, `rexpfit`, `weightjackfourier`, `perlsrc/`.

**Fidelity contract:** unchanged from Phase 1. Same LAPACK calls on the same
matrices → bit-identical. Any algorithmic shortcut must be gated behind a
flag the legacy C binary already exposes (`newmode:`, `oldmode:`, etc.) so the
default path reproduces legacy output to `%12.6f` / `%15.9g`.

---

## 1. What landed in Phase 1

| Crate | State |
|---|---|
| `admx-core` | Types, jackknife (`wjackest`, `wjackvest`), block assignment, `linalg` with `dspev`/`dpotrf`/`dpotri`/`dsyr`/`dgemm`/`pdinv`. |
| `admx-io` | Packed ANCESTRYMAP + PLINK BED + EIGENSTRAT readers, `.snp`/`.ind` parsers, `getpars` (with `+++` includes and `dostrsub`), `.gz` transparent decompression. |
| `admx-fstats` | Full C `dofstats` WLS flow (`qpsubs.c:4665`). Canonical basis (anchored f3), per-fstat block-jackknife, WLS consensus onto basis, `.fstats` writer. |
| `admx-cli/qpfstats` | End-to-end. Parfile → `.fstats` file. |

The `.fstats` file is the handoff surface between Phase 1 and Phase 2:
qpWave and qpAdm both consume it as their primary input, and the writer is
done. The `loadfstats` reader (on the consumer side) is still a stub —
Phase 2 M2.5.

Four bugs worth remembering (see `admixtools-rs/changelog.md` for the full
list):

- Empty jackknife blocks produced NaN when block assignment left holes.
  `admx-fstats` now threads a single `active_blocks` filter through both
  the per-fstat and basis-level jackknife passes. **Replicate this pattern
  in qpAdm's `calcevar`.**
- `pdinv`'s triangle-mirror direction was initially reversed, silently
  wiping LAPACK's inverse. Tests did not catch it; a 3-pop fixture did.
  **Build a pdinv-on-random-SPD-matrices unit test in Phase 2 M2.0.**
- `wjackvest` originally only returned the covariance, so callers used
  the raw full-sample estimate as the output mean. It now returns
  `(jackest, var)`. **qpAdm bootstrap machinery must use `jackest`, not
  the raw mean.**
- The prior basis pipeline was a simplified unweighted version of C's
  WLS consensus. The current driver reproduces C's full flow. **Do not
  reintroduce the shortcut as a "perf win" without gating it behind a
  flag and a golden-log parity test.**

---

## 2. Phase 2 plan — milestones

Each milestone ends with a golden-log diff against the C reference on
fixtures grown out from `admixtools-rs/tests/fixtures/test1/` (3 pops, 17
SNPs, smoke-only) to at least one realistic case (~20 pops, ≥100k SNPs).
M2.0 builds the fixture; everything after gates on it.

### M2.0 — fidelity harness

Before writing qpWave/qpAdm math, stand up the test infrastructure they
need:

- Promote `tests/fixtures/test1` from a smoke to a real golden. Add a
  mid-size fixture (~20 pops, 100–200k SNPs — drop a slice of one of the
  datasets under `examples/`) with C reference logs for `qpfstats`,
  `qpWave`, and `qpAdm` checked in.
- `admx-cli/tests/golden_log.rs` harness: normalize timestamps, paths,
  and the `memory:` lines; diff `summ:`, `f4info:`, `best coefficients:`,
  `error covariance`, per-fstat rows at `%12.6f` / `%15.9g`.
- `cargo test --workspace` runs the harness in CI.
- **pdinv SPD regression test.** Inverse of a random SPD matrix, compare
  against `ndarray`/`nalgebra` dense inverse to `1e-10`. Non-negotiable —
  the latent pdinv mirror bug from Phase 1 would have been caught by
  this in five minutes.

### M2.1 — `admx-rank`

Port `src/f4rank.c` (~750 LOC). Two rank-test routines:

- **`ranktest`** — alternating least squares. The 6-deep scalar loop
  that builds the coefficient matrix has Kronecker structure
  (`Phase 1 §3.4`); replace with a fixed count of `dgemm` calls. Force
  reduction order to stay deterministic across BLAS thread counts.
- **`ranktestfix`** — same solver on a constrained system (rows of `A`
  zeroed per fix pattern). Reuse the GEMM kernel.
- **`normab`** — normalization + sign convention. Printed in the `A:` /
  `B:` blocks, so byte-for-byte matters.
- **`newmode: YES`** closed-form shortcut — Cholesky-whitened weighted
  SVD. Not default. Gated behind the legacy flag; separate golden
  fixture since output differs below the printed precision.

**Test:** `f4info:` blocks (rank, dof, chisq, tail, dofdiff, chisqdiff,
taildiff) diff against C on every candidate rank.

### M2.2 — `admx-qpwave`

Port `src/qpWave.c` (~1k LOC).

- **`doq4vecb`** — per-block f4 sums/ratios. `allsnps` short-circuit,
  `fancyf4` toggle, and `addscaldiag(yvar, diagvarplus, nl*nr)` call
  site preserved. `yscale = 1e-4` regularizer is a Phase 1 invariant;
  do not change.
- **`loadymv`** — the `fstatsname:` path. This is where the Phase 1
  stub `admx-fstats::fstats_io::load_fstats` becomes load-bearing;
  land that reader first (see M2.5 — may pull forward).
- **`checkmv`** — error handling with exact C exit codes.
- Rank-test driver loop from `0..n_left - 1`, chisq + tail
  probabilities printed in the legacy column widths.

**Test:** full `qpWave.log` diff.

### M2.3 — `admx-qpadm`

Port `src/qpAdm.c` (~2.5k LOC) — the biggest piece of Phase 2. Depends
on M2.1 and M2.2.

- **`calcadm` / `calcadmfix`** — weight solver with sum-to-one
  constraint. NNLS-style.
- **`setktable` + fix-pattern enumeration** — `2^nl` independent
  `ranktestfix` calls. Parallelize across fix patterns with rayon; do
  not parallelize inside the solver.
- **Nested-model p-values**, `summ:` line, `hires` output, `hiprec_covar`
  path.
- **`calcevar`** (delete-block jackknife) and **`calcevarboot`**
  (numboot = 1000 default). Reuse the active-blocks pattern from
  Phase 1 `admx-fstats` — never pass empty blocks into the jackknife.
- **PRNG.** `admx-qpadm/src/prng.rs` has ~43 lines of partial work.
  Finish the `ranmar` port; it must reproduce the C stream bit-for-bit
  so bootstrap replicate order matches (Phase 1 §7 Q2). The answer to
  that open question is **exact match** — bootstrap covariance is a
  printed field and downstream tooling parses it.
- **`doratio`, `gendstat`, `worstb`, `phiint`, `mcest`** — port as-is
  behind their existing parfile flags.
- **Kill `system("qpfstats ...")`.** The Phase 1 plan flagged this as
  the #1 layering violation (§3.6). Replace with a direct
  `admx_fstats::driver::run_qpfstats` call. Zero algorithmic change;
  removes process-startup + tempdir brittleness.

**Test:** full `qpAdm.log` diff — `details:` block, `best pat:` rows,
and the `error covariance` matrix under both `calcevar` and
`calcevarboot`.

### M2.4 — perf pass

This is Phase 1 M8 rolled forward. Gate on M2.0–M2.3 all green.

- Rayon over blocks in `dofstats` / `countpops` accumulation. Sort
  block-wise reductions by block id before final sum → deterministic
  across thread counts.
- Rayon over qpAdm fix patterns (already called out in M2.3).
- SIMD packed-genotype decode (AVX2 / NEON) behind a runtime
  `is_x86_feature_detected!` check. Scalar path stays as the fallback.
- Pin `OPENBLAS_NUM_THREADS=1` when parallelizing over SNPs — avoids
  oversubscription and locks in reduction order.
- Benchmark matrix: `(numeg ∈ {20, 50, 100}) × (ncols ∈ {100k, 600k, 1.2M})`
  against current C binaries.
- **Target:** ≥10× on the `allsnps` pipeline (dominated by §3.1 and
  §3.4), ≥3× on legacy `allsnps: NO`. Golden-log gate stays green
  throughout — any perf change that moves output by more than a ULP
  needs a flag.

### M2.5 — robustness tail

Cleanup items that are not blocking but need to land before the first
external release:

- **`load_fstats`** reader. Currently errors "not yet fully implemented".
  Pull forward to before M2.2 if `loadymv` needs it.
- **Recursive `dostrsub`** — C applies substitution in a loop until
  fixed-point. Phase 1 ships iterative-until-convergence, but we've
  only tested flat cases. Confirm matches on nested macros.
- **`log` subscriber.** Route at WARN to stderr by default. **Never let
  `log::warn!` leak to stdout** — downstream Perl in `perlsrc/` parses
  stdout and extra lines will break it.
- **CI.** `cargo check`, `cargo test --workspace`, golden-log harness,
  and `cargo clippy -- -D warnings` on a minimal lint set.

---

## 3. Invariants (unchanged, do not touch)

Same as Phase 1 §5 — restated for completeness:

- All printed field widths and precisions exactly as the C code emits
  them. Logs are machine-parsed.
- Regularizer constants: `yscale = 0.0001` (qpWave), qpAdm's
  `yscale / diagvarplus` fallback chain, `addscaldiag` call sites,
  `diag = 1e-8` in `dofstats`.
- `blgsize = 0.05` Morgans default; `tagnumber` assignment order.
- `oldmode = YES` default in qpAdm, `fancyf4 = YES` default,
  `numboot = 1000` default.
- Exit codes + `fatalx` strings when observable through exit status.
- `fstats` file format, including the `basepop:` header line parsed
  by `fstats2popl` on readback.

---

## 4. Open questions

1. **Bootstrap PRNG (revived from Phase 1 §7 Q2).** Commit to bit-exact
   `ranmar`. Anything else means `error covariance` output drifts and
   we can never close the golden-log gate.
2. **M2.1 `newmode: YES` closed-form path.** Ship in Phase 2 or defer?
   Default path is already iterative-LS + GEMM, which covers the
   fidelity gate. The SVD shortcut is a pure perf win and has its own
   golden-log surface. **Recommendation:** ship in Phase 2 but gate
   behind `newmode: YES` and a separate fixture.
3. **Structured JSON side-channel.** Phase 1 §7 Q4 punted. Some
   downstream users are starting to parse logs with brittle regex.
   Land a `--json-stdout` side-channel in M2.5 — zero cost if
   nobody sets it, and future-proofs against log format churn.
4. **Test fixture size.** `tests/fixtures/test1` is 3 pops / 17 SNPs.
   Fine for smoke, not for measuring anything. Pick one of the
   `examples/` datasets, subset to ~20 pops / 100k SNPs, commit raw
   genotypes + C reference logs. Decide in M2.0.

---

## 5. Dependencies and sequencing

```
M2.0 (harness + pdinv regression)
  └── M2.1 (admx-rank)
        └── M2.2 (admx-qpwave)              ← also depends on M2.5 load_fstats
              └── M2.3 (admx-qpadm)
                    └── M2.4 (perf pass)
                          └── M2.5 (robustness tail, CI)
```

`load_fstats` (M2.5) is a hard dep for M2.2 `loadymv`; pull forward
if needed — it's ~200 LOC of parsing, not a blocker.

M2.4 can start opportunistically against M2.1 once the rank engine
is bit-parity with C; no need to wait for M2.3 to finish.
