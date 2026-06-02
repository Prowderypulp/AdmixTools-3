# Real Benchmark

Date: 2026-06-02

Workspace used for the real-data check:
`/home/drtex/Disk2/QP/Work/maitrus/seq`

Input configuration used:
- `parfile.txt`
- direct genotype path:
  - `genotypename: final.geno`
  - `snpname: final.snp`
  - `indivname: final.ind`
  - `popleft: left1.txt`
  - `popright: right1.txt`
  - `details: YES`
  - `allsnps: YES`
  - `summary: YES`
  - `inbreed: NO`

Original DReichLab output checked:
- `/home/drtex/Disk2/QP/Work/maitrus/seq/out.log`
- footer in that file:
  - `##end of qpAdm:        8.421 seconds cpu        2.953 Mbytes in use`

Rust output checked:
- direct run of:
  - `/home/drtex/Projects/admixtools-rs/target/release/qpAdm -p parfile.txt`

Observed comparison on the same config (BEFORE fix):
- Original qpAdm base model:
  - `chisq: 8.544`
  - `tail: 0.287109382`
  - coefficients: `0.311 0.492 0.125 0.071`
- Rust qpAdm base model:
  - `chisq: 7.692`
  - `tail: 0.360528`
  - coefficients: `0.255 0.580 0.125 0.040`

## Resolution (2026-06-02)

Root cause: the C `hashets == 0` per-fstat sigma adjustment
(`qpsubs.c:4864-4882`) was not ported. For each f-statistic, if a non-inbred
population with **no observed heterozygous genotype call** (C `counthets`,
i.e. any individual genotype `g == 1` — distinct from the model-based hetrate)
appears on a diagonal term (its index occurs >1× among the stat's four
`(a,b,c,d)` indices), C inflates `jsig = sqrt(jsig^2 + 100)`, heavily
downweighting that fstat in the WLS consensus fit. Pseudo-haploid aDNA pops
(genotypes only `0`/`2`) have a positive hetrate but zero het *calls*, so they
trigger it. The adjustment is gated on the inbreed flag (under `inbreed:YES`
every pop is inbred and C skips it) — which is exactly why the original
`inbreed:YES` fidelity pass never exercised this path.

Fix: `admx-fstats/src/driver.rs` now tracks `het_seen` per pop in the genotype
loop and applies the `sqrt(jsig^2 + 100)` inflation between the per-fstat
jackknife and the flatten step. qpfstats now prints the real `adjusted sigs`
count. qpWave's direct path shares the same `run_qpfstats`, so it is covered by
the identical fix (not separately ported). The basis-anchor pop ordering in both
qpAdm and qpWave was also aligned to C's `mkfstats` (right pops first), though
that is numerically a no-op (the rank/f4 pipeline is anchor-invariant).

After fix — Rust matches C on this run:
- qpAdm base model: `chisq: 8.544`, coefficients `0.311 0.492 0.125 0.071` —
  exact match.
- qpfstats `adjusted sigs: 1440` and `lambdascale: 3.160` both match C.
- qpWave (same parfile) matches C across the full rank ladder, e.g.
  codimension-1 `f4rank: 3 chisq: 8.386 tail: 0.2998` (identical); remaining
  differences are last-digit rounding (e.g. rank-0 chisq 8699.256 vs 8699.263).
- Caveat: the bootstrap std-errors differ slightly (e.g. 0.094 vs C's ~0.088)
  because the legacy LCG bootstrap PRNG stream is not bit-identical. Headline
  statistics (chisq, tail, coefficients) match.

Conclusion:
- The Rust port now matches the original on this real-data direct-genotype run
  under `inbreed:NO`.

## Time & memory (2026-06-02)

Measured with `/usr/bin/time -v` on both binaries (same parfile, same machine).

NOTE: the C `##end of qpAdm:  8.421 seconds cpu  2.953 Mbytes in use` footer is
the in-binary allocator counter, NOT real usage — it undercounts wall time by
~25× and peak RSS by ~1000× (it ignores the 2 GB genotype file and libc
allocations). The real measured C run is recorded in `c.time` and is what the
table below compares against.

| Metric        | C binary (`c.time`) | Rust (release) | Verdict     |
|---------------|---------------------|----------------|-------------|
| Wall (real)   | 210.71 s            | 52.06 s        | ~4× faster  |
| CPU (user)    | 209.41 s            | 215.32 s       | ~same (+3%) |
| Sys           | 0.71 s              | 0.62 s         | ~same       |
| Peak RSS      | 3.06 GB (3058348 KB)| 2.29 GB (2402968 KB) | ~25% less |

Reading it:
- Total CPU work is essentially identical (215 s vs 209 s) — the algorithms do
  the same amount of computation, as expected from a faithful port.
- Wall-clock is ~4× faster because the Rust port parallelizes (215 s user /
  52 s wall ≈ 4 threads busy); the C binary is single-threaded (209 s user ≈
  210 s wall).
- Peak memory is ~25% lower. Both runs are dominated by the 2.08 GB `final.geno`.

Caveat for future benchmarking: always use external timing (`/usr/bin/time -v`),
never the C footer.

## Optimization: parallel f-stat scan (2026-06-02)

Profiling (`ADMX_PROFILE=1`) attributed ~35 s of the ~50 s qpAdm run to the
single-threaded per-SNP f-statistic scan. Parallelizing it over jackknife blocks
(disjoint contiguous SNP ranges → per-block sums keep serial SNP order → output
bit-identical) plus `lto=fat`:

| Phase                     | Before  | After   | Note                         |
|---------------------------|---------|---------|------------------------------|
| SNP f-stat scan           | 35.1 s  | ~8.1 s  | ~4.3×, rayon over blocks     |
| qpfstats end-to-end wall  | ~37 s   | 10.4 s  | bootstrap-free; ~12× CPU util|
| Peak RSS                  | 2.29 GB | 2.29 GB | mmap-dominated, unchanged    |

Output verified byte-identical: full qpAdm and qpfstats (incl. hetrate) diffed
against pre-optimization baselines, and the `golden_log` C-fidelity tests pass.

The qpAdm `numboot` bootstrap (~7–23 s, noisy) was left serial: its `doranktest`
calls enter OpenBLAS, and the system OpenBLAS build deadlocks under concurrent
entry from rayon workers. It is unchanged from before. See `changelog.md`.
