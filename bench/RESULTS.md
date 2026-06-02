# Benchmark: Rust vs C AdmixTools

Date: 2026-06-02. Machine: local (Linux 7.0.8-zen). Warm page cache for both.

## Dataset
`/home/drtex/Disk2/QP/fstatistic/final.{geno,snp,ind}` — 1,233,013 SNPs × 893 individuals (~1.1 GB packed).
9 populations (see `bench/poplist`). Params: `allsnps: YES`, `hires: YES`, `inbreed: YES`, `noxdata: YES`.

## qpfstats — timing (this is where the compute is; reads all SNPs)

Measured after all fidelity fixes below, so both binaries process the same SNP set
and produce the same numbers (fair comparison).

| Impl | User CPU | Wall | Max RSS |
|------|---------:|-----:|--------:|
| C `qpfstats`    | 78.1 s | 1:19 | 1.05 GB |
| Rust `qpfstats` |  4.46 s | 0:04 | 1.36 GB |

**Rust ≈ 17.5× faster on CPU (≈ 21× wall)**, at ~30% higher peak RSS.
Source of the win: matrix-based O(L·N²) accumulation vs C's per-statistic loop.

## qpfstats — fidelity (matches C)

After the fixes, on this inbreed:YES / allsnps:YES dataset:
- header `np: 9 n: 666 nbasis: 36 coefsize: 23976` — exact
- per-pop `hetrate` / `valid snps` — exact (e.g. Uzbekistan 0.182 / 776583)
- `lambdascale: 2.840` — exact
- all 36 basis statistics: max abs diff vs C = **7e-6** (fitted mean), **4e-6** (raw).
  Residual is float accumulation-order noise, not a bug.

Plus the small golden fixture (inbreed:NO) still passes bit-for-bit.

## qpWave / qpAdm — validated via the PRECOMPUTED-FSTATS path (the standard workflow)

Both Rust binaries read `bench/c.fstats` (9 real pops) and were diffed against the C
binaries reading the same file, on a 3-left / 4-right split:

- qpWave: rank-0 chisq **50.244**, rank-1 **8.772**, tails to all printed digits — exact.
- qpAdm: coefficients **0.752 / 0.248**, chisq **8.787**, all three fixed-pattern rows
  (00/01/10) and the nested-model p-value (**0.009228**) — exact.

Bootstrap covariance / std-errors are RNG-dependent and not yet bit-identical (placeholder
LCG in `admx-qpadm/src/prng.rs`, plus a differing seed).

### Genotype path (direct from .geno, inbreed:YES) — now MATCHES C
After the block-assignment fix (below) and `inbreed` plumbing:
- qpWave genotype path: chisq **51.110 / 7.640** — matches C (tails agree to ~6 sig figs;
  deep-tail digits are float noise).
- qpAdm genotype path: coefficients **0.757 / 0.243** — matches C.

(C's genotype path is just qpfstats over the split's pops followed by the same solver —
verified: "C qpfstats on 7 pops → precomputed qpWave" reproduces C's genotype result
exactly. There is no separate `doq4vecb` discrepancy.)

## Bugs found & fixed during this comparison
All in the f-stat engine, each verified against the C binary:
0. **Jackknife block assignment over the wrong SNP set** (`admx-fstats/src/driver.rs`).
   Blocks were assigned over ALL SNPs up front; C assigns over only the used
   (non-ignored + polymorphic) SNPs (`setblocks` on the post-`rmsnps` list). Monomorphic
   SNPs shifted block anchors, perturbing per-block membership → the jackknife *variance*
   (covariance) diverged ~1-2% while means matched. Fixed by computing the block index
   in-loop over used SNPs. This closed the last residual: qpfstats now matches C to 6
   decimals (means AND sigmas) and the genotype-path qpWave/qpAdm match C.
1. **anchor-presence gate** (`admx-fstats/src/driver.rs`) wrongly restricted the
   all-pairs FST accumulation → corrupted lambdascale & per-pop counts. Removed.
2. **lambdascale `fstest` denominator** must be the per-pair valid-SNP count
   (C dofstnumx), not the global SNP count. Added `pair_bcnt` accumulator.
3. **frequency- vs heterozygosity-validity** (`accumulator.rs` + `driver.rs`): a
   low-sample inbred pop has a valid frequency but no aax; C keeps it usable in any
   f-stat where it is off-diagonal (drops the stat only when the -999 aax hits a
   diagonal term, via `yy < -99`). Rust had been excluding the pop entirely, which
   shrank every multi-pop f3/f4 toward zero. Split into `mask` (freq) + `het_valid`
   (aax) and added the `yy < -99` drop in `eval_fstat`.
4. CLI `n`/`coefsize` header formula was wrong (np*(np+1)/2); now full_stats size.

(qpWave/qpAdm fixes — eigenvector transpose, two regularizers, solvitforcez panic —
are documented in the changelog.)
