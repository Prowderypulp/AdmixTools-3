# AdmixTools Rust Port (admixtools-rs)

This is the official high-performance Rust port of the legacy C `ADMIXTOOLS` toolkit. 
The codebase is currently in Phase 2 completion, meaning the core tools `qpWave` and `qpAdm` are fully functional and validated for numerical parity with the original C binaries.

## Implemented Tools

The following tools are available as binaries in this workspace:

- **qpWave**: Tests whether a set of target populations is consistent with being descended from $N$ waves of admixture from a set of source populations.
- **qpAdm**: Estimates the proportions of admixture in a target population from a set of source populations, leveraging $F_4$ statistics and block-jackknife for standard errors.

## Installation & Building

Make sure you have Rust (MSRV 1.75) installed. To build all tools with optimizations:

```bash
cargo build --release
```

The compiled binaries will be available in `target/release/`.

## Usage

Both `qpWave` and `qpAdm` are invoked by passing a parameter file (`parfile`). 
You can run them directly via Cargo or by executing the compiled binaries.

### Running qpWave

```bash
cargo run --release --bin qpWave -- -p <path_to_parfile>
```

**qpWave Parameter Options:**
- `fstatsname: <file>`: Use pre-computed $F$-statistics.
- `genotypename: <file>`, `snpname: <file>`, `indivname: <file>`: Compute statistics on the fly (direct-genotype fallback).
- `popleft: <list>`: List of target populations.
- `popright: <list>`: List of outgroup/reference populations.
- `allsnps: YES/NO`: (Default NO) Use all SNPs (mutually exclusive missingness).
- `fancyf4: YES/NO`: (Default YES) Use the $F_3$ canonical basis calculation for performance.

### Running qpAdm

```bash
cargo run --release --bin qpAdm -- -p <path_to_parfile>
```

**qpAdm Parameter Options:**
- `fstatsname: <file>`: Requires pre-computed $F$-statistics. Direct-genotype fallback is not yet supported.
- `popleft: <list>`: Target population (first element) and source populations.
- `popright: <list>`: List of outgroup/reference populations.
- `allsnps: YES/NO`: (Default NO) Use all SNPs.
- `numboot: <int>`: (Default 1000) Number of block-jackknife iterations for standard error computation.
- `seed: <int>`: Legacy PRNG seed for deterministic bootstrapping.

## Performance

`qpAdm` is benchmarked against the legacy C `qpAdm` binary (AdmixTools 8.0.1)
via `admx-bench`. End-to-end wallclock, 7 runs each after 2 warmup runs,
14-core host:

| tier  | n_snps  | n_inds | C median | Rust median | **speedup** |
|-------|--------:|-------:|---------:|------------:|------------:|
| small | 20,000  | 70     | 0.091 s  | 0.007 s     | **13.4×**   |
| med   | 100,000 | 200    | 0.264 s  | 0.014 s     | **18.3×**   |

Numerical output matches the C reference bit-for-bit on the golden-log harness
(`admx-cli/tests/golden_log.rs::test_qpadm_golden`). Two changes drove the
result, both implemented without altering the numerical pipeline:

1. **Single-threaded BLAS** at startup (`admx-core::linalg::set_blas_single_threaded`).
   qpAdm makes thousands of small (≤16×16) BLAS calls inside the bootstrap;
   OpenBLAS's OpenMP fan-out previously dominated runtime (profiling showed
   96.6% of samples in libgomp barriers).
2. **Rayon-parallel bootstrap** (`admx-qpadm/src/bootstrap.rs`). The 1000
   independent bootstrap iterations run on a dedicated thread pool whose
   workers each force single-threaded BLAS, so nested OMP threading does not
   regress the parallel section.

See `admx-qpadm/plan.md` for the full investigation, profile data, and design
rationale.

## Development and Phase 3

Phase 3 focus:
- `SIMD` optimizations for accelerated $F$-statistic accumulation in `qpfstats`
  (the §3.1 matrix-based accumulation is the next algorithmic win — see
  `admx-bench/plan.md`).
- Porting of `qpDstat` and `qp3Pop`.

For detailed architecture differences between the C and Rust codebases, refer to `RUST_PORT_PHASE2.md`.
