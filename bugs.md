# AdmixTools-rs — Bugs & Fidelity Gaps

Snapshot: 2026-04-25. Against C reference in `../src/` and `RUST_PORT_PHASE2.md` gates.

Severity: **P0** = will produce wrong numbers on real data; **P1** = blocks a milestone gate or contaminates output; **P2** = fidelity gap, no wrong final answer on typical inputs.

Resolved entries are in `changelog.md`.

---

## Open

### P0

**`ranktestfix` — wrong `vl` constraint set (admx-rank/src/ranktest.rs:403–413)**

The Rust `vl` construction has an extra filter `if vfix[j] == 0 { continue; }` that is absent in the C:

```c
// C (f4rank.c:276-289)
for (j = 0; j < m; ++j) {
    if (i == j) continue;
    vl[nf] = j * rank + k;   // ALL j != i, regardless of vfix[j]
    ++nf;
}
```

```rust
// Rust — WRONG
for j in 0..m {
    if i == j { continue; }
    if vfix[j] == 0 { continue; }  // extra filter, not in C
    vl.push(j * rank + k_idx);
}
```

Effect: fewer zero-constraints are applied to A during the fixed-rank ALS solve. For a model with any non-fixed populations (vfix[j]=0), the corresponding A column entries are not forced to zero. Produces wrong admixture-related A matrices and chi-square values for every `ranktestfix`/`doranktestfix` call. qpAdm fixed-population models will silently produce wrong p-values.

Fix: remove the `if vfix[j] == 0 { continue; }` guard.

---

### P1

**`ranktest` — LAPACK packed-matrix layout mismatch (admx-rank/src/ranktest.rs:235–243)**

`dspev_l` (uplo=`'L'`) expects the lower triangle packed **column-by-column** (LAPACK convention). The Rust packing loop iterates row-by-row:

```rust
// Rust — WRONG (row-wise lower triangle)
for i in 0..n {
    for j in 0..=i {
        wright_packed[idx] = -wright[i * n + j];
        idx += 1;
    }
}
```

For `n ≥ 3` this swaps elements: LAPACK reads `a[2,0]` from the slot that holds `a[1,1]` and vice versa (and similarly for larger `n`). Result: the initial B vectors used to seed the ALS are the eigenvectors of the wrong matrix.

For `n = 2` the row-wise and column-wise orderings coincide, so the bug is latent.

Effect: after 20 ALS iterations the solver usually finds the same fixed point, but for ill-conditioned inputs or a poor initialization it can converge to a wrong local minimum, producing an incorrect chi-square value. The correct packing loop is:

```rust
for j in 0..n {
    for i in j..n {
        wright_packed[idx] = -wright[i * n + j];  // column-wise
        idx += 1;
    }
}
```

**DEBUG print in `ranktest` rank=0 path (admx-rank/src/ranktest.rs:196)**

```rust
println!("DEBUG varinv Rust: {:20.12e} {:20.12e}", varinv[0], varinv[1]);
```

This line fires on every qpWave rank-0 test and appears verbatim in the log. It will cause the golden-log regression in `admx-cli/tests/golden_log.rs` to fail and contaminates production output. Remove before any validation run.

**DEBUG prints in qpwave driver (admx-qpwave/src/driver.rs:106–111)**

```rust
println!("DEBUG ymean: {:.6} {:.6}", ymean[0], ymean[1]);
print!("DEBUG yvar: ");
for i in 0..yvar.len() { print!("{:20.12e} ", yvar[i]); }
println!();
```

Same issue as above — contaminates all qpWave output. Remove.

**`admx-qpadm` milestone entirely unimplemented**

`admx-qpadm/src/driver.rs`, `calcadm.rs`, and `bootstrap.rs` are single-line doc stubs. `admx-cli/src/bin/qpAdm.rs` prints "qpAdm stub". The qpAdm milestone gate (`RUST_PORT_PHASE2.md` M7) is fully blocked.

Functions needed:
- `calcadm` / `calcadmfix` — weight solver from `qpAdm.c:calcadm`
- `calcevarboot` — bootstrap covariance from `qpAdm.c`
- `doq4vecb` / `loadymv` caller loop in the driver
- CLI: parfile parsing, output formatting

**`qpwave` direct-genotype path unimplemented (admx-qpwave/src/driver.rs:141–144)**

When `fstatsname:` is absent, the driver immediately returns an error:

```rust
Err(AdmxError::Fatal("qpWave without fstatsname not yet fully implemented".into()))
```

The C `qpWave` can run directly from genotype files (calling `doq4vecb` internally). This path is needed for any run that does not pre-compute an fstats file.

---

### P2

**`prng.rs` LCG constants are placeholders (admx-qpadm/src/prng.rs:26)**

```rust
// TODO(M7): Port the exact constants from the C LCG.
// Placeholder constants — replace with the real ones from nicksrc.
self.state = self.state.wrapping_mul(6364136223846793005).wrapping_add(1);
```

The `SRAND`/`ranmod` LCG used in `qpAdm.c` for `calcevarboot` has different multiplier/increment. Until ported, bootstrap covariance estimates will not be bit-identical to C for the same seed.

**`rtlchsq` extreme-tail branch absent (admx-rank/src/stats.rs)**

The C code switches to `rtlg` (log-space computation) when the tail probability drops below ~1e-6. The Rust `pochisq` is used throughout. For very small p-values (large chi-square or high df) the result may underflow to zero where the C would return a finite small value. Not a problem for typical runs but matters for reporting highly significant differences.

**`ranktestfix` uses uniform RNG instead of Gaussian for A/B initialisation (admx-rank/src/ranktest.rs:441–448)**

C uses `gaussa(A, adim)` (standard-normal draws from `nicklib`). Rust uses a simple LCG scaled to `[-1, 1]`. Initialization is uniform, not Gaussian. After 50 ALS iterations the final fixed point is normally the same, but for degenerate inputs convergence may be slower or land on a different local minimum.

April 25

1. High: the workspace does not currently build because both CLI
     binaries call a nonexistent ParFile::get_int(key, default)
     overload. ParFile only exposes get_int(&self, key) -> Option<i32>
     in admx-io/src/params.rs:133, but admx-cli/src/bin/qpAdm.rs:34 and
     admx-cli/src/bin/qpWave.rs:71 pass a second default argument and
     then cast the Option<i32>. cargo test --workspace fails here before
     any runtime validation.
  2. High: run_qpwave never applies the legacy diagonal regularizer to
     yvar before checkmv and doranktest. The Rust path builds yvar and
     immediately calls checkmv in admx-qpwave/src/driver.rs:181, but the
     C reference does if (diagvarplus < 0) diagvarplus = yscale;
     addscaldiag(yvar, diagvarplus, nl*nr); first in /home/drtex/Code/
     AdmixTools/src/qpWave.c:467. That changes both abort behavior and
     the matrix the rank test sees, so the Rust qpWave statistics are
     not fidelity-equivalent even if the solver itself is correct.
  3. Medium: both CLIs print the parameter banner in the wrong place, so
     logs will differ from the Phase 2 golden format. The Rust code
     emits ##PARAMETER NAME: VALUE inside the per-line loop in admx-cli/
     src/bin/qpWave.rs:49 and admx-cli/src/bin/qpAdm.rs:17, which
     repeats the banner once per parameter. The C binaries print that
     header once, then dump the full parameter set via writepars(ph)
     in /home/drtex/Code/AdmixTools/src/qpWave.c:586 and /home/drtex/
     Code/AdmixTools/src/qpAdm.c:1246.
  4. Medium: qpWave's direct-genotype path still ignores fancyf4. The
     flag is parsed into admx-qpwave/src/driver.rs:17, but the direct
     path just routes through run_qpfstats with a admx-qpwave/src/
     driver.rs:112 that has no fancyf4 equivalent, and config.fancyf4 is
     never used in the driver. In C, qpWave sets this explicitly before
     doq4vecb, so fancyf4: NO currently behaves like YES.
  5. Medium: qpAdm is still a stub, which remains a hard Phase 2
     blocker. admx-qpadm/src/driver.rs:12 unconditionally returns "qpAdm
     milestone entirely unimplemented", while the CLI in admx-cli/src/
     bin/qpAdm.rs:36 already advertises and parses parameters as if the
     feature exists.

  The bugs.md items around dspev packing, ranktestfix’s vl construction,
  and rtlchsq’s extreme-tail branch look fixed in the current Rust tree.
