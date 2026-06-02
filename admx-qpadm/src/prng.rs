//! Legacy linear-congruential PRNG.
//!
//! Byte-for-byte port of the C `SRAND` / `ranmod` used in `qpAdm.c` so that
//! `numboot:` bootstrap covariances are bit-identical to the legacy output.
//!
//! This is intentionally NOT a cryptographically secure RNG.  It exists solely
//! to reproduce the legacy sequence.

#[cfg(target_os = "linux")]
extern "C" {
    fn srandom(seed: u32);
    fn random() -> std::ffi::c_long;
}

#[cfg(not(target_os = "linux"))]
extern "C" {
    // macOS also has srandom and random
    fn srandom(seed: u32);
    fn random() -> std::ffi::c_long;
}

/// Serializes every test that touches the process-global glibc `random()`
/// stream.  Rust runs tests in parallel threads, so without this guard two RNG
/// tests interleave their `srandom`/`random` calls and corrupt each other's
/// sequence.  Hold the guard for the whole RNG-using body.
#[cfg(test)]
pub(crate) static RNG_TEST_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Legacy LCG state.  Matches the C global seed.
///
/// The underlying uniform stream is glibc's process-global `random()`, so this
/// holds no uniform state itself — only the Marsaglia-polar normal cache
/// (`iset`/`gset`) that C's `gauss()` keeps in function-static storage
/// (`nicksrc/gauss.c`).  Only one instance may be live at a time.
pub struct LegacyLcg {
    /// Mirror of C `gauss()`'s static `iset`: a cached normal is pending.
    iset: bool,
    /// Mirror of C `gauss()`'s static `gset`: the cached normal value.
    gset: f64,
}

impl LegacyLcg {
    /// Seed the generator (mirrors C `SRAND(seed)`).
    pub fn new(seed: u64) -> Self {
        unsafe { srandom(seed as u32); }
        Self { iset: false, gset: 0.0 }
    }

    /// Generate the next uniform in [0, 1).  Byte-for-byte port of C `drand2()`
    /// (`nicksrc/gds.c`): two `DRAND()` draws combined for full mantissa
    /// precision, `DRAND() = (random() % BIGINT) / BIGINT`, `BIGINT = INT_MAX`.
    pub fn next_f64(&mut self) -> f64 {
        unsafe {
            let maxran = 1.0 - f64::EPSILON;
            let maxran1 = (i32::MAX as f64 - 1.0) / (i32::MAX as f64);
            let eps = maxran - maxran1;

            let r1 = random() % (i32::MAX as std::ffi::c_long);
            let r2 = random() % (i32::MAX as std::ffi::c_long);

            let x = (r1 as f64) / (i32::MAX as f64);
            let y = (r2 as f64) / (i32::MAX as f64);

            x + y * eps
        }
    }

    /// Standard normal — exact port of C `gauss()` (`nicksrc/gauss.c`), the
    /// Marsaglia *polar* method (NR in C, pp. 289ff), NOT Box-Muller.  It draws
    /// a pair `(v1, v2)` uniform on the unit square, rejects until they fall in
    /// the open unit disc, returns `v2*fac` now and caches `v1*fac` in `gset`
    /// for the next call.  The cache is what makes the consumed-uniform count
    /// and the returned sequence match C bit-for-bit, so it must persist across
    /// every `gauss()` call within a run.
    pub fn gauss(&mut self) -> f64 {
        if self.iset {
            self.iset = false;
            return self.gset;
        }
        loop {
            let v1 = 2.0 * self.next_f64() - 1.0;
            let v2 = 2.0 * self.next_f64() - 1.0;
            let rsq = v1 * v1 + v2 * v2;
            if rsq < 1.0 && rsq != 0.0 {
                let fac = (-2.0 * rsq.ln() / rsq).sqrt();
                self.gset = v1 * fac;
                self.iset = true;
                return v2 * fac;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lcg_deterministic() {
        let _guard = RNG_TEST_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let mut rng1 = LegacyLcg::new(42);
        let val1 = rng1.next_f64();
        let mut rng2 = LegacyLcg::new(42);
        let val2 = rng2.next_f64();
        assert_eq!(val1.to_bits(), val2.to_bits());
    }
}
