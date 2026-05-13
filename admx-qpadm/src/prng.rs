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

/// Legacy LCG state.  Matches the C global seed.
pub struct LegacyLcg {
    _dummy: (),
}

impl LegacyLcg {
    /// Seed the generator (mirrors C `SRAND(seed)`).
    pub fn new(seed: u64) -> Self {
        unsafe { srandom(seed as u32); }
        Self { _dummy: () }
    }

    /// Generate the next pseudo-random value in [0, 1) using a single `random()` call.
    /// This matches the C `DRAND()` macro used by `gauss()`.
    pub fn next_f64_drand(&mut self) -> f64 {
        unsafe {
            let r = random() % (i32::MAX as std::ffi::c_long);
            (r as f64) / (i32::MAX as f64)
        }
    }

    /// Generate the next pseudo-random value in [0, 1).
    ///
    /// Must match the C `ranmod()` output sequence exactly.
    pub fn next_f64(&mut self) -> f64 {
        unsafe {
            let maxran = 1.0 - std::f64::EPSILON;
            let maxran1 = (i32::MAX as f64 - 1.0) / (i32::MAX as f64);
            let eps = maxran - maxran1;
            
            let r1 = random() % (i32::MAX as std::ffi::c_long);
            let r2 = random() % (i32::MAX as std::ffi::c_long);
            
            let x = (r1 as f64) / (i32::MAX as f64);
            let y = (r2 as f64) / (i32::MAX as f64);
            
            x + y * eps
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lcg_deterministic() {
        let mut rng1 = LegacyLcg::new(42);
        let val1 = rng1.next_f64();
        let mut rng2 = LegacyLcg::new(42);
        let val2 = rng2.next_f64();
        assert_eq!(val1.to_bits(), val2.to_bits());
    }
}
