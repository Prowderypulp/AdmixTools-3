//! PLINK `.bed` format reader (PACKEDPED, SNP-major).
//!
//! Ported from `convertf-rs/src/geno/packed_ped.rs`.
//!
//! PLINK uses a different 2-bit convention and LSB-first bit order.
//! We recode to canonical AdmixTools encoding via a 256-entry LUT.

use crate::{GenoReader, Layout, Storage, open_storage};
use std::io;
use std::path::Path;

pub const BED_MAGIC: [u8; 3] = [0x6c, 0x1b, 0x01];

/// LUT: PLINK byte → AdmixTools canonical byte.
static PLINK_TO_AM: [u8; 256] = build_plink_to_am();

const fn recode_plink_to_am(two: u8) -> u8 {
    match two & 0b11 { 0b00 => 0b10, 0b01 => 0b11, 0b10 => 0b01, _ => 0b00 }
}

const fn build_plink_to_am() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        let s0 = b & 0b11;
        let s1 = (b >> 2) & 0b11;
        let s2 = (b >> 4) & 0b11;
        let s3 = (b >> 6) & 0b11;
        let out = (recode_plink_to_am(s0) << 6)
                | (recode_plink_to_am(s1) << 4)
                | (recode_plink_to_am(s2) << 2)
                |  recode_plink_to_am(s3);
        t[i] = out;
        i += 1;
    }
    t
}

pub struct PackedPedReader {
    storage: Storage,
    nind: usize,
    nsnp: usize,
    plink_rec_bytes: usize,
    am_rec_bytes: usize,
    last_byte_valid_samples: usize,
    next_idx: usize,
}

impl PackedPedReader {
    pub fn open(path: &Path, nind: usize, nsnp: usize) -> io::Result<Self> {
        let storage = open_storage(path)?;

        if storage.len() < 3 || storage[0] != BED_MAGIC[0] || storage[1] != BED_MAGIC[1] {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("{}: not a PLINK .bed file (bad magic)", path.display())));
        }
        if storage[2] != 0x01 {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("{}: sample-major .bed not supported", path.display())));
        }

        let plink_rec_bytes = (nind + 3) / 4;
        let expected_len = 3 + plink_rec_bytes * nsnp;
        if storage.len() < expected_len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("{}: file size {} < expected {}", path.display(), storage.len(), expected_len)));
        }

        let am_rec_bytes = (nind * 2 + 7) / 8;
        debug_assert_eq!(plink_rec_bytes, am_rec_bytes);

        let last_byte_valid_samples = if nind > 0 { ((nind - 1) % 4) + 1 } else { 0 };

        Ok(Self {
            storage, nind, nsnp, plink_rec_bytes, am_rec_bytes,
            last_byte_valid_samples, next_idx: 0,
        })
    }
}

impl GenoReader for PackedPedReader {
    fn nind(&self) -> usize { self.nind }
    fn nsnp(&self) -> usize { self.nsnp }
    fn layout(&self) -> Layout { Layout::SnpMajor }
    fn record_bytes(&self) -> usize { self.am_rec_bytes }

    fn read_record(&mut self, dst: &mut [u8]) -> io::Result<bool> {
        if self.next_idx >= self.nsnp { return Ok(false); }
        if dst.len() != self.am_rec_bytes {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                format!("dst len {} != record_bytes {}", dst.len(), self.am_rec_bytes)));
        }
        let start = 3 + self.plink_rec_bytes * self.next_idx;
        let end = start + self.plink_rec_bytes;
        let src = &self.storage[start..end];

        for (i, &b) in src.iter().enumerate() {
            dst[i] = PLINK_TO_AM[b as usize];
        }

        // Mask trailing padding bits in the last byte.
        if self.nind % 4 != 0 {
            let valid = self.last_byte_valid_samples;
            let keep_bits = valid * 2;
            let mask = 0xFFu8 << (8 - keep_bits);
            let last = self.am_rec_bytes - 1;
            dst[last] &= mask;
        }

        self.next_idx += 1;
        Ok(true)
    }
}
