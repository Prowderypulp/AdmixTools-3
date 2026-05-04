//! PACKEDANCESTRYMAP format reader.
//!
//! Ported from `convertf-rs/src/geno/packed_am.rs`.
//!
//! - SNP-major. One record per SNP.
//! - `rlen = max(48, ceil(nind * 2 / 8))` bytes per record.
//! - First record is header: `"GENO %d %d %x %x"`.
//! - Subsequent records: 2-bit MSB-first canonical encoding.

use crate::{GenoReader, Layout, Storage, open_storage};
use std::io;
use std::path::Path;

pub const MIN_RLEN: usize = 48;

#[inline]
pub fn rlen_for(nind: usize) -> usize {
    std::cmp::max(MIN_RLEN, (nind * 2 + 7) / 8)
}

/// Memory-mapped reader for PACKEDANCESTRYMAP `.geno` files.
pub struct PackedAmReader {
    storage: Storage,
    nind: usize,
    nsnp: usize,
    rlen: usize,
    rec_bytes: usize,
    next_idx: usize,
    pub ihash: u32,
    pub shash: u32,
}

impl PackedAmReader {
    /// Open a PACKEDANCESTRYMAP `.geno` file. `nind` and `nsnp` must come
    /// from the companion `.ind` / `.snp` files.
    pub fn open(path: &Path, nind: usize, nsnp: usize) -> io::Result<Self> {
        let storage = open_storage(path)?;

        let rlen = rlen_for(nind);
        let expected_len = rlen * (nsnp + 1);
        if storage.len() < expected_len {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("PACKEDANCESTRYMAP {}: file size {} < expected {} (rlen={}, nsnp={}, nind={})",
                    path.display(), storage.len(), expected_len, rlen, nsnp, nind)));
        }
        if storage.len() > expected_len {
            log::warn!(
                "PACKEDANCESTRYMAP {}: {} trailing bytes past expected end",
                path.display(), storage.len() - expected_len
            );
        }

        let (hdr_nind, hdr_nsnp, ihash, shash) = parse_header(&storage[..rlen])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData,
                format!("parse header of {}: {}", path.display(), e)))?;

        if hdr_nind != nind {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("PACKEDANCESTRYMAP header nind {} != .ind count {}", hdr_nind, nind)));
        }
        if hdr_nsnp != nsnp {
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("PACKEDANCESTRYMAP header nsnp {} != .snp count {}", hdr_nsnp, nsnp)));
        }

        let rec_bytes = (nind * 2 + 7) / 8;
        Ok(Self { storage, nind, nsnp, rlen, rec_bytes, next_idx: 0, ihash, shash })
    }
}

impl GenoReader for PackedAmReader {
    fn nind(&self) -> usize { self.nind }
    fn nsnp(&self) -> usize { self.nsnp }
    fn layout(&self) -> Layout { Layout::SnpMajor }
    fn record_bytes(&self) -> usize { self.rec_bytes }
    fn header_hashes(&self) -> Option<(u32, u32)> { Some((self.ihash, self.shash)) }

    fn read_record(&mut self, dst: &mut [u8]) -> io::Result<bool> {
        if self.next_idx >= self.nsnp { return Ok(false); }
        if dst.len() != self.rec_bytes {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                format!("dst len {} != record_bytes {}", dst.len(), self.rec_bytes)));
        }
        let start = self.rlen * (1 + self.next_idx);  // skip header record
        let end = start + self.rec_bytes;
        dst.copy_from_slice(&self.storage[start..end]);
        self.next_idx += 1;
        Ok(true)
    }
}

const HEADER_MAGIC: &[u8] = b"GENO ";

/// Parse `"GENO %d %d %x %x"` from the header record.
fn parse_header(hdr: &[u8]) -> Result<(usize, usize, u32, u32), String> {
    if !hdr.starts_with(HEADER_MAGIC) {
        return Err("not a PACKEDANCESTRYMAP file (missing 'GENO ' magic)".to_string());
    }
    let end = hdr.iter().position(|&b| b == 0).unwrap_or(hdr.len());
    let text = std::str::from_utf8(&hdr[..end])
        .map_err(|e| format!("header not UTF-8: {e}"))?;
    let mut it = text.split_ascii_whitespace();
    let _geno = it.next();
    let nind: usize = it.next().ok_or("header: missing nind")?
        .parse().map_err(|e| format!("header nind: {e}"))?;
    let nsnp: usize = it.next().ok_or("header: missing nsnp")?
        .parse().map_err(|e| format!("header nsnp: {e}"))?;
    let ihash = u32::from_str_radix(
        it.next().ok_or("header: missing ihash")?, 16)
        .map_err(|e| format!("header ihash: {e}"))?;
    let shash = u32::from_str_radix(
        it.next().ok_or("header: missing shash")?, 16)
        .map_err(|e| format!("header shash: {e}"))?;
    Ok((nind, nsnp, ihash, shash))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rlen_min_48() {
        assert_eq!(rlen_for(1), 48);
        assert_eq!(rlen_for(192), 48);
        assert_eq!(rlen_for(193), 49);
        assert_eq!(rlen_for(4000), (4000 * 2 + 7) / 8);
    }

    #[test]
    fn parses_header() {
        let mut buf = vec![0u8; 48];
        let s = b"GENO 5 3 1a2b cafe";
        buf[..s.len()].copy_from_slice(s);
        let (nind, nsnp, ih, sh) = parse_header(&buf).unwrap();
        assert_eq!((nind, nsnp, ih, sh), (5, 3, 0x1a2b, 0xcafe));
    }

    #[test]
    fn rejects_wrong_magic() {
        let buf = vec![b'X'; 48];
        assert!(parse_header(&buf).is_err());
    }
}
