//! `admx-io` — Genotype, SNP, and individual file I/O.
//!
//! Ported from `convertf-rs` (which itself faithfully ports `mcio.c`).
//! Supports:
//! - Packed ANCESTRYMAP (`inpack`, `inpackt`)
//! - EIGENSTRAT text genotypes
//! - PLINK `.bed` (SNP-major only)
//! - `.ind` / `.snp` metadata
//! - Parameter file parsing (`getpars`)
//! - `.gz` transparent decompression via `flate2`

pub mod codec;
pub mod packed_am;
pub mod eigenstrat;
pub mod packed_ped;
pub mod snp;
pub mod indiv;
pub mod params;

use std::path::Path;
use std::fs::File;
use std::io::{self, Read};
use flate2::read::GzDecoder;

pub enum Storage {
    Mmap(memmap2::Mmap),
    Buffer(Vec<u8>),
}

impl std::ops::Deref for Storage {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        match self {
            Storage::Mmap(m) => m,
            Storage::Buffer(v) => v,
        }
    }
}

pub fn open_storage<P: AsRef<Path>>(path: P) -> io::Result<Storage> {
    let path = path.as_ref();
    let file = File::open(path)?;
    if path.extension().map_or(false, |ext| ext == "gz") {
        let mut decoder = GzDecoder::new(file);
        let mut buffer = Vec::new();
        decoder.read_to_end(&mut buffer)?;
        Ok(Storage::Buffer(buffer))
    } else {
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        Ok(Storage::Mmap(mmap))
    }
}

/// Record orientation (SNP-major vs sample-major).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Layout {
    /// One record per SNP, covering all samples.
    SnpMajor,
    /// One record per sample, covering all SNPs (TGENO — out of scope for Phase 1).
    SampleMajor,
}

/// Streaming genotype reader trait.
///
/// All implementations normalize to canonical 2-bit MSB-first encoding:
/// `00=0 (hom ref), 01=1 (het), 10=2 (hom alt), 11=missing`.
pub trait GenoReader {
    fn nind(&self) -> usize;
    fn nsnp(&self) -> usize;
    fn layout(&self) -> Layout;

    /// Bytes per record in canonical encoding: `ceil(n * 2 / 8)`.
    fn record_bytes(&self) -> usize {
        let n = match self.layout() {
            Layout::SnpMajor => self.nind(),
            Layout::SampleMajor => self.nsnp(),
        };
        (n * 2 + 7) / 8
    }

    /// Hashes from the geno file header (PACKEDANCESTRYMAP).
    /// Returns `None` for header-less formats (EIGENSTRAT, PLINK BED).
    fn header_hashes(&self) -> Option<(u32, u32)> { None }

    /// Read the next record into `dst`. Returns `Ok(false)` at EOF.
    fn read_record(&mut self, dst: &mut [u8]) -> Result<bool, std::io::Error>;

    /// Random-access view of the canonical record bytes, for readers backed by a
    /// contiguous in-memory or mmap buffer (PACKEDANCESTRYMAP, in-memory BED).
    /// Returns `(bytes, off0, stride, reclen)` such that SNP `idx`'s record is
    /// `&bytes[off0 + idx*stride .. off0 + idx*stride + reclen]`. The returned
    /// slice is shareable across threads (`&[u8]` is `Sync`), enabling parallel
    /// SNP-major scans. Streaming/text readers (EIGENSTRAT) return `None` and
    /// callers fall back to the sequential `read_record` path.
    fn random_access(&self) -> Option<(&[u8], usize, usize, usize)> { None }
}

/// Quick check if a file is PACKEDANCESTRYMAP.
pub fn is_packed_am<P: AsRef<Path>>(path: P) -> bool {
    let path = path.as_ref();
    let res = (|| {
        let file = File::open(path)?;
        let mut reader: Box<dyn Read> = if path.extension().map_or(false, |ext| ext == "gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut buf = [0u8; 5];
        reader.read_exact(&mut buf)?;
        Ok::<_, io::Error>(&buf == b"GENO ")
    })();
    res.unwrap_or(false)
}

/// Quick check if a file is PLINK `.bed`.
pub fn is_bed<P: AsRef<Path>>(path: P) -> bool {
    let path = path.as_ref();
    let res = (|| {
        let file = File::open(path)?;
        let mut reader: Box<dyn Read> = if path.extension().map_or(false, |ext| ext == "gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
        let mut buf = [0u8; 3];
        reader.read_exact(&mut buf)?;
        Ok::<_, io::Error>(buf[0] == 0x6c && buf[1] == 0x1b)
    })();
    res.unwrap_or(false)
}

/// Quick check if a file is EIGENSTRAT text.
/// 
/// Heuristic: first non-empty line has no high-bit chars and matches nind.
pub fn is_eigenstrat<P: AsRef<Path>>(path: P, nind: usize) -> bool {
    let storage = open_storage(path).ok();
    if let Some(s) = storage {
        let first_line = split_lines(&s).next();
        if let Some(line) = first_line {
            line.len() == nind && line.iter().all(|&b| b.is_ascii_digit())
        } else {
            false
        }
    } else {
        false
    }
}

/// Split mmap bytes at `\n` without allocating. Strips trailing `\r`.
pub(crate) fn split_lines(bytes: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut start = 0usize;
    let len = bytes.len();
    std::iter::from_fn(move || {
        if start >= len { return None; }
        let rest = &bytes[start..];
        match memchr::memchr(b'\n', rest) {
            Some(off) => {
                let mut end = start + off;
                let line_start = start;
                start = end + 1;
                if end > line_start && bytes[end - 1] == b'\r' { end -= 1; }
                Some(&bytes[line_start..end])
            }
            None => {
                let line_start = start;
                start = len;
                Some(&bytes[line_start..len])
            }
        }
    })
}
