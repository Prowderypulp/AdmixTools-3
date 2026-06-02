//! EIGENSTRAT text genotype reader.
//!
//! Ported from `convertf-rs/src/geno/eigenstrat.rs`.
//!
//! - SNP-major. One line per SNP, one ASCII char per sample, `\n` terminator.
//! - Chars: `'0'`, `'1'`, `'2'`, `'9'` (missing). No header record.

use crate::{GenoReader, Layout, Storage, open_storage};
use std::io;
use std::path::Path;

/// LUT: ASCII byte → 2-bit genotype. Non-{0,1,2} → 0b11 (missing).
static ASCII_TO_2BIT: [u8; 256] = build_ascii_lut();

const fn build_ascii_lut() -> [u8; 256] {
    let mut t = [0b11u8; 256];
    t[b'0' as usize] = 0b00;
    t[b'1' as usize] = 0b01;
    t[b'2' as usize] = 0b10;
    t
}

pub struct EigenstratReader {
    storage: Storage,
    nind: usize,
    nsnp: usize,
    rec_bytes: usize,
    line_starts: Vec<usize>,
    next_idx: usize,
}

impl EigenstratReader {
    pub fn open(path: &Path, nind: usize, nsnp: usize) -> io::Result<Self> {
        let storage = open_storage(path)?;
        if storage.is_empty() {
            if nsnp == 0 {
                return Ok(Self {
                    storage,
                    nind, nsnp, rec_bytes: 0,
                    line_starts: Vec::new(), next_idx: 0,
                });
            }
            return Err(io::Error::new(io::ErrorKind::InvalidData,
                format!("EIGENSTRAT {}: empty file but .snp has {} rows", path.display(), nsnp)));
        }

        let line_starts = index_lines(&storage, nind, nsnp)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData,
                format!("{}: {}", path.display(), e)))?;

        let rec_bytes = (nind * 2 + 7) / 8;
        Ok(Self { storage, nind, nsnp, rec_bytes, line_starts, next_idx: 0 })
    }
}

fn index_lines(bytes: &[u8], nind: usize, expected_nsnp: usize) -> Result<Vec<usize>, String> {
    let mut starts = Vec::with_capacity(expected_nsnp);
    let mut cursor = 0usize;
    let len = bytes.len();

    while cursor < len {
        if bytes[cursor] == b'\n' { cursor += 1; continue; }
        if bytes[cursor] == b'\r' && cursor + 1 < len && bytes[cursor + 1] == b'\n' {
            cursor += 2; continue;
        }

        let line_start = cursor;
        let nl = memchr::memchr(b'\n', &bytes[cursor..])
            .map(|off| cursor + off)
            .unwrap_or(len);

        let content_end = if nl > line_start && bytes[nl.saturating_sub(1)] == b'\r' {
            nl - 1
        } else {
            nl
        };
        let content_len = content_end - line_start;

        if content_len != nind {
            return Err(format!(
                "EIGENSTRAT line {} has {} chars, expected {} (one per sample)",
                starts.len() + 1, content_len, nind));
        }
        if starts.len() >= expected_nsnp {
            return Err(format!("EIGENSTRAT has more lines than .snp has rows ({})", expected_nsnp));
        }
        starts.push(line_start);
        cursor = nl + 1;
    }

    if starts.len() != expected_nsnp {
        return Err(format!(
            "EIGENSTRAT line count {} != .snp row count {}",
            starts.len(), expected_nsnp));
    }
    Ok(starts)
}

impl GenoReader for EigenstratReader {
    fn nind(&self) -> usize { self.nind }
    fn nsnp(&self) -> usize { self.nsnp }
    fn layout(&self) -> Layout { Layout::SnpMajor }
    fn record_bytes(&self) -> usize { self.rec_bytes }

    fn read_record(&mut self, dst: &mut [u8]) -> io::Result<bool> {
        if self.next_idx >= self.nsnp { return Ok(false); }
        if dst.len() != self.rec_bytes {
            return Err(io::Error::new(io::ErrorKind::InvalidInput,
                format!("dst len {} != record_bytes {}", dst.len(), self.rec_bytes)));
        }
        let start = self.line_starts[self.next_idx];
        let end = start + self.nind;
        let line = &self.storage[start..end];

        for b in dst.iter_mut() { *b = 0; }
        let full = self.nind / 4;
        for i in 0..full {
            let off = i * 4;
            let b = (ASCII_TO_2BIT[line[off    ] as usize] << 6)
                  | (ASCII_TO_2BIT[line[off + 1] as usize] << 4)
                  | (ASCII_TO_2BIT[line[off + 2] as usize] << 2)
                  |  ASCII_TO_2BIT[line[off + 3] as usize];
            dst[i] = b;
        }
        let tail = self.nind % 4;
        if tail > 0 {
            let off = full * 4;
            let mut b = 0u8;
            for k in 0..tail {
                b |= ASCII_TO_2BIT[line[off + k] as usize] << (6 - 2 * k);
            }
            dst[full] = b;
        }

        self.next_idx += 1;
        Ok(true)
    }
}
