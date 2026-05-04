//! AdmixTools `.snp` format parser.
//!
//! Ported from `convertf-rs/src/meta/snp.rs`.
//!
//! Six whitespace-separated columns:
//! `snp_id  chrom  gen_pos  phys_pos  allele1  allele2`

use crate::split_lines;
use memmap2::Mmap;
use std::fs::File;
use std::io;
use std::path::Path;

/// One parsed SNP row.
#[derive(Debug, Clone)]
pub struct SnpRow {
    pub id: String,
    pub chrom: u8,
    pub genetic_pos: f64,
    pub physical_pos: u64,
    pub allele1: u8,
    pub allele2: u8,
}

/// Parse a `.snp` file. Returns rows in file order.
pub fn read(path: &Path, numchrom: u32) -> io::Result<Vec<SnpRow>> {
    let file = File::open(path)?;
    if file.metadata()?.len() == 0 { return Ok(Vec::new()); }
    let mmap = unsafe { Mmap::map(&file) }?;

    let mut rows = Vec::new();
    for (lineno, line) in split_lines(&mmap).enumerate() {
        let lineno = lineno + 1;
        if line.iter().all(|&b| b.is_ascii_whitespace()) { continue; }

        let row = parse_snp_line(line, numchrom)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData,
                format!("{}:{}: {}", path.display(), lineno, e)))?;
        rows.push(row);
    }
    Ok(rows)
}

fn parse_snp_line(line: &[u8], numchrom: u32) -> Result<SnpRow, String> {
    let mut cols = line.split(|b: &u8| b.is_ascii_whitespace())
        .filter(|c| !c.is_empty());

    let id = cols.next().ok_or("missing snp id")?;
    let chrom_raw = cols.next().ok_or("missing chrom")?;
    let gen_pos = cols.next().ok_or("missing genetic pos")?;
    let phys_pos = cols.next().ok_or("missing physical pos")?;
    let a1 = cols.next();
    let a2 = cols.next();

    let id = std::str::from_utf8(id).map_err(|e| format!("bad id: {e}"))?;
    let chrom = parse_chrom(chrom_raw, numchrom)?;

    let gen_pos: f64 = std::str::from_utf8(gen_pos)
        .map_err(|e| format!("bad gen_pos: {e}"))?
        .parse().map_err(|e| format!("bad gen_pos: {e}"))?;
    let phys_pos: u64 = std::str::from_utf8(phys_pos)
        .map_err(|e| format!("bad phys_pos: {e}"))?
        .parse().map_err(|e| format!("bad phys_pos: {e}"))?;

    let (allele1, allele2) = match (a1, a2) {
        (Some(a), Some(b)) if a.len() == 1 && b.len() == 1 => (a[0], b[0]),
        (None, None) => (b'X', b'X'),
        _ => return Err("expected 4 or 6 whitespace columns".to_string()),
    };

    Ok(SnpRow {
        id: id.to_owned(),
        chrom, genetic_pos: gen_pos, physical_pos: phys_pos,
        allele1, allele2,
    })
}

fn parse_chrom(raw: &[u8], numchrom: u32) -> Result<u8, String> {
    let s = std::str::from_utf8(raw).map_err(|e| format!("bad chrom: {e}"))?;
    let up = s.to_ascii_uppercase();
    let v: u32 = match up.as_str() {
        "X"  => numchrom + 1,
        "Y"  => numchrom + 2,
        "MT" | "M" => numchrom + 3,
        "XY" => numchrom + 4,
        num  => num.parse().map_err(|e| format!("bad chrom: {e}"))?,
    };
    if v > u8::MAX as u32 { return Err(format!("chrom {} out of u8 range", v)); }
    Ok(v as u8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic_line() {
        let line = b"rs1 1 0.001 752566 A G";
        let row = parse_snp_line(line, 22).unwrap();
        assert_eq!(row.id, "rs1");
        assert_eq!(row.chrom, 1);
        assert_eq!(row.physical_pos, 752566);
        assert_eq!(row.allele1, b'A');
        assert_eq!(row.allele2, b'G');
    }

    #[test]
    fn handles_x_chrom() {
        let line = b"rsX X 0.0 1000 A T";
        let row = parse_snp_line(line, 22).unwrap();
        assert_eq!(row.chrom, 23);
    }
}
