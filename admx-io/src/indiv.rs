//! AdmixTools `.ind` format parser.
//!
//! Ported from `convertf-rs/src/meta/ind.rs`.
//!
//! Three whitespace-separated columns: `sample_id  sex  population`

use crate::split_lines;
use admx_core::types::Sex;
use memmap2::Mmap;
use std::fs::File;
use std::io;
use std::path::Path;

/// One parsed individual row.
#[derive(Debug, Clone)]
pub struct IndRow {
    pub id: String,
    pub sex: Sex,
    pub pop: String,
    /// `true` if population was literally "Ignore" (case-insensitive).
    pub ignore: bool,
}

/// Parse a `.ind` file. Returns rows in file order.
pub fn read(path: &Path) -> io::Result<Vec<IndRow>> {
    let file = File::open(path)?;
    if file.metadata()?.len() == 0 { return Ok(Vec::new()); }
    let mmap = unsafe { Mmap::map(&file) }?;

    let mut rows = Vec::new();
    for (lineno, line) in split_lines(&mmap).enumerate() {
        let lineno = lineno + 1;
        if line.iter().all(|&b| b.is_ascii_whitespace()) { continue; }

        let row = parse_ind_line(line)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData,
                format!("{}:{}: {}", path.display(), lineno, e)))?;
        rows.push(row);
    }
    Ok(rows)
}

fn parse_ind_line(line: &[u8]) -> Result<IndRow, String> {
    let mut cols = line.split(|b: &u8| b.is_ascii_whitespace())
        .filter(|c| !c.is_empty());

    let id = cols.next().ok_or("missing sample id")?;
    let sex_raw = cols.next().ok_or("missing sex")?;
    let pop = cols.next().ok_or("missing population")?;

    let id = std::str::from_utf8(id).map_err(|e| format!("bad id: {e}"))?.to_owned();
    let pop = std::str::from_utf8(pop).map_err(|e| format!("bad pop: {e}"))?.to_owned();
    let sex = Sex::from_char(*sex_raw.first().unwrap_or(&b'U'));
    let ignore = pop.eq_ignore_ascii_case("Ignore");

    Ok(IndRow { id, sex, pop, ignore })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_basic_line() {
        let line = b"S_French-1.DG M French.DG";
        let row = parse_ind_line(line).unwrap();
        assert_eq!(row.id, "S_French-1.DG");
        assert_eq!(row.sex, Sex::Male);
        assert_eq!(row.pop, "French.DG");
        assert!(!row.ignore);
    }

    #[test]
    fn detects_ignore() {
        let line = b"sample1 U Ignore";
        let row = parse_ind_line(line).unwrap();
        assert!(row.ignore);
    }
}
