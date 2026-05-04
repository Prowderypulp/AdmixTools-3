//! Core types mirroring the C structs in `admutils.h`, `mcio.h`, and `qpsubs.h`.
//!
//! Field names intentionally match the C originals so the codebase can be
//! cross-referenced against the legacy source line-by-line.

/// Maximum length of an identifier string (matches C `IDSIZE = 129`).
pub const IDSIZE: usize = 129;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
    Unknown,
}

impl Sex {
    pub fn from_char(c: u8) -> Self {
        match c {
            b'M' | b'm' | b'1' => Sex::Male,
            b'F' | b'f' | b'2' => Sex::Female,
            _ => Sex::Unknown,
        }
    }

    pub fn as_ind_char(self) -> char {
        match self { Sex::Male => 'M', Sex::Female => 'F', Sex::Unknown => 'U' }
    }
}

/// Individual / sample record.  Mirrors `Indiv` in `admutils.h`.
#[derive(Debug, Clone)]
pub struct Indiv {
    /// Sample identifier (max `IDSIZE` chars).
    pub id: String,
    /// Population / ethnic-group label (the `egroup` pointer in C).
    pub egroup: String,
    /// Gender: Male, Female, or Unknown.
    pub sex: Sex,
    /// Numeric index assigned during loading (`idnum` in C).
    pub idnum: usize,
    /// Affected status (case/control).
    pub affstatus: i32,
    /// If `true`, this individual is excluded from analysis.
    pub ignore: bool,
    /// Population-group code used in per-SNP loops (`xtypes[i]` in C).
    /// Negative ⇒ ignore.
    pub gkode: i32,
}

/// SNP (single-nucleotide polymorphism) record.  Mirrors `SNP` in `admutils.h`.
#[derive(Debug, Clone)]
pub struct Snp {
    /// SNP identifier.
    pub id: String,
    /// Chromosome number.
    pub chrom: i32,
    /// Chromosome as string (e.g., "X", "MT").
    pub cchrom: String,
    /// Genetic position (Morgans).
    pub genpos: f64,
    /// Physical position (bp).
    pub physpos: f64,
    /// Reference and alternate alleles.
    pub alleles: [char; 2],
    /// Whether this SNP is ignored.
    pub ignore: bool,
    /// Jackknife block assignment (`tagnumber` in C).
    pub tagnumber: i32,
    /// Weight for this SNP.
    pub weight: f64,
}

/// Per-population allele counts for a single SNP.
///
/// In the C code these are computed by `loadaa` / `countpops` and stored as
/// `int counts[ncols][numeg][2]`.  We use a struct for clarity.
#[derive(Debug, Clone, Copy, Default)]
pub struct PopCounts {
    /// Count of reference alleles observed.
    pub nref: u32,
    /// Count of alternate alleles observed.
    pub nalt: u32,
}

/// Detailed per-population genotype counts for a single SNP.
/// 
/// Matches the 5-element array used in `loadaa` in `qpsubs.c`:
/// - `cc[0]`: hom ref females (diploid)
/// - `cc[1]`: het females (diploid)
/// - `cc[2]`: hom alt females (diploid)
/// - `cc[3]`: ref males (haploid / X)
/// - `cc[4]`: alt males (haploid / X)
#[derive(Debug, Clone, Copy, Default)]
pub struct DetailedPopCounts {
    pub f0: u32,
    pub f1: u32,
    pub f2: u32,
    pub m0: u32,
    pub m1: u32,
}

impl DetailedPopCounts {
    /// Convert to basic PopCounts (total allele counts).
    pub fn to_pop_counts(&self) -> PopCounts {
        PopCounts {
            nref: 2 * self.f0 + self.f1 + self.m0,
            nalt: 2 * self.f2 + self.f1 + self.m1,
        }
    }
}

impl PopCounts {
    /// Total non-missing allele count.
    #[inline]
    pub fn total(&self) -> u32 {
        self.nref + self.nalt
    }

    /// Allele frequency of the alternate allele.
    /// Returns `None` if no data.
    #[inline]
    pub fn freq(&self) -> Option<f64> {
        let t = self.total();
        if t == 0 {
            None
        } else {
            Some(self.nalt as f64 / t as f64)
        }
    }
}

/// Genotype file format enumeration.  Mirrors `outputmodetype` in `mcio.h`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenoFormat {
    AncestryMap,
    Eigenstrat,
    Ped,
    PackedPed,
    PackedAncestryMap,
    TransposePacked,
}

/// A single jackknife block.
#[derive(Debug, Clone)]
pub struct Block {
    /// Block index (0-based).
    pub index: usize,
    /// SNP column indices belonging to this block.
    pub snp_indices: Vec<usize>,
}

/// F-statistic basis entry identifier.
///
/// Each entry in the canonical basis is either an f2, f3, or f4,
/// indexed by population indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FStatKind {
    F2(usize, usize),
    F3(usize, usize, usize),
    F4(usize, usize, usize, usize),
}

/// Packed result of the `ranktest` / `ranktestfix` family.
/// Mirrors `F4INFO` in `f4rank.h`.
#[derive(Debug, Clone)]
pub struct F4Info {
    pub nl: usize,
    pub nr: usize,
    pub rank: usize,
    /// Degrees of freedom from jackknife.
    pub dof_jack: f64,
    /// Degrees of freedom from approximate chi-square.
    pub dof: f64,
    pub dof_diff: f64,
    pub chisq: f64,
    pub chisq_diff: f64,
    /// Left factor matrix (nl × rank), column-major.
    pub a: Vec<f64>,
    /// Right factor matrix (nr × rank), column-major.
    pub b: Vec<f64>,
    /// Mean matrix (nl × nr), column-major.
    pub mean: Vec<f64>,
    /// Residual matrix (nl × nr), column-major.
    pub resid: Vec<f64>,
}
