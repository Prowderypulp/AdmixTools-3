//! Block jackknife utility.

use crate::types::Snp;

/// Assigns SNPs to jackknife blocks based on genetic distance.
/// 
/// Starts a new block when:
/// 1. Chromosome changes.
/// 2. Genetic distance from the block start exceeds `blgsize`.
/// 
/// Returns a vector of block indices for each SNP. SNPs with `-1` are ignored.
pub fn assign_blocks(snps: &[Snp], blgsize: f64) -> Vec<i32> {
    let mut block_assignments = vec![-1; snps.len()];
    if snps.is_empty() { return block_assignments; }

    let mut current_block = 0;
    let mut block_start_pos = -1.0;
    let mut current_chrom = -1;

    for (i, snp) in snps.iter().enumerate() {
        if snp.ignore { continue; }

        if snp.chrom != current_chrom || block_start_pos < 0.0 || (snp.genpos - block_start_pos).abs() >= blgsize {
            current_block += 1;
            current_chrom = snp.chrom;
            block_start_pos = snp.genpos;
        }
        block_assignments[i] = current_block;
    }

    block_assignments
}
