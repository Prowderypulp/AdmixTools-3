//! Deterministic synthetic EIGENSTRAT generator.
//!
//! Per-SNP allele-freq drift model: ancestral p ~ Uniform(0.05, 0.95),
//! per-pop p_k = clamp(p + N(0, drift_k), 0.01, 0.99), genotype = Binomial(ploidy, p_k).
//! Reproducible from a single u64 seed.

use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use rand::distributions::{Distribution, Uniform};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha20Rng;

#[derive(Clone, Copy, Debug)]
pub enum Sex { M, F, U }

impl Sex {
    fn as_char(self) -> char { match self { Sex::M => 'M', Sex::F => 'F', Sex::U => 'U' } }
}

#[derive(Clone, Debug)]
pub struct Population {
    pub name: String,
    pub n_indivs: usize,
    /// Per-SNP drift std-dev relative to the ancestral allele freq.
    pub drift: f64,
}

#[derive(Clone, Debug)]
pub struct DatasetSpec {
    pub stem: String,
    pub seed: u64,
    pub n_snps: usize,
    /// Fraction of SNPs placed on chromosome X (rest on autosomes 1..=22, round-robin).
    pub frac_x: f64,
    /// Fraction of genotypes randomly set to missing ('9').
    pub frac_missing: f64,
    pub populations: Vec<Population>,
    /// Sex assignment is round-robin over (M, F) per population, except the last 1/4 marked U.
    pub assign_sex: bool,
}

impl DatasetSpec {
    /// Standard tiers used by the parity + time benches.
    pub fn tier(name: &str, seed: u64) -> Option<Self> {
        let stem = name.to_string();
        let pops = |sizes: &[(&str, usize, f64)]| {
            sizes.iter().map(|(n, k, d)| Population {
                name: (*n).to_string(), n_indivs: *k, drift: *d,
            }).collect()
        };
        Some(match name {
            "tiny"  => DatasetSpec { stem, seed, n_snps:    2_000, frac_x: 0.0,  frac_missing: 0.00,
                populations: pops(&[("AG",4,0.05),("BG",4,0.08),("zz",4,0.10),("qq",4,0.12),("QQ",4,0.10),("LM",4,0.09),("ex",4,0.07)]),
                assign_sex: true },
            "small" => DatasetSpec { stem, seed, n_snps:   20_000, frac_x: 0.05, frac_missing: 0.02,
                populations: pops(&[("AG",10,0.05),("BG",10,0.08),("zz",10,0.10),("qq",10,0.12),("QQ",10,0.10),("LM",10,0.09),("ex",10,0.07)]),
                assign_sex: true },
            "med"   => DatasetSpec { stem, seed, n_snps:  100_000, frac_x: 0.05, frac_missing: 0.03,
                populations: pops(&[
                    ("AG",20,0.04),("BG",20,0.06),("zz",20,0.08),("qq",20,0.10),
                    ("QQ",20,0.08),("LM",20,0.07),("ex",20,0.05),("p8",20,0.06),
                    ("p9",20,0.07),("p10",20,0.05),
                ]),
                assign_sex: true },
            "large" => DatasetSpec { stem, seed, n_snps:  500_000, frac_x: 0.05, frac_missing: 0.03,
                populations: pops(&[
                    ("AG",30,0.04),("BG",30,0.06),("zz",30,0.08),("qq",30,0.10),
                    ("QQ",30,0.08),("LM",30,0.07),("ex",30,0.05),("p8",30,0.06),
                    ("p9",30,0.07),("p10",30,0.05),("p11",30,0.06),("p12",30,0.07),
                ]),
                assign_sex: true },
            _ => return None,
        })
    }

    pub fn n_indivs(&self) -> usize { self.populations.iter().map(|p| p.n_indivs).sum() }
}

/// Write `<out_dir>/{stem}.geno`, `.snp`, `.ind` plus left/right pop lists and a `par.qpfstats`.
pub fn write(spec: &DatasetSpec, out_dir: &Path) -> std::io::Result<DatasetPaths> {
    fs::create_dir_all(out_dir)?;
    let mut rng = ChaCha20Rng::seed_from_u64(spec.seed);

    // ---- individuals
    let mut indivs: Vec<(String, Sex, String)> = Vec::new();
    for (pop_idx, pop) in spec.populations.iter().enumerate() {
        for i in 0..pop.n_indivs {
            let sex = if !spec.assign_sex { Sex::U }
                else if (pop_idx * 7 + i) % 4 == 3 { Sex::U }
                else if i % 2 == 0 { Sex::M } else { Sex::F };
            indivs.push((format!("{}_{:03}", pop.name, i), sex, pop.name.clone()));
        }
    }
    let n_ind = indivs.len();

    // ---- SNPs (chrom round-robin: X assigned to first frac_x of indices, rest spread over 1..=22)
    let n_x = (spec.n_snps as f64 * spec.frac_x).round() as usize;
    let snp_paths = out_dir.join(format!("{}.snp", spec.stem));
    let mut snp_w = BufWriter::new(File::create(&snp_paths)?);
    let geno_path = out_dir.join(format!("{}.geno", spec.stem));
    let mut geno_w = BufWriter::new(File::create(&geno_path)?);

    let unif = Uniform::new(0.05_f64, 0.95_f64);
    let pop_index: Vec<usize> = indivs.iter().enumerate().flat_map(|(i, _)| std::iter::once(i)).collect();
    // Map indiv index -> pop index for quick lookup.
    let indiv_pop: Vec<usize> = {
        let mut v = Vec::with_capacity(n_ind);
        for (pi, p) in spec.populations.iter().enumerate() {
            for _ in 0..p.n_indivs { v.push(pi); }
        }
        v
    };
    let _ = pop_index;

    let mut row = vec![b'0'; n_ind + 1];
    row[n_ind] = b'\n';

    for s in 0..spec.n_snps {
        let chrom: u32 = if s < n_x { 23 } else { (s as u32 % 22) + 1 };
        let bp: u64 = 1_000_000 + (s as u64) * 137;
        let cm: f64 = bp as f64 / 1e8;
        writeln!(snp_w, "rs{:07}\t{}\t{:.6}\t{}\tA\tC", s, chrom, cm, bp)?;

        let p_anc: f64 = unif.sample(&mut rng);
        // Per-pop freq draws.
        let p_pop: Vec<f64> = spec.populations.iter().map(|pop| {
            let z: f64 = sample_normal(&mut rng);
            (p_anc + z * pop.drift).clamp(0.01, 0.99)
        }).collect();

        for (i, (_, sex, _)) in indivs.iter().enumerate() {
            let p = p_pop[indiv_pop[i]];
            let haploid = chrom == 23 && matches!(sex, Sex::M);
            let g: u8 = if rng.gen::<f64>() < spec.frac_missing {
                b'9'
            } else if haploid {
                if rng.gen::<f64>() < p { b'2' } else { b'0' }
            } else {
                let a = (rng.gen::<f64>() < p) as u8;
                let b = (rng.gen::<f64>() < p) as u8;
                b'0' + a + b
            };
            row[i] = g;
        }
        geno_w.write_all(&row)?;
    }
    snp_w.flush()?; geno_w.flush()?;

    // ---- ind file
    let ind_path = out_dir.join(format!("{}.ind", spec.stem));
    let mut ind_w = BufWriter::new(File::create(&ind_path)?);
    for (name, sex, pop) in &indivs {
        writeln!(ind_w, "{}\t{}\t{}", name, sex.as_char(), pop)?;
    }
    ind_w.flush()?;

    // ---- pop lists. By convention: outpop = first pop; left = pops[0..3]; right = pops[3..].
    let pops: Vec<&str> = spec.populations.iter().map(|p| p.name.as_str()).collect();
    let n_left = 3.min(pops.len().saturating_sub(2));
    let popleft: Vec<&str> = pops.iter().take(n_left).copied().collect();
    let popright: Vec<&str> = pops.iter().skip(n_left).copied().collect();

    let poplist = out_dir.join("poplist");
    fs::write(&poplist, pops.join("\n") + "\n")?;
    let popleft_p = out_dir.join("popleft");
    fs::write(&popleft_p, popleft.join("\n") + "\n")?;
    let popright_p = out_dir.join("popright");
    fs::write(&popright_p, popright.join("\n") + "\n")?;

    // ---- parfiles
    let par_fstats = out_dir.join("par.qpfstats");
    fs::write(&par_fstats, format!(
        "genotypename: {stem}.geno\nsnpname:      {stem}.snp\nindivname:    {stem}.ind\n\
         poplistname:  poplist\nfstatsoutname: {stem}.fstats\n\
         allsnps:      YES\nhires:        YES\ninbreed:      NO\nnoxdata:      NO\n",
        stem = spec.stem,
    ))?;
    let par_qpwave = out_dir.join("par.qpwave");
    fs::write(&par_qpwave, format!(
        "fstatsname: {stem}.fstats\npopleft:    popleft\npopright:   popright\nallsnps:    YES\n",
        stem = spec.stem,
    ))?;
    let par_qpadm = out_dir.join("par.qpadm");
    fs::write(&par_qpadm, format!(
        "fstatsname: {stem}.fstats\npopleft:    popleft\npopright:   popright\nallsnps:    YES\nseed:       314159265\n",
        stem = spec.stem,
    ))?;

    Ok(DatasetPaths { dir: out_dir.to_path_buf(), stem: spec.stem.clone() })
}

#[derive(Debug, Clone)]
pub struct DatasetPaths {
    pub dir: PathBuf,
    pub stem: String,
}

/// Box–Muller via the underlying RNG. Avoids pulling in rand_distr.
fn sample_normal<R: Rng>(rng: &mut R) -> f64 {
    let u1: f64 = rng.gen::<f64>().max(1e-300);
    let u2: f64 = rng.gen::<f64>();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}
