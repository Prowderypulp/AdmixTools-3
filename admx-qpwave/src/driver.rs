use admx_core::error::{AdmxError, AdmxResult};
use admx_core::types::{F4Info, Snp, Indiv};
use admx_fstats::fstats_io::load_fstats;
use admx_rank::ranktest::doranktest;
use admx_rank::checkmv::checkmv;
use ndarray::{Array1, Array2};
use std::path::Path;

use admx_fstats::driver::{run_qpfstats, QpfstatsConfig};
use admx_fstats::basis::FBasis;

pub struct QpWaveConfig {
    pub fstatsname: Option<String>,
    pub popleft: Vec<String>,
    pub popright: Vec<String>,
    pub allsnps: bool,
    pub fancyf4: bool,
    pub yscale: f64,
    pub genotypename: Option<String>,
    pub snpname: Option<String>,
    pub indivname: Option<String>,
    pub numchrom: i32,
    pub noxdata: bool,
    pub blgsize: f64,
}

impl Default for QpWaveConfig {
    fn default() -> Self {
        Self {
            fstatsname: None,
            popleft: vec![],
            popright: vec![],
            allsnps: false,
            fancyf4: true,
            yscale: 0.0001,
            genotypename: None,
            snpname: None,
            indivname: None,
            numchrom: 22,
            noxdata: true,
            blgsize: 0.05,
        }
    }
}

pub fn run_qpwave(config: &QpWaveConfig) -> AdmxResult<Vec<F4Info>> {
    let (means, covar, pop_labels, basis_indices, anchor) = if let Some(ref fstats_file) = config.fstatsname {
        load_fstats(Path::new(fstats_file))?
    } else {
        // Direct-genotype path
        if !config.fancyf4 {
            return Err(AdmxError::Fatal("fancyf4: NO is not supported in the direct-genotype path. Pre-compute fstats or use fancyf4: YES.".into()));
        }
        if config.genotypename.is_none() || config.snpname.is_none() || config.indivname.is_none() {
             return Err(AdmxError::Fatal("Must provide either fstatsname or (genotypename, snpname, indivname)".into()));
        }

        let genotypename = config.genotypename.as_ref().unwrap();
        let snpname = config.snpname.as_ref().unwrap();
        let indivname = config.indivname.as_ref().unwrap();

        let mut pop_list = Vec::new();
        pop_list.push(config.popleft[0].clone());
        for p in &config.popleft[1..] {
            if !pop_list.contains(p) { pop_list.push(p.clone()); }
        }
        for p in &config.popright {
            if !pop_list.contains(p) { pop_list.push(p.clone()); }
        }

        let snp_rows = admx_io::snp::read(Path::new(snpname), config.numchrom as u32).map_err(AdmxError::Io)?;
        let ind_rows = admx_io::indiv::read(Path::new(indivname)).map_err(AdmxError::Io)?;

        let mut snps: Vec<Snp> = snp_rows.iter().map(|r| Snp {
            id: r.id.clone(),
            chrom: r.chrom as i32,
            cchrom: String::new(), 
            genpos: r.genetic_pos,
            physpos: r.physical_pos as f64,
            alleles: [r.allele1 as char, r.allele2 as char],
            ignore: false,
            tagnumber: -1,
            weight: 1.0,
        }).collect();

        for snp in snps.iter_mut() {
            if config.noxdata && snp.chrom == (config.numchrom + 1) { snp.ignore = true; }
            if snp.chrom > (config.numchrom + 1) { snp.ignore = true; }
            if snp.chrom == 0 { snp.ignore = true; }
        }

        let indivs: Vec<Indiv> = ind_rows.iter().map(|r| Indiv {
            id: r.id.clone(),
            egroup: r.pop.clone(),
            sex: r.sex,
            idnum: 0,
            affstatus: 0,
            ignore: r.ignore,
            gkode: 0,
        }).collect();

        let mut reader: Box<dyn admx_io::GenoReader> = if admx_io::is_packed_am(&genotypename) {
            Box::new(admx_io::packed_am::PackedAmReader::open(Path::new(&genotypename), indivs.len(), snps.len())
                .map_err(AdmxError::Io)?)
        } else if admx_io::is_bed(&genotypename) {
            Box::new(admx_io::packed_ped::PackedPedReader::open(Path::new(&genotypename), indivs.len(), snps.len())
                .map_err(AdmxError::Io)?)
        } else if admx_io::is_eigenstrat(&genotypename, indivs.len()) {
            Box::new(admx_io::eigenstrat::EigenstratReader::open(Path::new(&genotypename), indivs.len(), snps.len())
                .map_err(AdmxError::Io)?)
        } else {
            return Err(AdmxError::Fatal("Unsupported genotype format".to_string()));
        };

        let qpf_config = QpfstatsConfig {
            blgsize: config.blgsize,
            inbreed: false,
            hires: true,
            allsnps: config.allsnps,
            noxdata: config.noxdata,
            numchrom: config.numchrom,
            doscale: true,
            anchor_pop: pop_list[0].clone(),
        };

        let result = run_qpfstats(reader.as_mut(), &snps, &indivs, &pop_list, &qpf_config)?;
        
        let basis = FBasis::new(0, pop_list.len()); 
        let basis_indices: Vec<(usize, usize)> = basis.pops.iter()
            .flat_map(|&a| basis.pops.iter().filter(move |&&b| b >= a).map(move |&b| (a, b)))
            .collect();

        (result.means, result.covar, pop_list.clone(), basis_indices, pop_list[0].clone())
    };

    let nl = config.popleft.len() - 1;
    let nr = config.popright.len() - 1;
    
    let mut ymean = vec![0.0; nl * nr];
    let mut yvar = vec![0.0; (nl * nr) * (nl * nr)];
    
    let get_pop_idx = |name: &str| -> AdmxResult<usize> {
        pop_labels.iter().position(|p| p == name)
            .ok_or_else(|| AdmxError::Fatal(format!("Population {} not found in fstats", name)))
    };

    let base_idx = get_pop_idx(&anchor)?;
    
    let get_f3_idx = |x: usize, y: usize| -> Option<usize> {
        if x == base_idx || y == base_idx { return None; }
        let (a, b) = if x < y { (x, y) } else { (y, x) };
        basis_indices.iter().position(|&(i, j)| i == a && j == b)
    };

    let mut fsindex = vec![[-1isize; 4]; nl * nr];
    let a = get_pop_idx(&config.popleft[0])?;
    let c = get_pop_idx(&config.popright[0])?;
    
    for i in 1..=nl {
        let b = get_pop_idx(&config.popleft[i])?;
        for j in 1..=nr {
            let d = get_pop_idx(&config.popright[j])?;
            
            let idx = (i - 1) * nr + (j - 1);
            let c1 = get_f3_idx(a, c);
            let c2 = get_f3_idx(b, d);
            let c3 = get_f3_idx(a, d);
            let c4 = get_f3_idx(b, c);
            
            fsindex[idx][0] = c1.map(|x| x as isize).unwrap_or(-1);
            fsindex[idx][1] = c2.map(|x| x as isize).unwrap_or(-1);
            fsindex[idx][2] = c3.map(|x| x as isize).unwrap_or(-1);
            fsindex[idx][3] = c4.map(|x| x as isize).unwrap_or(-1);
            
            let mut mean = 0.0;
            if let Some(x) = c1 { mean += means[x]; }
            if let Some(x) = c2 { mean += means[x]; }
            if let Some(x) = c3 { mean -= means[x]; }
            if let Some(x) = c4 { mean -= means[x]; }
            ymean[idx] = mean;
        }
    }
    
    for i1 in 0..nl * nr {
        for i2 in 0..nl * nr {
            let fs1 = fsindex[i1];
            let fs2 = fsindex[i2];
            let mut v = 0.0;
            
            for k1 in 0..4 {
                if fs1[k1] < 0 { continue; }
                let sign1 = if k1 < 2 { 1.0 } else { -1.0 };
                for k2 in 0..4 {
                    if fs2[k2] < 0 { continue; }
                    let sign2 = if k2 < 2 { 1.0 } else { -1.0 };
                    v += sign1 * sign2 * covar[[fs1[k1] as usize, fs2[k2] as usize]];
                }
            }
            yvar[i1 * (nl * nr) + i2] = v;
        }
    }
    let diagvarplus = config.yscale;
    let mut trace = 0.0;
    for i in 0..(nl * nr) {
        trace += yvar[i * (nl * nr) + i];
    }
    let y = diagvarplus * trace;
    for i in 0..(nl * nr) {
        yvar[i * (nl * nr) + i] += y;
    }

    let ret = checkmv(&ymean, &yvar, nl, nr);

    if ret == -1 {
        return Err(AdmxError::Fatal("f4 stats all zero. Rank 0!. Aborting run".into()));
    }
    if ret == -2 {
        return Err(AdmxError::Fatal("f4 variance absursly small. Aborting run".into()));
    }

    let maxrank = std::cmp::min(nl, nr);
    let mut results: Vec<F4Info> = Vec::new();
    
    for x in 0..=maxrank {
        let mut f4pt = F4Info {
            nl, nr, rank: x, dof_jack: 0.0, dof: 0.0, dof_diff: 0.0,
            chisq: 0.0, chisq_diff: 0.0, a: vec![0.0; nl * x], b: vec![0.0; nr * x],
            mean: vec![0.0; nl * nr], resid: vec![0.0; nl * nr],
        };
        
        doranktest(&ymean, &yvar, nl, nr, x, config.yscale, &mut f4pt);
        
        if x > 0 {
            let prev = &results[x - 1];
            f4pt.dof_diff = prev.dof - f4pt.dof;
            f4pt.chisq_diff = prev.chisq - f4pt.chisq;
        }
        results.push(f4pt);
    }
    
    Ok(results)
}
