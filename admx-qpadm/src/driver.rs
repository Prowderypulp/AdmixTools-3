use std::path::Path;
use admx_core::error::{AdmxError, AdmxResult};
use admx_core::types::{F4Info, Snp, Indiv};
use admx_fstats::basis::FBasis;
use admx_fstats::fstats_io::load_fstats;
use admx_fstats::driver::{run_qpfstats, QpfstatsConfig};
use admx_rank::ranktest::{doranktest, doranktestfix};
use admx_rank::checkmv::checkmv;
use admx_rank::stats::rtlchsq;
use crate::calcadm::{calcadm, calcadmfix};
use crate::bootstrap::calcevarboot;
use crate::prng::LegacyLcg;

pub struct QpAdmConfig {
    pub fstatsname: Option<String>,
    pub popleft: Vec<String>,
    pub popright: Vec<String>,
    pub allsnps: bool,
    pub fancyf4: bool,
    pub inbreed: bool,
    pub yscale: f64,
    pub numboot: usize,
    pub seed: u32,
    pub genotypename: Option<String>,
    pub snpname: Option<String>,
    pub indivname: Option<String>,
    pub numchrom: i32,
    pub noxdata: bool,
    pub blgsize: f64,
}

fn binary_string(mut val: usize, len: usize) -> String {
    let mut s = String::new();
    for _ in 0..len {
        if val % 2 == 1 { s.push('1'); } else { s.push('0'); }
        val /= 2;
    }
    s
}

pub fn run_qpadm(config: &QpAdmConfig) -> AdmxResult<()> {
    let (means, covar, pop_labels, basis_indices, anchor): (ndarray::Array1<f64>, ndarray::Array2<f64>, Vec<String>, Vec<(usize, usize)>, String) = if let Some(ref fstats_file) = config.fstatsname {
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
            inbreed: config.inbreed,
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

    let ret = checkmv(&ymean, &yvar, nl, nr);
    if ret == -1 {
        return Err(AdmxError::Fatal("f4 stats all zero. Rank 0!. Aborting run".into()));
    }
    if ret == -2 {
        return Err(AdmxError::Fatal("f4 variance absursly small. Aborting run".into()));
    }

    // C qpAdm.c:266-268,597-598 — for qpAdm the default is diagvarplus = 0.0
    // (NOT yscale, as in qpWave), so addscaldiag(yvar, diagvarplus) adds
    // nothing to yvar here. The yscale regularization is applied only to
    // varinv inside ranktest. (checkmv above runs on the raw yvar, matching C.)
    let diagvarplus = 0.0_f64;
    if diagvarplus != 0.0 {
        let d = nl * nr;
        let trace: f64 = (0..d).map(|i| yvar[i * d + i]).sum();
        let yadd = diagvarplus * trace;
        for i in 0..d {
            yvar[i * d + i] += yadd;
        }
    }

    let mut results: Vec<F4Info> = Vec::new();
    let maxrank = nl;
    let rank = nl - 1;

    for x in rank..=maxrank {
        let mut f4pt = F4Info {
            nl, nr, rank: x, dof_jack: 0.0, dof: 0.0, dof_diff: 0.0,
            chisq: 0.0, chisq_diff: 0.0, a: vec![0.0; nl * x], b: vec![0.0; nr * x],
            mean: vec![0.0; nl * nr], resid: vec![0.0; nl * nr],
        };
        doranktest(&ymean, &yvar, nl, nr, x, config.yscale, &mut f4pt);
        if x > 0 && !results.is_empty() {
            let prev = &results[0]; // results[0] is rank nl - 1
            f4pt.dof_diff = prev.dof - f4pt.dof;
            f4pt.chisq_diff = prev.chisq - f4pt.chisq;
        }
        results.push(f4pt);
    }

    let mut pval = 0.0;
    for (i, x) in (rank..=maxrank).enumerate() {
        if x == rank { println!("codimension 1"); }
        if x == maxrank { println!("full rank"); }
        let f4pt = &results[i];
        println!("f4info: rank: {} dof: {:.3} chisq: {:.3}", f4pt.rank, f4pt.dof, f4pt.chisq);
        pval = rtlchsq(f4pt.dof_diff.round() as usize, f4pt.chisq_diff);
        println!();
    }

    let f4pt = &results[0];
    let mut ww = vec![0.0; nl];
    calcadm(&mut ww, &f4pt.a, nl).map_err(|_| AdmxError::Fatal("calcadm failed".into()))?;

    println!("best coefficients: ");
    for w in &ww { print!("{:9.3} ", w); }
    println!();

    let mut jmean = vec![0.0; nl];
    let mut var = vec![0.0; nl * nl];
    let mut lcg = LegacyLcg::new(config.seed as u64);
    calcevarboot(&mut jmean, &mut var, &ymean, &yvar, nl, nr, nl * nr, config.numboot, &mut lcg).map_err(|_| AdmxError::Fatal("calcevarboot failed".into()))?;

    println!("bootstrap saimpling of F-coeffs: {}", config.numboot);
    println!("zzjmean ");
    for w in &jmean { print!("{:9.3} ", w); }
    println!();

    println!("      std. errors: ");
    for i in 0..nl { print!("{:9.3} ", var[i * nl + i].sqrt()); }
    println!();
    println!();

    let mut jvar = vec![0; nl * nl];
    for i in 0..nl*nl {
        jvar[i] = (var[i] * 1.0e6).round() as i64;
    }
    println!("error covariance (* 1,000,000)");
    for i in 0..nl {
        for j in 0..nl {
            print!("{:6} ", jvar[i * nl + j]);
        }
        println!();
    }
    println!();
    println!();

    if nl == 1 {
        println!("single source. terminating");
        println!("run qpDstat if you want Z scores for f4 stats");
        println!("## end of run");
        return Ok(());
    }

    let xmax = (1 << nl) - 1;
    let mut nfeas = vec![0; nl];
    let mut yscbest = vec![-1.0e5; nl];
    let mut sbest = vec![0; nl];
    let mut ychi = vec![0.0; nl];

    println!(" {:12}  wt  dof     chisq       tail prob", "fixed pat");
    for wt in 0..nl {
        for k in 0..=xmax {
            let mut vfix = vec![0i32; nl];
            let mut t = 0;
            for i in 0..nl {
                vfix[i] = (k >> i) & 1;
                t += vfix[i];
            }
            if t != wt as i32 { continue; }

            let mut f4pt_fix = F4Info {
                nl, nr, rank: nl - 1, dof_jack: 0.0, dof: 0.0, dof_diff: 0.0,
                chisq: 0.0, chisq_diff: 0.0, a: vec![0.0; nl * (nl - 1)], b: vec![0.0; nr * (nl - 1)],
                mean: vec![0.0; nl * nr], resid: vec![0.0; nl * nr],
            };

            let mut ww_fix = vec![0.0; nl];

            if wt == 0 {
                f4pt_fix = results[0].clone();
                calcadm(&mut ww_fix, &f4pt_fix.a, nl).map_err(|_| AdmxError::Fatal("calcadm failed".into()))?;
            } else {
                doranktestfix(&ymean, &yvar, nl, nr, nl - 1, config.yscale, &mut f4pt_fix, &vfix);
                calcadmfix(&mut ww_fix, &f4pt_fix.a, nl, &vfix).map_err(|_| AdmxError::Fatal("calcadmfix failed".into()))?;
            }

            let mut yfeas = ww_fix[0];
            for w in &ww_fix {
                if *w < yfeas { yfeas = *w; }
            }

            let dof = f4pt_fix.dof.round() as i32;
            let mut tail = 0.0;
            let mut ttail = 0.0;
            if dof > 0 {
                tail = rtlchsq(dof as usize, f4pt_fix.chisq);
                ttail = tail;
            }

            let mut isfeas = false;
            if yfeas >= -0.001 {
                isfeas = true;
                nfeas[wt as usize] += 1;
            }

            if !isfeas { tail = 0.0; }
            if ttail < 1.0e-30 { ttail = 0.0; }

            if nfeas[wt as usize] == 1 {
                yscbest[wt as usize] = tail;
                sbest[wt as usize] = k;
                ychi[wt as usize] = f4pt_fix.chisq;
            }

            if nfeas[wt as usize] == 0 && ttail > yscbest[wt as usize] {
                yscbest[wt as usize] = ttail;
                sbest[wt as usize] = k;
                ychi[wt as usize] = f4pt_fix.chisq;
            }

            if tail > yscbest[wt as usize] {
                yscbest[wt as usize] = tail;
                sbest[wt as usize] = k;
                ychi[wt as usize] = f4pt_fix.chisq;
            }

            print!(" {:12}  {}   {:3} {:9.3} {:15.6} ", binary_string(k as usize, nl), wt, dof, f4pt_fix.chisq, ttail);
            for w in &ww_fix { print!("{:9.3} ", w); }
            if !isfeas { print!(" infeasible"); }
            println!();
        }
    }

    let mut qcoeffs = vec![0.0; nl];
    for w in 0..nl {
        let best_k = sbest[w];
        let mut vfix = vec![0; nl];
        for i in 0..nl { vfix[i] = (best_k >> i) & 1; }
        
        let mut f4pt_fix = F4Info {
            nl, nr, rank: nl - 1, dof_jack: 0.0, dof: 0.0, dof_diff: 0.0,
            chisq: 0.0, chisq_diff: 0.0, a: vec![0.0; nl * (nl - 1)], b: vec![0.0; nr * (nl - 1)],
            mean: vec![0.0; nl * nr], resid: vec![0.0; nl * nr],
        };

        if w == 0 {
            f4pt_fix = results[0].clone();
            calcadm(&mut qcoeffs, &f4pt_fix.a, nl).map_err(|_| AdmxError::Fatal("calcadm failed".into()))?;
        } else {
            doranktestfix(&ymean, &yvar, nl, nr, nl - 1, config.yscale, &mut f4pt_fix, &vfix);
            calcadmfix(&mut qcoeffs, &f4pt_fix.a, nl, &vfix).map_err(|_| AdmxError::Fatal("calcadmfix failed".into()))?;
        }

        let dof_diff = if w > 0 {
            let prev_k = sbest[w - 1];
            let mut vfix_prev = vec![0; nl];
            for i in 0..nl { vfix_prev[i] = (prev_k >> i) & 1; }
            let mut prev_f4pt = F4Info {
                nl, nr, rank: nl - 1, dof_jack: 0.0, dof: 0.0, dof_diff: 0.0,
                chisq: 0.0, chisq_diff: 0.0, a: vec![0.0; nl * (nl - 1)], b: vec![0.0; nr * (nl - 1)],
                mean: vec![0.0; nl * nr], resid: vec![0.0; nl * nr],
            };
            if w - 1 == 0 {
                prev_f4pt = results[0].clone();
            } else {
                doranktestfix(&ymean, &yvar, nl, nr, nl - 1, config.yscale, &mut prev_f4pt, &vfix_prev);
            }
            f4pt_fix.dof - prev_f4pt.dof
        } else {
            0.0
        };

        let chisq_diff = if w > 0 {
            let prev_k = sbest[w - 1];
            let mut vfix_prev = vec![0; nl];
            for i in 0..nl { vfix_prev[i] = (prev_k >> i) & 1; }
            let mut prev_f4pt = F4Info {
                nl, nr, rank: nl - 1, dof_jack: 0.0, dof: 0.0, dof_diff: 0.0,
                chisq: 0.0, chisq_diff: 0.0, a: vec![0.0; nl * (nl - 1)], b: vec![0.0; nr * (nl - 1)],
                mean: vec![0.0; nl * nr], resid: vec![0.0; nl * nr],
            };
            if w - 1 == 0 {
                prev_f4pt = results[0].clone();
            } else {
                doranktestfix(&ymean, &yvar, nl, nr, nl - 1, config.yscale, &mut prev_f4pt, &vfix_prev);
            }
            f4pt_fix.chisq - prev_f4pt.chisq
        } else {
            0.0
        };

        let p_nested = if w > 0 && dof_diff > 0.0 {
            rtlchsq(dof_diff.round() as usize, chisq_diff)
        } else {
            0.0
        };

        println!("best pat: {:12} {:15.6}          -  - ", binary_string(best_k as usize, nl), yscbest[w]);
        if w > 0 {
            println!("best pat: {:12} {:15.6}  p-value for nested model: {:15.6}", binary_string(best_k as usize, nl), yscbest[w], p_nested);
        }
    }

    println!("summ: {} {:3} {:12.6} ", config.popleft[0], nl, pval);
    for q in &qcoeffs { print!("{:9.3} ", q); }
    println!();
    
    Ok(())
}
