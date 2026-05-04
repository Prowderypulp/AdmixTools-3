use std::fs::read_to_string;
use std::path::Path;
use admx_core::error::{AdmxError, AdmxResult};
use admx_core::types::{Snp, Indiv, FStatKind};
use admx_io::params::ParFile;
use admx_fstats::driver::{run_qpfstats, QpfstatsConfig};
use admx_fstats::fstats_io::dump_fstats;
use admx_fstats::basis::FBasis;

fn main() -> AdmxResult<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 || args[1] != "-p" {
        eprintln!("Usage: qpfstats -p <parfile>");
        std::process::exit(1);
    }

    let par = ParFile::load(&args[2])?;

    let genotypename = par.get_string("genotypename:").ok_or_else(|| AdmxError::Fatal("genotypename missing".to_string()))?;
    let snpname = par.get_string("snpname:").ok_or_else(|| AdmxError::Fatal("snpname missing".to_string()))?;
    let indivname = par.get_string("indivname:").ok_or_else(|| AdmxError::Fatal("indivname missing".to_string()))?;
    let outpop = par.get_string("outpop:").unwrap_or_else(|| "NULL".to_string());
    let outputname = par.get_string("fstatsoutname:").or_else(|| par.get_string("output:")).ok_or_else(|| AdmxError::Fatal("output name missing".to_string()))?;
    
    let poplistname = par.get_string("poplistname:").or_else(|| par.get_string("poplist:"));
    let mut pop_list: Vec<String> = if let Some(path) = poplistname {
        read_to_string(path).map_err(AdmxError::Io)?
            .lines()
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        return Err(AdmxError::Fatal("poplist missing".to_string()));
    };

    // Reorder pop_list so outpop is at index 0 (matching C fstats2popl behavior)
    if outpop != "NULL" {
        if let Some(pos) = pop_list.iter().position(|p| p == &outpop) {
            let p = pop_list.remove(pos);
            pop_list.insert(0, p);
        } else {
            return Err(AdmxError::Fatal(format!("outpop population {} not in poplist", outpop)));
        }
    }

    let blgsize = par.get_dbl("blgsize:").unwrap_or(0.05);
    let inbreed = par.get_int("inbreed:").unwrap_or(0) != 0;
    let hires = par.get_int("hires:").unwrap_or(1) != 0;
    let allsnps = par.get_int("allsnps:").unwrap_or(0) != 0;
    let noxdata = par.get_int("noxdata:").unwrap_or(1) != 0;
    let numchrom = par.get_int("numchrom:").unwrap_or(22);
    let doscale = par.get_int("doscale:").or_else(|| par.get_int("scale:")).unwrap_or(1) != 0;

    println!("parameter file: {}", args[2]);
    par.write_pars();

    // Load metadata
    let snp_rows = admx_io::snp::read(Path::new(&snpname), numchrom as u32).map_err(AdmxError::Io)?;
    let ind_rows = admx_io::indiv::read(Path::new(&indivname)).map_err(AdmxError::Io)?;

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
        if noxdata && snp.chrom == (numchrom + 1) { snp.ignore = true; }
        if snp.chrom > (numchrom + 1) { snp.ignore = true; }
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

    // Initialize reader
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

    let config = QpfstatsConfig {
        blgsize,
        inbreed,
        hires,
        allsnps,
        noxdata,
        numchrom,
        doscale,
        anchor_pop: outpop.clone(),
    };

    let result = run_qpfstats(reader.as_mut(), &snps, &indivs, &pop_list, &config)?;
    let (means, covar) = (result.means, result.covar);

    let pop_map: std::collections::HashMap<String, usize> = pop_list.iter()
        .enumerate()
        .map(|(i, name)| (name.clone(), i))
        .collect();

    // Diagnostic prints to match C logs
    let np = pop_list.len();
    let nbasis = np * (np - 1) / 2;
    let n = np * (np + 1) / 2;
    let coefsize = n * nbasis;
    println!("np: {} n: {} nbasis: {} coefsize: {}", np, n, nbasis, coefsize);
    println!("outpop: {} ", outpop);
    println!("valid snps: {}", result.valid_snps);

    println!("                     pop : size   INBREED");
    let mut pop_sizes = vec![0usize; np];
    let mut pop_samples = vec![0usize; np];
    for ind in &indivs {
        if ind.ignore { continue; }
        if let Some(&idx) = pop_map.get(&ind.egroup) {
            pop_sizes[idx] += 1;
            pop_samples[idx] += 1;
        }
    }

    for i in 0..np {
        println!("{:3} {:20} : {:4}    {}", i, pop_list[i], pop_sizes[i], if inbreed { "YES" } else { "NO" });
    }

    println!("before setwt numsnps: {}", result.before_setwt_snps);
    println!("setwt numsnps: {}", result.num_snps);
    println!("number of blocks for moving block jackknife: {}", result.num_blocks);
    let total_indivs = pop_sizes.iter().sum::<usize>();
    println!("snps: {}  indivs: {}", result.num_snps, total_indivs);

    for (i, hetrate, valid_snps) in result.pop_stats {
        println!("pop: {:20}  hetrate: {:9.3} valid snps: {:9} samples: {:4}", 
                 pop_list[i], hetrate, valid_snps, pop_samples[i]);
    }

    if doscale {
        println!("lambdascale: {:9.3}", result.lambdascale);
    }

    println!("statistics multiplied by 1000");
    let short_names: Vec<String> = pop_list.iter().map(|s| s[..1].to_string()).collect();
    
    println!("fst:");
    print!("     ");
    for name in &short_names { print!("{:>5}", name); }
    println!(" ");
    for i in 0..np {
        print!("{:>3}  ", short_names[i]);
        for j in 0..np {
            let val = if i == j { 0.0 } else { result.fst[[i, j]] * 1000.0 };
            print!(" {:4.0}", val);
        }
        println!();
    }
    println!();

    println!("f2:");
    print!("     ");
    for name in &short_names { print!("{:>5}", name); }
    println!(" ");
    for i in 0..np {
        print!("{:>3}  ", short_names[i]);
        for j in 0..np {
            let val = if i == j { 0.0 } else { result.f2[[i, j]] * 1000.0 };
            print!(" {:4.0}", val);
        }
        println!();
    }
    println!();

    println!("adjusted sigs: 0");

    // Basis rows
    let basis = FBasis::new(0, np); 
    let basis_stats = basis.stats();
    for i in 0..nbasis {
        let kind = basis_stats[i];
        if let FStatKind::F3(_, a, b) = kind {
            let i1 = pop_list.iter().position(|p| p == &pop_list[a]).unwrap();
            let i2 = pop_list.iter().position(|p| p == &pop_list[b]).unwrap();
            
            println!("basis: {:3} {:3} {:12.6} {:12.6} :: {:12.6} {:12.6}",
                     i1, i2, result.w3[i], result.wls_stderr[i], means[i], covar[[i, i]].sqrt());
        }
    }

    // Dump basis (anchor is at index 0 after reordering)
    let basis_indices: Vec<(usize, usize)> = basis.pops.iter()
        .flat_map(|&a| basis.pops.iter().filter(move |&&b| b >= a).map(move |&b| (a, b)))
        .collect();

    dump_fstats(&outputname, &means, &covar, &pop_list, &basis_indices, &pop_list[0], hires)?;

    println!("fstats file: {} written", outputname);

    Ok(())
}
