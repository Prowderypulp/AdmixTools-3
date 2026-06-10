use std::env;
use std::path::Path;
use std::fs::read_to_string;
use admx_core::error::AdmxResult;
use admx_io::params::ParFile;
use admx_qpadm::driver::{run_qpadm, QpAdmConfig};

fn main() -> AdmxResult<()> {
    let args: Vec<String> = env::args().collect();
    let par_idx = args.iter().position(|a| a == "-p").expect("Missing -p parfile");
    let par_path = Path::new(&args[par_idx + 1]);
    let pf = ParFile::load(par_path)?;

    println!("qpAdm: parameter file: {}", par_path.display());
    println!("### THE INPUT PARAMETERS");
    println!("##PARAMETER NAME: VALUE");
    
    let content = read_to_string(par_path).unwrap();
    for line in content.lines() {
        if !line.trim().is_empty() && !line.trim().starts_with('#') {
            let parts: Vec<&str> = line.split(':').collect();
            if parts.len() == 2 {
                println!("{}: {}", parts[0].trim(), parts[1].trim());
            }
        }
    }
    println!("## -------------");

    let fstatsname = pf.get_string("fstatsname:");
    let popleft = pf.get_list("popleft:").unwrap_or_default();
    let popright = pf.get_list("popright:").unwrap_or_default();
    let allsnps = pf.get_bool("allsnps:", false);
    let fancyf4 = pf.get_bool("fancyf4:", true);
    let inbreed = pf.get_bool("inbreed:", false);
    let numboot = pf.get_int("numboot:").unwrap_or(1000) as usize;
    let seed = pf.get_int("seed:").unwrap_or(0) as u32;

    let genotypename = pf.get_string("genotypename:");
    let snpname = pf.get_string("snpname:");
    let indivname = pf.get_string("indivname:");
    let numchrom = pf.get_int("numchrom:").unwrap_or(22);
    let noxdata = pf.get_int("noxdata:").unwrap_or(1) != 0;
    let blgsize = pf.get_dbl("blgsize:").unwrap_or(0.05);

    let config = QpAdmConfig {
        fstatsname,
        popleft: popleft.clone(),
        popright: popright.clone(),
        allsnps,
        fancyf4,
        inbreed,
        yscale: 0.0001,
        numboot,
        seed,
        genotypename,
        snpname,
        indivname,
        numchrom,
        noxdata,
        blgsize,
    };
    
    println!("left pops:");
    for p in &popleft {
        println!("{:>20}", p);
    }
    println!("\nright pops:");
    for p in &popright {
        println!("{:>20}", p);
    }
    println!();

    run_qpadm(&config)?;

    Ok(())
}
