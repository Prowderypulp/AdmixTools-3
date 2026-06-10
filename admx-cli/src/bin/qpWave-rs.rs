use std::env;
use std::path::Path;
use admx_core::error::AdmxResult;
use admx_io::params::ParFile;
use admx_qpwave::driver::{run_qpwave, QpWaveConfig};
use admx_rank::stats::rtlchsq;
use std::fs::read_to_string;

fn printstrmat(ss: &[String], mat: &[f64], m: usize, n: usize) {
    if m == 0 || n == 0 { return; }
    
    let mut scale = vec![0.0; n];
    for col in 0..n {
        let mut sum_sq = 0.0;
        for row in 0..m {
            let val = mat[row * n + col];
            sum_sq += val * val;
        }
        if sum_sq >= 0.0 {
            let y = sum_sq.sqrt() / (m as f64).sqrt();
            scale[col] = 1.0 / y;
        }
    }
    
    print!("{:>15} ", "scale");
    for s in &scale {
        print!("{:9.3} ", s);
    }
    println!();
    
    for row in 0..m {
        print!("{:>15} ", ss[row]);
        for col in 0..n {
            print!("{:9.3} ", mat[row * n + col] * scale[col]);
        }
        println!();
    }
}

fn main() -> AdmxResult<()> {
    let args: Vec<String> = env::args().collect();
    let par_idx = args.iter().position(|a| a == "-p").expect("Missing -p parfile");
    let par_path = Path::new(&args[par_idx + 1]);
    let pf = ParFile::load(par_path)?;

    println!("qpWave: parameter file: {}", par_path.display());
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
    let _details = pf.get_bool("details:", false);

    let genotypename = pf.get_string("genotypename:");
    let snpname = pf.get_string("snpname:");
    let indivname = pf.get_string("indivname:");
    let numchrom = pf.get_int("numchrom:").unwrap_or(22);
    let noxdata = pf.get_int("noxdata:").unwrap_or(1) != 0;
    let blgsize = pf.get_dbl("blgsize:").unwrap_or(0.05);

    let config = QpWaveConfig {
        fstatsname,
        popleft: popleft.clone(),
        popright: popright.clone(),
        allsnps,
        fancyf4,
        yscale: 0.0001,
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

    let results = run_qpwave(&config)?;

    for res in &results {
        println!("f4info: ");
        print!("f4rank: {} dof: {:6.0} chisq: {:9.3} tail: {:20.9e} ", 
            res.rank, res.dof, res.chisq, rtlchsq(res.dof as usize, res.chisq));
        println!("dofdiff: {:6.0} chisqdiff: {:9.3} taildiff: {:20.9e}", 
            res.dof_diff, res.chisq_diff, rtlchsq(res.dof_diff as usize, res.chisq_diff));
        
        if res.rank > 0 {
            println!("B:");
            // B is nr x rank
            // Let's create BT = B^T
            let mut bt = vec![0.0; res.nr * res.rank];
            for r in 0..res.rank {
                for c in 0..res.nr {
                    bt[c * res.rank + r] = res.b[r * res.nr + c];
                }
            }
            printstrmat(&popright[1..], &bt, res.nr, res.rank);
            
            println!("A:");
            printstrmat(&popleft[1..], &res.a, res.nl, res.rank);
        }
        println!();
    }

    Ok(())
}

