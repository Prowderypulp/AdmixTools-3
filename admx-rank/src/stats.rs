use core::f64::consts::PI;

const LOG_SQRT_PI: f64 = 0.5723649429247000870717135; // log(sqrt(PI))
const I_SQRT_PI: f64 = 0.5641895835477562869480795;   // 1 / sqrt(PI)
const BIGX: f64 = 20.0;

fn ntail(zval: f64) -> f64 {
    if zval == 0.0 {
        return 0.5;
    }
    if zval < 0.0 {
        return 1.0 - ntail(-zval);
    }
    if zval < 20.0 {
        let t = zval / std::f64::consts::SQRT_2;
        return libm::erfc(t) / 2.0;
    }
    
    let t = (-0.5 * zval * zval).exp();
    t / ((2.0 * PI).sqrt() * zval)
}

fn pochisq(x: f64, df: usize) -> f64 {
    if x <= 0.0 || df < 1 {
        return 1.0;
    }

    let a = 0.5 * x;
    let even = (df % 2) == 0;
    let mut y = 0.0;
    if df > 1 {
        y = (-a).exp();
    }
    let mut s = if even { y } else { 2.0 * ntail(x.sqrt()) };

    if df > 2 {
        let x_lim = 0.5 * (df as f64 - 1.0);
        let mut z = if even { 1.0 } else { 0.5 };
        
        if a > BIGX {
            let mut e = if even { 0.0 } else { LOG_SQRT_PI };
            let c = a.ln();
            while z <= x_lim {
                e = z.ln() + e;
                s += (c * z - a - e).exp();
                z += 1.0;
            }
            return s;
        } else {
            let mut e = if even { 1.0 } else { I_SQRT_PI / a.sqrt() };
            let mut c = 0.0;
            while z <= x_lim {
                e = e * (a / z);
                c = c + e;
                z += 1.0;
            }
            return c * y + s;
        }
    }
    
    s
}

fn ltlg1(a: f64, x: f64) -> f64 {
    let tiny = 1.0e-14;
    let mut s = 1.0 / a;
    let mut r = s;
    for k in 1..=60 {
        let yk = k as f64;
        r *= x / (a + yk);
        s += r;
        if (r / s).abs() < tiny {
            break;
        }
    }
    let xam = (a * x.ln()) - x;
    let mut y1 = xam + s.ln();
    let (lgamma_a, _) = libm::lgamma_r(a);
    y1 -= lgamma_a;
    y1.exp()
}

fn rtlg2(a: f64, x: f64) -> f64 {
    let mut t0 = 0.0;
    for k in (1..=60).rev() {
        let yk = k as f64;
        let top = yk - a;
        let mut bot = yk / (x + t0);
        bot += 1.0;
        t0 = top / bot;
    }
    let xam = (a * x.ln()) - x;
    let mut y1 = xam - (x + t0).ln();
    let (lgamma_a, _) = libm::lgamma_r(a);
    y1 -= lgamma_a;
    y1.exp()
}

fn rtlg(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    if x <= 1.0 + a {
        return 1.0 - ltlg1(a, x);
    }
    rtlg2(a, x)
}

pub fn rtlchsq(df: usize, z: f64) -> f64 {
    if df == 1 {
        return 2.0 * ntail(z.sqrt());
    }
    if df == 2 {
        return (-0.5 * z).exp();
    }
    if df == 0 {
        return 1.0;
    }
    let y = pochisq(z, df);
    if y < 1.0e-6 {
        let a = 0.5 * (df as f64);
        let x = 0.5 * z;
        return rtlg(a, x);
    }
    y
}
