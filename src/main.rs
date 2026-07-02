use nbody_sim::analysis;
use nbody_sim::body::State;
use nbody_sim::integrator::{self, Method};
use nbody_sim::io::{self, CsvWriter};
use nbody_sim::plot;
use std::path::PathBuf;
use std::time::Instant;

struct Args {
    input: PathBuf,
    output: PathBuf,
    dt: f64,
    years: f64,
    softening: f64,
    method: Method,
    snapshot_every: usize,
    plot: bool,
}

impl Args {
    fn parse() -> Result<Args, String> {
        let mut args = Args {
            input: PathBuf::from("data/solar_system.json"),
            output: PathBuf::from("output/sim.csv"),
            dt: 0.001,
            years: 100.0,
            softening: 1e-4,
            method: Method::Rk4,
            snapshot_every: 100,
            plot: false,
        };
        let mut it = std::env::args().skip(1);
        while let Some(flag) = it.next() {
            let mut value = |name: &str| {
                it.next().ok_or_else(|| format!("{name} needs a value"))
            };
            match flag.as_str() {
                "--input" => args.input = PathBuf::from(value("--input")?),
                "--output" => args.output = PathBuf::from(value("--output")?),
                "--dt" => args.dt = value("--dt")?.parse().map_err(|e| format!("--dt: {e}"))?,
                "--years" => {
                    args.years = value("--years")?.parse().map_err(|e| format!("--years: {e}"))?
                }
                "--softening" => {
                    args.softening = value("--softening")?
                        .parse()
                        .map_err(|e| format!("--softening: {e}"))?
                }
                "--integrator" => args.method = value("--integrator")?.parse()?,
                "--snapshot-every" => {
                    args.snapshot_every = value("--snapshot-every")?
                        .parse()
                        .map_err(|e| format!("--snapshot-every: {e}"))?
                }
                "--plot" => args.plot = true,
                "--help" | "-h" => {
                    println!("{USAGE}");
                    std::process::exit(0);
                }
                other => return Err(format!("unknown argument '{other}'\n{USAGE}")),
            }
        }
        if args.dt <= 0.0 || args.years <= 0.0 {
            return Err("--dt and --years must be positive".into());
        }
        Ok(args)
    }
}

const USAGE: &str = "\
nbody-sim — N-body gravitational simulation (units: AU, yr, M☉)

USAGE: nbody-sim [OPTIONS]

OPTIONS:
    --input <FILE>          initial conditions JSON [default: data/solar_system.json]
    --output <FILE>         trajectory CSV [default: output/sim.csv]
    --dt <YEARS>            time step [default: 0.001]
    --years <YEARS>         total duration [default: 100]
    --softening <AU>        softening length ε [default: 1e-4]
    --integrator <NAME>     rk4 | euler [default: rk4]
    --snapshot-every <N>    write every N-th step [default: 100]
    --plot                  write orbit + energy-drift SVGs next to the output CSV";

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse()?;
    let bodies = io::load_bodies(&args.input)?;
    let masses: Vec<f64> = bodies.iter().map(|b| b.mass).collect();

    let mut state = State::from_bodies(&bodies);
    state.to_com_frame(&masses);

    let steps = (args.years / args.dt).round() as usize;
    let e0 = analysis::conserved(&state, &masses, args.softening).total_energy;

    println!(
        "{} bodies | {} | dt = {} yr | {} yr = {} steps | ε = {} AU",
        bodies.len(),
        match args.method {
            Method::Rk4 => "RK4",
            Method::Euler => "Euler",
        },
        args.dt,
        args.years,
        steps,
        args.softening
    );

    let mut writer = CsvWriter::create(&args.output)?;
    let mut trajectories: Vec<(String, Vec<(f64, f64)>)> = bodies
        .iter()
        .map(|b| (b.name.clone(), Vec::new()))
        .collect();
    let mut drift_series: Vec<(f64, f64)> = Vec::new();

    let started = Instant::now();
    for step in 0..=steps {
        let t = step as f64 * args.dt;

        if step % args.snapshot_every == 0 || step == steps {
            let c = analysis::conserved(&state, &masses, args.softening);
            let drift = analysis::relative_drift(c.total_energy, e0);
            writer.snapshot(t, &bodies, &state, &c, drift)?;
            for (i, (_, pts)) in trajectories.iter_mut().enumerate() {
                pts.push((state.positions[i].x, state.positions[i].y));
            }
            drift_series.push((t, drift));
        }

        if step < steps {
            state = integrator::step(args.method, &state, &masses, args.dt, args.softening);
        }
    }
    let elapsed = started.elapsed();
    writer.finish()?;

    let c = analysis::conserved(&state, &masses, args.softening);
    let final_drift = analysis::relative_drift(c.total_energy, e0);
    println!(
        "done in {:.2} s | final |ΔE/E0| = {:.3e} | |L| = {:.6e}",
        elapsed.as_secs_f64(),
        final_drift,
        c.angular_momentum.norm()
    );
    println!("wrote {}", args.output.display());

    if args.plot {
        let dir = args.output.parent().unwrap_or_else(|| ".".as_ref());
        let orbit_path = dir.join("orbit_plot.svg");
        let drift_path = dir.join("energy_drift.svg");
        plot::orbit_plot(&trajectories, &orbit_path, "Orbit trajectories (xy plane)")?;
        plot::drift_plot(
            &[("this run".to_string(), drift_series)],
            &drift_path,
            "Relative energy drift",
        )?;
        println!("wrote {} and {}", orbit_path.display(), drift_path.display());
    }
    Ok(())
}
