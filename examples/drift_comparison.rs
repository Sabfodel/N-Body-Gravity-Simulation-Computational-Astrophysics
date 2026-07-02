//! RK4 vs Euler energy drift on the full solar system — generates the
//! comparison figure shown in the README.
//!
//! Run: cargo run --release --example drift_comparison

use nbody_sim::analysis;
use nbody_sim::body::State;
use nbody_sim::integrator::{self, Method};
use nbody_sim::io;
use nbody_sim::plot;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bodies = io::load_bodies(Path::new("data/solar_system.json"))?;
    let masses: Vec<f64> = bodies.iter().map(|b| b.mass).collect();
    let mut initial = State::from_bodies(&bodies);
    initial.to_com_frame(&masses);

    let dt = 0.001;
    let years = 100.0;
    let softening = 1e-4;
    let steps = (years / dt) as usize;
    let sample_every = 100;

    let mut series = Vec::new();
    for (label, method) in [("RK4", Method::Rk4), ("Euler", Method::Euler)] {
        let mut state = initial.clone();
        let e0 = analysis::conserved(&state, &masses, softening).total_energy;
        let mut drift = Vec::new();
        for step in 0..=steps {
            if step % sample_every == 0 {
                let e = analysis::conserved(&state, &masses, softening).total_energy;
                drift.push((step as f64 * dt, analysis::relative_drift(e, e0)));
            }
            if step < steps {
                state = integrator::step(method, &state, &masses, dt, softening);
            }
        }
        let final_drift = drift.last().unwrap().1;
        println!("{label}: final |ΔE/E0| = {final_drift:.3e}");
        series.push((label.to_string(), drift));
    }

    let out = Path::new("assets/energy_drift.svg");
    plot::drift_plot(
        &series,
        out,
        "Energy drift — RK4 vs Euler, solar system, dt = 0.001 yr",
    )?;
    println!("wrote {}", out.display());
    Ok(())
}
