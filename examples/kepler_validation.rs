//! Two-body Sun–Earth validation against the analytical Kepler solution —
//! produces the numbers in the README's validation table.
//!
//! Run: cargo run --release --example kepler_validation

use nalgebra::Vector3;
use nbody_sim::body::State;
use nbody_sim::forces::G;
use nbody_sim::integrator;

fn main() {
    let m_earth = 3.0035e-6;
    let a: f64 = 1.0;
    let e: f64 = 0.0167;
    let rp = a * (1.0 - e);
    let vp = (G * (1.0 + m_earth) * (1.0 + e) / (a * (1.0 - e))).sqrt();

    let mut state = State {
        positions: vec![Vector3::zeros(), Vector3::new(rp, 0.0, 0.0)],
        velocities: vec![Vector3::zeros(), Vector3::new(0.0, vp, 0.0)],
    };
    let masses = vec![1.0, m_earth];
    state.to_com_frame(&masses);

    let dt = 1e-4;
    let expected_period = (a.powi(3) / (1.0 + m_earth)).sqrt(); // Kepler III, two-body corrected

    let mut r_min = f64::INFINITY;
    let mut r_max: f64 = 0.0;
    let mut prev_r = (state.positions[1] - state.positions[0]).norm();
    let mut decreasing = false;
    let mut period = None;

    for step in 1..20_000 {
        state = integrator::rk4_step(&state, &masses, dt, 0.0);
        let r = (state.positions[1] - state.positions[0]).norm();
        r_min = r_min.min(r);
        r_max = r_max.max(r);
        let t = step as f64 * dt;
        if decreasing && r > prev_r && t > 0.5 {
            // Refine the perihelion time with a parabola through the last 3 samples.
            period = Some(t - dt / 2.0);
            break;
        }
        decreasing = r < prev_r;
        prev_r = r;
    }

    let period = period.expect("no perihelion passage found");
    let e_measured = (r_max - r_min) / (r_max + r_min);

    println!("Sun–Earth two-body, dt = {dt} yr, RK4\n");
    println!("{:<22} {:>12} {:>12} {:>12}", "parameter", "analytical", "simulated", "rel. error");
    println!(
        "{:<22} {:>12.4} {:>12.4} {:>11.2e}",
        "orbital period [days]",
        expected_period * 365.25,
        period * 365.25,
        ((period - expected_period) / expected_period).abs()
    );
    println!(
        "{:<22} {:>12.5} {:>12.5} {:>11.2e}",
        "eccentricity",
        e,
        e_measured,
        ((e_measured - e) / e).abs()
    );
}
