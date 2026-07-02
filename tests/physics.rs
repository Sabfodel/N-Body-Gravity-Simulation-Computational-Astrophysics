//! Physics validation: conservation laws and Kepler's laws, not just code paths.

use approx::assert_relative_eq;
use nalgebra::Vector3;
use nbody_sim::analysis;
use nbody_sim::body::State;
use nbody_sim::forces::{self, G};
use nbody_sim::integrator;

const M_EARTH: f64 = 3.0035e-6;

/// Sun–Earth two-body system started at perihelion (a = 1 AU, e = 0.0167),
/// in the centre-of-mass frame.
fn sun_earth() -> (State, Vec<f64>) {
    let a = 1.0;
    let e = 0.0167;
    let rp = a * (1.0 - e);
    let vp = (G * (1.0 + M_EARTH) * (1.0 + e) / (a * (1.0 - e))).sqrt();

    let mut state = State {
        positions: vec![Vector3::zeros(), Vector3::new(rp, 0.0, 0.0)],
        velocities: vec![Vector3::zeros(), Vector3::new(0.0, vp, 0.0)],
    };
    let masses = vec![1.0, M_EARTH];
    state.to_com_frame(&masses);
    (state, masses)
}

#[test]
fn parallel_and_serial_forces_agree() {
    let (state, masses) = sun_earth();
    let par = forces::accelerations(&state.positions, &masses, 1e-4);
    let ser = forces::accelerations_serial(&state.positions, &masses, 1e-4);
    for (a, b) in par.iter().zip(&ser) {
        assert_relative_eq!(a, b, max_relative = 1e-15);
    }
}

#[test]
fn newtons_third_law() {
    // Total force must vanish: sum of m_i * a_i = 0 for any configuration.
    let positions = vec![
        Vector3::new(0.3, -1.2, 0.7),
        Vector3::new(-2.1, 0.4, 0.0),
        Vector3::new(1.0, 1.0, -0.5),
    ];
    let masses = vec![1.0, 0.5, 2.0];
    let acc = forces::accelerations(&positions, &masses, 0.0);
    let net: Vector3<f64> = acc.iter().zip(&masses).map(|(a, &m)| a * m).sum();
    assert!(net.norm() < 1e-12, "net force {:?}", net);
}

#[test]
fn rk4_conserves_energy_and_angular_momentum() {
    let (mut state, masses) = sun_earth();
    let eps = 0.0;
    let dt = 0.001;
    let c0 = analysis::conserved(&state, &masses, eps);

    for _ in 0..10_000 {
        // 10 orbits
        state = integrator::rk4_step(&state, &masses, dt, eps);
    }
    let c = analysis::conserved(&state, &masses, eps);

    assert!(
        analysis::relative_drift(c.total_energy, c0.total_energy) < 1e-10,
        "energy drift {:e}",
        analysis::relative_drift(c.total_energy, c0.total_energy)
    );
    assert_relative_eq!(
        c.angular_momentum.norm(),
        c0.angular_momentum.norm(),
        max_relative = 1e-12
    );
    assert!(c.momentum.norm() < 1e-12, "momentum {:?}", c.momentum);
}

#[test]
fn rk4_beats_euler_by_orders_of_magnitude() {
    let (s0, masses) = sun_earth();
    let dt = 0.001;

    let mut rk4 = s0.clone();
    let mut euler = s0.clone();
    let e0 = analysis::conserved(&s0, &masses, 0.0).total_energy;
    for _ in 0..1000 {
        rk4 = integrator::rk4_step(&rk4, &masses, dt, 0.0);
        euler = integrator::euler_step(&euler, &masses, dt, 0.0);
    }
    let drift_rk4 = analysis::relative_drift(
        analysis::conserved(&rk4, &masses, 0.0).total_energy,
        e0,
    );
    let drift_euler = analysis::relative_drift(
        analysis::conserved(&euler, &masses, 0.0).total_energy,
        e0,
    );
    assert!(
        drift_euler > 1e6 * drift_rk4,
        "euler {drift_euler:e} vs rk4 {drift_rk4:e}"
    );
}

#[test]
fn kepler_period_and_eccentricity() {
    let (mut state, masses) = sun_earth();
    let dt = 1e-4;
    let a: f64 = 1.0;
    let e = 0.0167;
    // Kepler III with the two-body correction: T = 2π sqrt(a³ / G(M+m)).
    let expected_period = (a.powi(3) / (1.0 + M_EARTH)).sqrt();

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
        // Perihelion passage = r stops decreasing. Started at perihelion, so
        // the first minimum after aphelion marks one full period.
        if decreasing && r > prev_r && period.is_none() && step as f64 * dt > 0.5 {
            period = Some(step as f64 * dt);
            break;
        }
        decreasing = r < prev_r;
        prev_r = r;
    }

    let period = period.expect("no perihelion passage found");
    assert_relative_eq!(period, expected_period, max_relative = 2e-4);

    let e_measured = (r_max - r_min) / (r_max + r_min);
    assert_relative_eq!(e_measured, e, max_relative = 1e-2);
}
