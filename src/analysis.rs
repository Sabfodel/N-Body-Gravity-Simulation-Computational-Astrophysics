use crate::body::State;
use crate::forces::G;
use nalgebra::Vector3;

/// Conserved quantities of the system at one instant.
#[derive(Debug, Clone, Copy)]
pub struct Conserved {
    pub kinetic: f64,
    pub potential: f64,
    pub total_energy: f64,
    pub momentum: Vector3<f64>,
    pub angular_momentum: Vector3<f64>,
}

/// Energy, momentum and angular momentum. The potential uses the *softened*
/// form -G m_i m_j / sqrt(r² + ε²), whose gradient is exactly the force in
/// `forces.rs` — so E = K + U is conserved by the softened dynamics too.
pub fn conserved(state: &State, masses: &[f64], softening: f64) -> Conserved {
    let eps2 = softening * softening;

    let kinetic: f64 = state
        .velocities
        .iter()
        .zip(masses)
        .map(|(v, &m)| 0.5 * m * v.norm_squared())
        .sum();

    let mut potential = 0.0;
    for i in 0..state.n() {
        for j in (i + 1)..state.n() {
            let r2 = (state.positions[i] - state.positions[j]).norm_squared() + eps2;
            potential -= G * masses[i] * masses[j] / r2.sqrt();
        }
    }

    let momentum: Vector3<f64> = state
        .velocities
        .iter()
        .zip(masses)
        .map(|(v, &m)| v * m)
        .sum();

    let angular_momentum: Vector3<f64> = state
        .positions
        .iter()
        .zip(&state.velocities)
        .zip(masses)
        .map(|((r, v), &m)| m * r.cross(v))
        .sum();

    Conserved {
        kinetic,
        potential,
        total_energy: kinetic + potential,
        momentum,
        angular_momentum,
    }
}

/// Relative energy drift |(E - E0) / E0| — the standard integrator quality metric.
pub fn relative_drift(e: f64, e0: f64) -> f64 {
    ((e - e0) / e0).abs()
}
