use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// A point mass. JSON fields match JPL Horizons state-vector conventions:
/// heliocentric ecliptic position [AU] and velocity [AU/yr].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Body {
    pub name: String,
    /// Mass in solar masses.
    pub mass: f64,
    /// Position [x, y, z] in AU.
    pub position: [f64; 3],
    /// Velocity [vx, vy, vz] in AU/yr.
    pub velocity: [f64; 3],
}

impl Body {
    pub fn position_vec(&self) -> Vector3<f64> {
        Vector3::from(self.position)
    }

    pub fn velocity_vec(&self) -> Vector3<f64> {
        Vector3::from(self.velocity)
    }
}

/// Dynamical state of the whole system, separated from the (constant) masses
/// so the integrator can treat it as a single vector to combine and scale.
#[derive(Debug, Clone)]
pub struct State {
    pub positions: Vec<Vector3<f64>>,
    pub velocities: Vec<Vector3<f64>>,
}

impl State {
    pub fn from_bodies(bodies: &[Body]) -> Self {
        State {
            positions: bodies.iter().map(Body::position_vec).collect(),
            velocities: bodies.iter().map(Body::velocity_vec).collect(),
        }
    }

    pub fn n(&self) -> usize {
        self.positions.len()
    }

    /// self + k * scale, element-wise. The vector-space operation RK4 is built on.
    pub fn add_scaled(&self, k: &Derivative, scale: f64) -> State {
        State {
            positions: self
                .positions
                .iter()
                .zip(&k.dpos)
                .map(|(p, dp)| p + dp * scale)
                .collect(),
            velocities: self
                .velocities
                .iter()
                .zip(&k.dvel)
                .map(|(v, dv)| v + dv * scale)
                .collect(),
        }
    }

    /// Shift to the centre-of-mass frame: total momentum and the COM itself
    /// become exactly zero, so the system does not drift out of the plot.
    pub fn to_com_frame(&mut self, masses: &[f64]) {
        let m_total: f64 = masses.iter().sum();
        let com: Vector3<f64> = self
            .positions
            .iter()
            .zip(masses)
            .map(|(r, &m)| r * m)
            .sum::<Vector3<f64>>()
            / m_total;
        let v_com: Vector3<f64> = self
            .velocities
            .iter()
            .zip(masses)
            .map(|(v, &m)| v * m)
            .sum::<Vector3<f64>>()
            / m_total;
        for r in &mut self.positions {
            *r -= com;
        }
        for v in &mut self.velocities {
            *v -= v_com;
        }
    }
}

/// Time derivative of a `State`: d(position)/dt and d(velocity)/dt.
#[derive(Debug, Clone)]
pub struct Derivative {
    pub dpos: Vec<Vector3<f64>>,
    pub dvel: Vec<Vector3<f64>>,
}
