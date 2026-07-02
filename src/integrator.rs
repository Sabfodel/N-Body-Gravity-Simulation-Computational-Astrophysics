use crate::body::{Derivative, State};
use crate::forces;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    Rk4,
    Euler,
}

impl std::str::FromStr for Method {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "rk4" => Ok(Method::Rk4),
            "euler" => Ok(Method::Euler),
            other => Err(format!("unknown integrator '{other}' (use rk4 or euler)")),
        }
    }
}

fn derivative(state: &State, masses: &[f64], softening: f64) -> Derivative {
    Derivative {
        dpos: state.velocities.clone(),
        dvel: forces::accelerations(&state.positions, masses, softening),
    }
}

/// One step of the classical 4th-order Runge-Kutta method:
///   y' = f(y),  y(t+dt) = y + dt/6 (k1 + 2k2 + 2k3 + k4)
/// Four force evaluations per step buy O(dt⁴) local accuracy — the energy
/// drift over fixed time shrinks ~16× when dt is halved.
pub fn rk4_step(state: &State, masses: &[f64], dt: f64, softening: f64) -> State {
    let k1 = derivative(state, masses, softening);
    let k2 = derivative(&state.add_scaled(&k1, dt / 2.0), masses, softening);
    let k3 = derivative(&state.add_scaled(&k2, dt / 2.0), masses, softening);
    let k4 = derivative(&state.add_scaled(&k3, dt), masses, softening);

    let n = state.n();
    let mut next = state.clone();
    for i in 0..n {
        next.positions[i] += (k1.dpos[i] + 2.0 * k2.dpos[i] + 2.0 * k3.dpos[i] + k4.dpos[i])
            * (dt / 6.0);
        next.velocities[i] += (k1.dvel[i] + 2.0 * k2.dvel[i] + 2.0 * k3.dvel[i] + k4.dvel[i])
            * (dt / 6.0);
    }
    next
}

/// Explicit (forward) Euler — first order. Included as the baseline the
/// README's energy-drift comparison is measured against.
pub fn euler_step(state: &State, masses: &[f64], dt: f64, softening: f64) -> State {
    let k = derivative(state, masses, softening);
    state.add_scaled(&k, dt)
}

pub fn step(method: Method, state: &State, masses: &[f64], dt: f64, softening: f64) -> State {
    match method {
        Method::Rk4 => rk4_step(state, masses, dt, softening),
        Method::Euler => euler_step(state, masses, dt, softening),
    }
}
