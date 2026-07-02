//! N-body gravitational simulation.
//!
//! Units throughout: length in AU, time in years, mass in solar masses.
//! In these units G = 4π² AU³ yr⁻² M☉⁻¹, so a body on a circular orbit of
//! semi-major axis 1 AU around 1 M☉ has period 1 yr (Kepler's third law).

pub mod analysis;
pub mod body;
pub mod forces;
pub mod integrator;
pub mod io;
pub mod plot;
