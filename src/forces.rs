use nalgebra::Vector3;
use rayon::prelude::*;

/// Gravitational constant in AU³ yr⁻² M☉⁻¹.
pub const G: f64 = 4.0 * std::f64::consts::PI * std::f64::consts::PI;

/// Softened gravitational acceleration on every body, O(N²), one rayon task
/// per body. The softening replaces |r|³ with (|r|² + ε²)^(3/2) so close
/// encounters cannot produce unbounded forces.
pub fn accelerations(
    positions: &[Vector3<f64>],
    masses: &[f64],
    softening: f64,
) -> Vec<Vector3<f64>> {
    let eps2 = softening * softening;
    positions
        .par_iter()
        .enumerate()
        .map(|(i, ri)| acceleration_on(i, ri, positions, masses, eps2))
        .collect()
}

/// Serial reference implementation — identical physics, no rayon. Kept for
/// benchmarking the parallel speedup and for debugging.
pub fn accelerations_serial(
    positions: &[Vector3<f64>],
    masses: &[f64],
    softening: f64,
) -> Vec<Vector3<f64>> {
    let eps2 = softening * softening;
    positions
        .iter()
        .enumerate()
        .map(|(i, ri)| acceleration_on(i, ri, positions, masses, eps2))
        .collect()
}

fn acceleration_on(
    i: usize,
    ri: &Vector3<f64>,
    positions: &[Vector3<f64>],
    masses: &[f64],
    eps2: f64,
) -> Vector3<f64> {
    let mut acc = Vector3::zeros();
    for (j, rj) in positions.iter().enumerate() {
        if j == i {
            continue;
        }
        let d = rj - ri;
        let r2 = d.norm_squared() + eps2;
        acc += d * (G * masses[j] / (r2 * r2.sqrt()));
    }
    acc
}
