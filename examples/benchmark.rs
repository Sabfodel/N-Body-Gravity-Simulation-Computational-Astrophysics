//! Serial vs Rayon-parallel force computation benchmark.
//!
//! Run: cargo run --release --example benchmark

use nalgebra::Vector3;
use nbody_sim::forces;
use std::time::Instant;

/// Deterministic pseudo-random cluster of N bodies (xorshift, no rand dependency).
fn make_cluster(n: usize) -> (Vec<Vector3<f64>>, Vec<f64>) {
    let mut seed = 0x2545F4914F6CDD1Du64;
    let mut next = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    let positions = (0..n)
        .map(|_| Vector3::new(next() * 20.0, next() * 20.0, next() * 20.0))
        .collect();
    let masses = (0..n).map(|_| 1e-6 + next().abs() * 1e-4).collect();
    (positions, masses)
}

fn time<F: FnMut()>(mut f: F, reps: usize) -> f64 {
    // Warm-up (thread pool spin-up, caches), then best-of-reps wall time.
    f();
    (0..reps)
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::INFINITY, f64::min)
}

fn main() {
    println!(
        "{} logical cores | best of 20 runs per case\n",
        std::thread::available_parallelism().map_or(0, |n| n.get())
    );
    println!("{:>6} {:>14} {:>14} {:>9}", "N", "serial [ms]", "rayon [ms]", "speedup");
    for n in [100, 500, 1000, 2000] {
        let (positions, masses) = make_cluster(n);
        let serial = time(
            || {
                std::hint::black_box(forces::accelerations_serial(&positions, &masses, 1e-4));
            },
            20,
        );
        let parallel = time(
            || {
                std::hint::black_box(forces::accelerations(&positions, &masses, 1e-4));
            },
            20,
        );
        println!(
            "{n:>6} {:>14.3} {:>14.3} {:>8.1}x",
            serial * 1e3,
            parallel * 1e3,
            serial / parallel
        );
    }
}
