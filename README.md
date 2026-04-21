# N-Body-Gravity-Simulation-Computational-Astrophysics

Rust License: MIT Build 
A high-performance N-body gravitational simulation written in Rust. Implements 4th-order Runge-Kutta (RK4) numerical integration
with data-parallel force computation via Rayon. Developed as a computational astrophysics research project 

Features
·	RK4 integration — significantly better energy conservation compared to Euler method
·	Parallel force computation — O(N²) loop distributed across all CPU cores via rayon
·	Softening parameter — prevents numerical instability during close encounters: |r|³ → (|r|² + ε²)^(3/2)
·	Energy & momentum tracking — total E = K + U and angular momentum L computed at every step for drift analysis
·	JSON input — initial conditions loaded from NASA JPL Horizons-compatible format
·	CSV output — positions, velocities, and energy written every N steps

Installation
Requirements: Rust 1.75+ (rustup.rs)
git clone https://github.com/USERNAME/nbody-sim
cd nbody-sim
cargo build --release


Usage
cargo run --release -- --input data/solar_system.json --dt 0.001 --years 100

Argument	Description	Default
--input	Initial conditions JSON file	data/solar_system.json
--dt	Time step (years)	0.001
--years	Total simulation duration (years)	100
--softening	Softening length ε (AU)	1e-4
--output	CSV output file	output/sim.csv


Project Structure
nbody-sim/
├── src/
│   ├── main.rs          # Entry point, CLI argument parsing, simulation loop
│   ├── body.rs          # struct Body { mass, position: Vector3, velocity, acceleration }
│   ├── integrator.rs    # RK4 algorithm
│   ├── forces.rs        # Gravitational force computation (Rayon parallel)
│   ├── io.rs            # JSON reader, CSV writer
│   └── analysis.rs      # Energy, momentum, and drift computation
├── data/
│   └── solar_system.json   # NASA JPL Horizons initial conditions
├── assets/
│   ├── orbit_plot.png
│   └── energy_drift.png
├── benches/
│   └── nbody_bench.rs      # Criterion.rs benchmarks
├── Cargo.toml
└── README.md


Results
Kepler Validation
Two-body Sun-Earth system compared against analytical solutions:
Parameter	      Analytical	 Simulation	Error
Orbital period	365.25 days 	 —    	  —
Eccentricity (e)	0.0167	     —	      —

Fill in your own measurements after running the simulation.

Energy Drift Analysis

log|(E(t) - E₀) / E₀| — RK4 vs Euler comparison across time:
Orbital Plots
Benchmark — Rust vs Python
500 bodies, 10-year evolution, Δt = 0.001 yr:
Method	           Time (s)	Speedup
Python scipy.odeint   	—	1×
Rust serial	            —	—×
Rust + Rayon           	—	—×

Fill in your own benchmark results.

Dependencies
Crate	Version	Purpose
nalgebra	0.32	3D vector and matrix operations
rayon	1.8	Data parallelism
serde + serde_json	1.0	JSON deserialization
plotters	0.3	Plot generation


Theoretical Background
  The equation of motion for body i in an N-body system:
\ddot{\mathbf{r}}_i = -G \sum_{j \neq i} \frac{m_j (\mathbf{r}_i - \mathbf{r}_j)}{|\mathbf{r}_i - \mathbf{r}_j|^3}
Softening correction applied to avoid singularities at small separations:
|\mathbf{r}|^3 \rightarrow \left(|\mathbf{r}|^2 + \varepsilon^2\right)^{3/2}
The simulation integrates this system using the 4th-order Runge-Kutta method with a fixed timestep \Delta t. Total energy E = K + U and angular momentum \mathbf{L} are tracked throughout to validate numerical accuracy.
For full theoretical background and methodology, see the project report.

References
·	Carroll, B. W. & Ostlie, D. A. (2017). An Introduction to Modern Astrophysics. Cambridge University Press.
·	Hairer, E., Nørsett, S. P. & Wanner, G. (1993). Solving Ordinary Differential Equations I. Springer.
·	NASA JPL Horizons System: https://ssd.jpl.nasa.gov/horizons/
·	Rust Programming Language — Official Documentation: https://doc.rust-lang.org

License
MIT © 
