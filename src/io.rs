use crate::analysis::Conserved;
use crate::body::{Body, State};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub fn load_bodies(path: &Path) -> Result<Vec<Body>, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let bodies: Vec<Body> = serde_json::from_str(&text)
        .map_err(|e| format!("cannot parse {}: {e}", path.display()))?;
    if bodies.len() < 2 {
        return Err("need at least two bodies".into());
    }
    Ok(bodies)
}

/// Writes two CSV files: `<output>` with one row per body per snapshot
/// (t, name, x, y, z, vx, vy, vz) and `<stem>_energy.csv` with one row per
/// snapshot (t, kinetic, potential, total, |L|, relative drift).
pub struct CsvWriter {
    trajectories: BufWriter<File>,
    energy: BufWriter<File>,
}

impl CsvWriter {
    pub fn create(output: &Path) -> Result<Self, Box<dyn std::error::Error>> {
        if let Some(dir) = output.parent() {
            if !dir.as_os_str().is_empty() {
                std::fs::create_dir_all(dir)?;
            }
        }
        let energy_path = output.with_file_name(format!(
            "{}_energy.csv",
            output.file_stem().unwrap_or_default().to_string_lossy()
        ));

        let mut trajectories = BufWriter::new(File::create(output)?);
        writeln!(trajectories, "t,name,x,y,z,vx,vy,vz")?;
        let mut energy = BufWriter::new(File::create(&energy_path)?);
        writeln!(energy, "t,kinetic,potential,total,angular_momentum,drift")?;

        Ok(CsvWriter {
            trajectories,
            energy,
        })
    }

    pub fn snapshot(
        &mut self,
        t: f64,
        bodies: &[Body],
        state: &State,
        c: &Conserved,
        drift: f64,
    ) -> std::io::Result<()> {
        for (i, body) in bodies.iter().enumerate() {
            let r = state.positions[i];
            let v = state.velocities[i];
            writeln!(
                self.trajectories,
                "{t:.6},{},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9}",
                body.name, r.x, r.y, r.z, v.x, v.y, v.z
            )?;
        }
        writeln!(
            self.energy,
            "{t:.6},{:.12e},{:.12e},{:.12e},{:.12e},{:.6e}",
            c.kinetic,
            c.potential,
            c.total_energy,
            c.angular_momentum.norm(),
            drift
        )
    }

    pub fn finish(mut self) -> std::io::Result<()> {
        self.trajectories.flush()?;
        self.energy.flush()
    }
}
