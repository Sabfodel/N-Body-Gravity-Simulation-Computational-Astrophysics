use plotters::prelude::*;
use std::path::Path;

/// XY-plane orbit trajectories, one polyline per body, equal axis scaling.
pub fn orbit_plot(
    trajectories: &[(String, Vec<(f64, f64)>)],
    path: &Path,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(dir) = path.parent() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir)?;
        }
    }

    let extent = trajectories
        .iter()
        .flat_map(|(_, pts)| pts.iter())
        .map(|&(x, y)| x.abs().max(y.abs()))
        .fold(1.0_f64, f64::max)
        * 1.1;

    let root = SVGBackend::new(path, (900, 900)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(45)
        .y_label_area_size(60)
        .build_cartesian_2d(-extent..extent, -extent..extent)?;
    chart
        .configure_mesh()
        .x_desc("x [AU]")
        .y_desc("y [AU]")
        .draw()?;

    for (i, (name, pts)) in trajectories.iter().enumerate() {
        let color = Palette99::pick(i).to_rgba();
        chart
            .draw_series(LineSeries::new(pts.iter().copied(), &color))?
            .label(name.clone())
            .legend(move |(x, y)| PathElement::new([(x, y), (x + 18, y)], color));
    }
    chart
        .configure_series_labels()
        .border_style(BLACK)
        .background_style(WHITE.mix(0.85))
        .draw()?;
    root.present()?;
    Ok(())
}

/// log10 of relative energy drift vs time for one or more labelled runs
/// (e.g. RK4 vs Euler on the same system).
pub fn drift_plot(
    series: &[(String, Vec<(f64, f64)>)],
    path: &Path,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(dir) = path.parent() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir)?;
        }
    }

    // Drift is plotted as log10; zero drift (t=0) is clamped to the floor.
    const FLOOR: f64 = -16.0;
    let log_series: Vec<(String, Vec<(f64, f64)>)> = series
        .iter()
        .map(|(label, pts)| {
            let logs = pts
                .iter()
                .map(|&(t, d)| (t, if d > 0.0 { d.log10().max(FLOOR) } else { FLOOR }))
                .collect();
            (label.clone(), logs)
        })
        .collect();

    let t_max = log_series
        .iter()
        .flat_map(|(_, pts)| pts.iter())
        .map(|&(t, _)| t)
        .fold(0.0_f64, f64::max);
    let (y_min, y_max) = log_series
        .iter()
        .flat_map(|(_, pts)| pts.iter())
        .fold((0.0_f64, FLOOR), |(lo, hi), &(_, y)| (lo.min(y), hi.max(y)));

    let root = SVGBackend::new(path, (900, 600)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(title, ("sans-serif", 28))
        .margin(20)
        .x_label_area_size(45)
        .y_label_area_size(60)
        .build_cartesian_2d(0.0..t_max, (y_min - 0.5)..(y_max + 0.5))?;
    chart
        .configure_mesh()
        .x_desc("t [yr]")
        .y_desc("log10 |ΔE / E0|")
        .draw()?;

    for (i, (label, pts)) in log_series.iter().enumerate() {
        let color = Palette99::pick(i).to_rgba();
        chart
            .draw_series(LineSeries::new(pts.iter().copied(), color.stroke_width(2)))?
            .label(label.clone())
            .legend(move |(x, y)| PathElement::new([(x, y), (x + 18, y)], color.stroke_width(2)));
    }
    chart
        .configure_series_labels()
        .position(SeriesLabelPosition::LowerRight)
        .border_style(BLACK)
        .background_style(WHITE.mix(0.85))
        .draw()?;
    root.present()?;
    Ok(())
}
