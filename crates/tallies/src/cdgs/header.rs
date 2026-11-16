//! A CDGS mesh's header: the lines from `mesh_id` to `source_data`.

use aci_tallies_geom::vec3::{cross, normalize, perpendicular_unit, Vec3};

use super::CdgsInfo;
use crate::common::float::{self, fixed_floats};
use crate::common::{scan, Axis, LabelKind, Warning};
use crate::meshtal::{Geometry, Transform};

/// The value after `name` on a `name value` line (`line` trimmed).
pub(super) fn keyword<'a>(line: &'a str, name: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(name)?;
    (rest.is_empty() || rest.starts_with(char::is_whitespace)).then(|| rest.trim())
}

/// The mesh number on a `mesh_id` line.
pub(super) fn mesh_id(line: &str) -> Result<u32, String> {
    let v = keyword(line.trim(), "mesh_id").unwrap_or("");
    v.parse()
        .map_err(|_| format!("unreadable mesh_id {:?}", line.trim()))
}

/// Lines of a header, read forwards.
struct Cursor<'a> {
    lines: &'a [&'a str],
    /// The next line to read.
    at: usize,
    /// Lines of the file before the header.
    line0: u64,
}

impl<'a> Cursor<'a> {
    /// The value of the next line that starts with one of `keys`.
    fn value(&mut self, keys: &[&str]) -> Result<&'a str, String> {
        while self.at < self.lines.len() {
            let line = self.lines[self.at].trim();
            self.at += 1;
            if let Some(v) = keys.iter().find_map(|k| keyword(line, k)) {
                return Ok(v);
            }
        }
        Err(format!("no {} line", keys[0]))
    }

    /// File line number (from 1) of the line last read.
    fn line_no(&self) -> u64 {
        self.line0 + self.at as u64
    }

    /// `count` numbers from the next line on, over as many lines as they
    /// take. With `overflow`, a field of `*` (a number too large for it) is
    /// NaN. Numbers run into each other (`1.0E+00-2.0E+00`) are split at
    /// the fixed widths D1SUNED writes (12 and 13 columns).
    fn numbers(&mut self, count: usize, overflow: bool) -> Result<Vec<f64>, String> {
        let mut out = Vec::with_capacity(count);
        while out.len() < count {
            let Some(line) = self.lines.get(self.at) else {
                return Err(format!(
                    "the header ends {} numbers short (line {})",
                    count - out.len(),
                    self.line_no()
                ));
            };
            self.at += 1;
            numbers_of(line, overflow, &mut out)
                .ok_or_else(|| format!("line {}: unreadable numbers {line:?}", self.line_no()))?;
        }
        if out.len() != count {
            return Err(format!(
                "line {}: {} numbers where {count} were due",
                self.line_no(),
                out.len()
            ));
        }
        Ok(out)
    }
}

/// The numbers of `line`, appended to `out`; None when it has something
/// else. See [`Cursor::numbers`].
pub(super) fn numbers_of(line: &str, overflow: bool, out: &mut Vec<f64>) -> Option<()> {
    let start = out.len();
    for tok in scan::tokens(line.as_bytes()) {
        match float::parse(tok) {
            Some(x) => out.push(x),
            None if overflow && tok.iter().all(|&b| b == b'*') => out.push(f64::NAN),
            None => {
                out.truncate(start);
                let row = [12, 13].into_iter().find_map(|w| fixed_floats(line, w))?;
                out.extend(row);
                return Some(());
            }
        }
    }
    Some(())
}

/// The header of the mesh whose `mesh_id` line is the first of `lines`,
/// line `line0 + 1` of the file.
pub(super) fn parse(lines: &[&str], line0: u64) -> Result<CdgsInfo, String> {
    let id = mesh_id(lines.first().copied().unwrap_or(""))?;
    let comment = lines.get(1).map_or("", |l| l.trim()).to_string();
    let mut c = Cursor {
        lines,
        at: 2,
        line0,
    };
    let number = |s: &str, name: &str| {
        float::parse(s.as_bytes()).ok_or_else(|| format!("unreadable {name} {s:?}"))
    };
    let mut warnings = Vec::new();
    let cooling_time = number(c.value(&["cooling_time"])?, "cooling_time")?;
    // The mesh's source: `mesh_source` (D1SUNED), `total_source` (F4Enix).
    let strength = number(c.value(&["mesh_source", "total_source"])?, "mesh_source")?;
    let energy_type = c.value(&["energy_type"])?;
    let count = c.value(&["energy_boundaries"])?;
    let n_energy: usize = count
        .parse()
        .map_err(|_| format!("unreadable energy_boundaries {count:?}"))?;
    let energies = c.numbers(n_energy, true)?;
    if energies.iter().any(|e| e.is_nan()) {
        warnings.push(Warning::at(
            c.line_no(),
            "an energy bound overflowed its field (******): it is NaN",
        ));
    }
    let energy = match energy_type {
        "bins" => Axis::mcnp(energies),
        "line" | "lines" => Axis::labelled(LabelKind::Line, energies),
        other => return Err(format!("unknown energy_type {other:?}")),
    };
    let mesh_type = c.value(&["mesh_type"])?;
    let dims = c.value(&["mesh_boundaries"])?;
    let n: Vec<usize> = dims
        .split_whitespace()
        .map(|w| w.parse())
        .collect::<Result<_, _>>()
        .map_err(|_| format!("unreadable mesh_boundaries {dims:?}"))?;
    let [na, nb, nc] = n[..] else {
        return Err(format!("mesh_boundaries takes 3 numbers, got {dims:?}"));
    };
    let frame = c.numbers(9, false)?;
    let a = c.numbers(na, false)?;
    let b = c.numbers(nb, false)?;
    let cc = c.numbers(nc, false)?;
    let vector = |k: usize| [frame[3 * k], frame[3 * k + 1], frame[3 * k + 2]];
    let geometry = match mesh_type {
        "rec" => {
            let transform = rectangular_frame(vector(0), vector(1), vector(2))?;
            if transform.is_some() {
                warnings.push(Warning::new(
                    "the mesh is placed by a TR card: its frame is read as origin, \
                     x axis and y axis, which no file has confirmed yet",
                ));
            }
            Geometry::Rectangular {
                x_bounds: a,
                y_bounds: b,
                z_bounds: cc,
                transform,
            }
        }
        // Printed r, theta, z: this crate's order is r, z, theta.
        "cyl" => {
            let axis = normalize(vector(1)).ok_or("the cylinder's axis is zero")?;
            let vec = perpendicular_unit(vector(2), axis)
                .ok_or("the cylinder's VEC is parallel to its axis")?;
            Geometry::Cylindrical {
                r_bounds: a,
                z_bounds: cc,
                theta_bounds: b,
                origin: vector(0),
                axis,
                vec,
            }
        }
        other => return Err(format!("unknown mesh_type {other:?}")),
    };
    geometry
        .validate()
        .map_err(|e| format!("invalid mesh: {e}"))?;
    Ok(CdgsInfo {
        id,
        comment,
        cooling_time,
        strength,
        energy,
        geometry,
        warnings,
    })
}

/// A rectangular mesh's frame, origin and its x and y axes: None for the
/// global one, else the transform (the axes made orthonormal, as they are
/// printed with 6 digits).
fn rectangular_frame(origin: Vec3, x: Vec3, y: Vec3) -> Result<Option<Transform>, String> {
    let near = |a: Vec3, b: Vec3| a.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-9);
    if near(origin, [0.0; 3]) && near(x, [1.0, 0.0, 0.0]) && near(y, [0.0, 1.0, 0.0]) {
        return Ok(None);
    }
    let x = normalize(x).ok_or("the mesh's x axis is zero")?;
    let z = normalize(cross(x, y)).ok_or("the mesh's x and y axes are parallel")?;
    let y = cross(z, x);
    Ok(Some(Transform {
        rotation: [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]],
        translation: origin,
    }))
}
