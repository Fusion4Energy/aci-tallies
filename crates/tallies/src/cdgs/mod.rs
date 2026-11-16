//! D1SUNED Common Decay Gamma Source (CDGS) files.
//!
//! A CDGS file holds decay gamma sources on meshes: for each mesh, a header
//! (cooling time, source strength, energy bins or lines, the mesh and its
//! orientation), then the voxels with a source, each with its cells' spectra.
//! The layout is D1SUNED 4.1's (`cdgsV2_print`, `cuv_cdgs_print` in its
//! `fmesh_mod.F90`), which F4Enix writes too:
//!
//! ```text
//! num_meshes 1
//! total_source 4.8108E-06
//! mesh_id 1
//! <comment>
//! cooling_time 0.000E+00
//! mesh_source 4.8108E-06
//! energy_type bins                (or lines)
//! energy_boundaries 4
//!  0.000   0.500   1.000   2.000
//! mesh_type rec                   (or cyl)
//! mesh_boundaries 4 3 3           (bounds per axis: x y z; r theta z)
//! <origin, x axis, y axis>        (cyl: origin, AXS, VEC)
//! <bounds of the first axis> ...  (theta in revolutions)
//! source_data
//! <voxel> <source> <volume> <cells>
//! <cell> <volume fraction> <source>   (per cell)
//! <source per energy bin>
//! <relative error per energy bin>
//! end_source_data
//! ```
//!
//! A voxel's source is the sum of its cells' (each cell's is its part of the
//! voxel's, not a density). Voxels are numbered from 1 over the axes in the
//! order printed, the last fastest: x, y, z; for a cylinder r, theta, z,
//! which [`Cdgs::read`] puts in this crate's (r, z, theta) order. Voxels with
//! no source are left out of the file and read as 0.

use std::io;
use std::path::{Path, PathBuf};

use ndarray::{Array3, Array5};

use crate::common::{float, Axis, LabelKind, ReadError, Warning};
use crate::meshtal::{Geometry, OutFormat, Tally, TallyInfo, Transform};

/// Whether `head`, a file's first bytes, is a CDGS file's.
pub fn is_cdgs(head: &[u8]) -> bool {
    String::from_utf8_lossy(head)
        .lines()
        .find(|l| !l.trim().is_empty())
        .is_some_and(|l| l.trim_start().starts_with("num_meshes"))
}

/// What a mesh of a CDGS file is, without its voxels.
#[derive(Debug, Clone)]
pub struct CdgsInfo {
    /// The `mesh_id`.
    pub id: u32,
    /// The line after `mesh_id` (what wrote the file, and when).
    pub comment: String,
    /// Cooling time (s).
    pub cooling_time: f64,
    /// The mesh's source strength (`mesh_source`, photons/s).
    pub strength: f64,
    /// Energy bins (`energy_type bins`, MeV), or lines (labels, MeV), with
    /// a Total when there is more than one.
    pub energy: Axis,
    pub geometry: Geometry,
    pub warnings: Vec<Warning>,
}

struct Mesh {
    info: CdgsInfo,
    /// Lines of the voxel data, `source_data` and `end_source_data` excluded.
    data: std::ops::Range<usize>,
}

/// An indexed CDGS file. It is read whole when opened (CDGS files list only
/// the voxels with a source).
pub struct Cdgs {
    path: Option<PathBuf>,
    lines: Vec<String>,
    total_source: Option<f64>,
    meshes: Vec<Mesh>,
    failed: Vec<(u32, String)>,
}

impl Cdgs {
    pub fn open(path: impl AsRef<Path>) -> io::Result<Cdgs> {
        let bytes = std::fs::read(path.as_ref())?;
        let mut c = Cdgs::from_bytes(&bytes);
        c.path = Some(path.as_ref().to_path_buf());
        Ok(c)
    }

    /// A CDGS file from its contents.
    pub fn from_bytes(bytes: &[u8]) -> Cdgs {
        let lines: Vec<String> = String::from_utf8_lossy(bytes)
            .lines()
            .map(str::to_string)
            .collect();
        let mut c = Cdgs {
            path: None,
            lines,
            total_source: None,
            meshes: Vec::new(),
            failed: Vec::new(),
        };
        c.index();
        c
    }

    pub fn path(&self) -> Option<&Path> {
        self.path.as_deref()
    }

    /// The `total_source` of every mesh (photons/s).
    pub fn total_source(&self) -> Option<f64> {
        self.total_source
    }

    pub fn tallies(&self) -> impl Iterator<Item = &CdgsInfo> {
        self.meshes.iter().map(|m| &m.info)
    }

    pub fn tally(&self, id: u32) -> Result<&CdgsInfo, ReadError> {
        self.mesh(id).map(|m| &m.info)
    }

    /// Meshes that could not be read, and why.
    pub fn failed(&self) -> &[(u32, String)] {
        &self.failed
    }

    fn mesh(&self, id: u32) -> Result<&Mesh, ReadError> {
        if let Some((_, why)) = self.failed.iter().find(|(f, _)| *f == id) {
            return Err(ReadError::Unsupported {
                id,
                reason: why.clone(),
            });
        }
        self.meshes
            .iter()
            .find(|m| m.info.id == id)
            .ok_or_else(|| ReadError::NoTally {
                id: Some(id),
                tallies: self.meshes.iter().map(|m| m.info.id).collect(),
            })
    }

    fn index(&mut self) {
        let mut i = 0;
        while i < self.lines.len() {
            let line = self.lines[i].trim();
            if let Some(v) = keyword(line, "total_source") {
                self.total_source = float::parse(v.as_bytes());
            }
            if let Some(v) = keyword(line, "mesh_id") {
                let id = v.parse().unwrap_or(0);
                let end = (i + 1..self.lines.len())
                    .find(|&n| self.lines[n].trim() == "end_source_data")
                    .unwrap_or(self.lines.len());
                match self.header(i, end, id) {
                    Ok((info, data)) => self.meshes.push(Mesh { info, data }),
                    Err(why) => self.failed.push((id, why)),
                }
                i = end;
            }
            i += 1;
        }
    }

    /// The header of the mesh whose `mesh_id` is line `start`, its block
    /// ending at line `end`; and the lines of its voxels.
    fn header(
        &self,
        start: usize,
        end: usize,
        id: u32,
    ) -> Result<(CdgsInfo, std::ops::Range<usize>), String> {
        let lines = &self.lines[..end];
        let mut warnings = Vec::new();
        let comment = lines.get(start + 1).map_or("", |l| l.trim()).to_string();
        let mut at = start + 2;
        let mut get = |name: &str| -> Result<(usize, String), String> {
            while at < lines.len() {
                let n = at;
                at += 1;
                if let Some(v) = keyword(lines[n].trim(), name) {
                    return Ok((n, v.to_string()));
                }
            }
            Err(format!("no {name} line"))
        };
        let number = |s: &str, name: &str| {
            float::parse(s.trim().as_bytes()).ok_or_else(|| format!("unreadable {name}: {s:?}"))
        };
        let (n, cooling) = get("cooling_time")?;
        let cooling_time = number(&cooling, "cooling_time")?;
        // The mesh's source: `mesh_source` (D1SUNED), `total_source` (F4Enix).
        let source = lines.get(n + 1).map_or("", |l| l.trim());
        let strength = ["mesh_source", "total_source"]
            .iter()
            .find_map(|k| keyword(source, k))
            .ok_or_else(|| format!("no mesh_source line after cooling_time: {source:?}"))
            .and_then(|v| number(v, "mesh_source"))?;
        let (_, energy_type) = get("energy_type")?;
        let (n, count) = get("energy_boundaries")?;
        let n_energy: usize = count
            .trim()
            .parse()
            .map_err(|_| format!("unreadable energy_boundaries {count:?}"))?;
        let (energies, next) = numbers(lines, n + 1, n_energy, true)?;
        if energies.iter().any(|e| e.is_nan()) {
            warnings.push(Warning::at(
                (n + 2) as u64,
                "an energy bound overflowed its field (******): it is NaN",
            ));
        }
        let energy = match energy_type.trim() {
            "bins" => Axis::mcnp(energies),
            "line" | "lines" => Axis::labelled(LabelKind::Line, energies),
            other => return Err(format!("unknown energy_type {other:?}")),
        };
        let mut at2 = next;
        let mut get2 = |name: &str| -> Result<(usize, String), String> {
            while at2 < lines.len() {
                let n = at2;
                at2 += 1;
                if let Some(v) = keyword(lines[n].trim(), name) {
                    return Ok((n, v.to_string()));
                }
            }
            Err(format!("no {name} line"))
        };
        let (_, mesh_type) = get2("mesh_type")?;
        let (n, dims) = get2("mesh_boundaries")?;
        let dims: Vec<usize> = dims
            .split_whitespace()
            .map(|w| {
                w.parse()
                    .map_err(|_| format!("unreadable mesh_boundaries {dims:?}"))
            })
            .collect::<Result<_, _>>()?;
        let [na, nb, nc] = dims[..] else {
            return Err(format!("mesh_boundaries takes 3 numbers, got {dims:?}"));
        };
        let (frame, next) = numbers(lines, n + 1, 9, false)?;
        let (a, next) = numbers(lines, next, na, false)?;
        let (b, next) = numbers(lines, next, nb, false)?;
        let (c, next) = numbers(lines, next, nc, false)?;
        let vector = |k: usize| [frame[3 * k], frame[3 * k + 1], frame[3 * k + 2]];
        let geometry = match mesh_type.trim() {
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
                    z_bounds: c,
                    transform,
                }
            }
            "cyl" => {
                let axis = normalized(vector(1)).ok_or("the cylinder's axis is zero")?;
                let v = vector(2);
                let along = dot(v, axis);
                let vec = normalized([
                    v[0] - along * axis[0],
                    v[1] - along * axis[1],
                    v[2] - along * axis[2],
                ])
                .ok_or("the cylinder's VEC is parallel to its axis")?;
                Geometry::Cylindrical {
                    r_bounds: a,
                    z_bounds: c,
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
        let data_start = (next..end)
            .find(|&n| lines[n].trim() == "source_data")
            .ok_or("no source_data line")?
            + 1;
        Ok((
            CdgsInfo {
                id,
                comment,
                cooling_time,
                strength,
                energy,
                geometry,
                warnings,
            },
            data_start..end,
        ))
    }

    /// Mesh `id` as a meshtal tally: each voxel's source per energy bin (or
    /// line), the sum of its cells', and the Total. Errors are combined as
    /// if the cells were independent; a voxel's Total has the error of the
    /// sum of its bins. `volume` is each voxel's, as the file gives it.
    pub fn read(&self, id: u32) -> Result<Tally, ReadError> {
        let mesh = self.mesh(id)?;
        let info = &mesh.info;
        let unsupported = |reason: String| ReadError::Unsupported { id, reason };
        let geometry = &info.geometry;
        let (ni, nj, nk) = geometry.dims();
        let voxels = ni * nj * nk;
        let n_bins = info.energy.n_bins();
        let ne = info.energy.len();
        // Dimensions in the file's order.
        let file_dims = match geometry {
            Geometry::Rectangular { .. } => [ni, nj, nk],
            Geometry::Cylindrical { .. } => [ni, nk, nj],
        };
        let flat = |n: usize| -> usize {
            let c = n % file_dims[2];
            let b = (n / file_dims[2]) % file_dims[1];
            let a = n / (file_dims[2] * file_dims[1]);
            match geometry {
                Geometry::Rectangular { .. } => (a * nj + b) * nk + c,
                Geometry::Cylindrical { .. } => (a * nj + c) * nk + b,
            }
        };
        let mut value = vec![0.0; ne * voxels];
        let mut var = vec![0.0; ne * voxels];
        let mut volume: Vec<f64> = geometry.voxel_volumes();
        let mut warnings = info.warnings.clone();
        let lines = &self.lines;
        let mut n = mesh.data.start;
        let mut seen = vec![false; voxels];
        while n < mesh.data.end {
            let line = lines[n].trim();
            if line.is_empty() {
                n += 1;
                continue;
            }
            let w: Vec<&str> = line.split_whitespace().collect();
            let bad =
                |what: &str| unsupported(format!("line {}: unreadable {what}: {line:?}", n + 1));
            if w.len() < 4 {
                return Err(bad("voxel line"));
            }
            let index: usize = w[0].parse().map_err(|_| bad("voxel number"))?;
            let vol = float::parse(w[w.len() - 2].as_bytes()).ok_or_else(|| bad("volume"))?;
            let n_cells: usize = w[w.len() - 1].parse().map_err(|_| bad("cell count"))?;
            if index == 0 || index > voxels {
                return Err(bad("voxel number (out of the mesh)"));
            }
            let v = flat(index - 1);
            if seen[v] {
                warnings.push(Warning::at(
                    (n + 1) as u64,
                    format!("voxel {index} listed twice"),
                ));
            }
            seen[v] = true;
            volume[v] = vol;
            n += 1;
            let mut total = 0.0;
            for _ in 0..n_cells {
                let cell = lines.get(n).map_or("", |l| l.trim());
                let cw: Vec<&str> = cell.split_whitespace().collect();
                let cell_total = cw
                    .get(2)
                    .and_then(|t| float::parse(t.as_bytes()))
                    .ok_or_else(|| {
                        unsupported(format!("line {}: unreadable cell line {cell:?}", n + 1))
                    })?;
                total += cell_total;
                let (values, next) = numbers(lines, n + 1, n_bins, false).map_err(unsupported)?;
                let (errors, next) = numbers(lines, next, n_bins, false).map_err(unsupported)?;
                for e in 0..n_bins {
                    value[e * voxels + v] += values[e];
                    var[e * voxels + v] += (values[e] * errors[e]).powi(2);
                }
                n = next;
            }
            if let Some(t) = info.energy.total_index() {
                value[t * voxels + v] = total;
                var[t * voxels + v] = (0..n_bins).map(|e| var[e * voxels + v]).sum();
            }
        }
        let rel: Vec<f64> = value
            .iter()
            .zip(&var)
            .map(|(&x, &s2)| if x != 0.0 { s2.sqrt() / x.abs() } else { 0.0 })
            .collect();
        let shape = (ne, 1, ni, nj, nk);
        let tally_info = TallyInfo {
            id,
            particle: "photon".to_string(),
            comment: Some(info.comment.clone()),
            geometry: geometry.clone(),
            geometry_from_input: false,
            energy_weighted: None,
            energy: info.energy.clone(),
            time: Axis::mcnp(vec![0.0, 1e36]),
            out_format: OutFormat::Cdgs,
            warnings: Vec::new(),
        };
        Ok(Tally::from_info(
            tally_info,
            Array5::from_shape_vec(shape, value).expect("one value per index"),
            Array5::from_shape_vec(shape, rel).expect("one value per index"),
            Some(Array3::from_shape_vec((ni, nj, nk), volume).expect("one volume per voxel")),
            warnings,
        ))
    }
}

/// The value after `name` on a `name value` line.
fn keyword<'a>(line: &'a str, name: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(name)?;
    (rest.is_empty() || rest.starts_with(char::is_whitespace)).then(|| rest.trim())
}

/// `count` numbers from line `start` on, over as many lines as they take,
/// and the line after them. With `overflow`, a field of `*` (a number too
/// large for it) is NaN. Numbers run into each other (`1.0E+00-2.0E+00`)
/// are split at the fixed widths D1SUNED writes (12 and 13 columns).
fn numbers(
    lines: &[String],
    start: usize,
    count: usize,
    overflow: bool,
) -> Result<(Vec<f64>, usize), String> {
    let mut out = Vec::with_capacity(count);
    let mut n = start;
    while out.len() < count {
        let line = lines.get(n).ok_or_else(|| {
            format!(
                "the file ends before {count} numbers from line {}",
                start + 1
            )
        })?;
        let mut row = Vec::new();
        let mut ok = true;
        for tok in line.split_whitespace() {
            match float::parse(tok.as_bytes()) {
                Some(x) => row.push(x),
                None if overflow && tok.bytes().all(|b| b == b'*') => row.push(f64::NAN),
                None => {
                    ok = false;
                    break;
                }
            }
        }
        if !ok {
            row = [12, 13]
                .into_iter()
                .find_map(|w| crate::cuv::fixed_floats(line, w))
                .ok_or_else(|| format!("line {}: unreadable numbers {line:?}", n + 1))?;
        }
        out.extend(row);
        n += 1;
    }
    if out.len() != count {
        return Err(format!(
            "line {}: {} numbers where {count} were due",
            n,
            out.len()
        ));
    }
    Ok((out, n))
}

/// A rectangular mesh's frame, origin and its x and y axes: None for the
/// global one, else the transform (the axes made orthonormal, as they are
/// printed with 6 digits).
fn rectangular_frame(
    origin: [f64; 3],
    x: [f64; 3],
    y: [f64; 3],
) -> Result<Option<Transform>, String> {
    let near = |a: [f64; 3], b: [f64; 3]| a.iter().zip(&b).all(|(p, q)| (p - q).abs() < 1e-9);
    if near(origin, [0.0; 3]) && near(x, [1.0, 0.0, 0.0]) && near(y, [0.0, 1.0, 0.0]) {
        return Ok(None);
    }
    let x = normalized(x).ok_or("the mesh's x axis is zero")?;
    let z = normalized(cross(x, y)).ok_or("the mesh's x and y axes are parallel")?;
    let y = cross(z, x);
    Ok(Some(Transform {
        rotation: [[x[0], y[0], z[0]], [x[1], y[1], z[1]], [x[2], y[2], z[2]]],
        translation: origin,
    }))
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn normalized(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = dot(a, a).sqrt();
    (n > 0.0).then(|| [a[0] / n, a[1] / n, a[2] / n])
}
