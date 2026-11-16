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
//! <cell> <volume fraction> <source> [<n>]   (per cell)
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
//!
//! Like the other readers, opening indexes the file in one parallel pass
//! ([`crate::common::scan`]: meshes start at their `mesh_id` line) and
//! parses the headers; a read cuts the mesh's voxels into segments of whole
//! voxels (found on its first read) and parses them in parallel. Only the
//! headers stay in memory.

mod header;

use std::io;
use std::ops::Range;
use std::path::Path;
use std::sync::{Arc, OnceLock};

use ndarray::{Array3, Array5};
use rayon::prelude::*;

use crate::common::scan::{self, Chunk, Segment};
use crate::common::source::{Opened, Source};
use crate::common::{float, Axis, ReadError, Warning};
use crate::meshtal::{Geometry, OutFormat, Tally, TallyInfo};

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
    /// The voxel lines (bytes), `source_data` and `end_source_data` excluded.
    data: Range<u64>,
    /// Lines of the file before `data`.
    data_line: u64,
    segments: OnceLock<Arc<Vec<Segment>>>,
}

/// An indexed CDGS file.
pub struct Cdgs {
    source: Source,
    chunks: Vec<Chunk>,
    total_source: Option<f64>,
    meshes: Vec<Mesh>,
    failed: Vec<(u32, String)>,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    MeshId,
    SourceData,
    EndSourceData,
}

impl Cdgs {
    /// Index the CDGS file at `path`.
    pub fn open(path: impl AsRef<Path>) -> io::Result<Cdgs> {
        Cdgs::build(Source::file(path.as_ref())?, scan::CHUNK_BYTES)
    }

    /// A CDGS file held in memory.
    pub fn from_bytes(bytes: impl Into<Arc<[u8]>>) -> Cdgs {
        Cdgs::build(Source::bytes(bytes), scan::CHUNK_BYTES).expect("reading memory doesn't fail")
    }

    pub(crate) fn build(source: Source, chunk: usize) -> io::Result<Cdgs> {
        let opened = source.open()?;
        let index = scan::index(&opened, chunk, |line| {
            let t = std::str::from_utf8(line).ok()?.trim();
            if header::keyword(t, "mesh_id").is_some() {
                Some(Kind::MeshId)
            } else if t == "source_data" {
                Some(Kind::SourceData)
            } else if t == "end_source_data" {
                Some(Kind::EndSourceData)
            } else {
                None
            }
        })?;
        let marks = &index.marks;
        let starts: Vec<usize> = (0..marks.len())
            .filter(|&n| marks[n].kind == Kind::MeshId)
            .collect();
        let preamble_end = starts.first().map_or(opened.len(), |&n| marks[n].start);
        let preamble = opened.read_upto(0..preamble_end.min(64 << 10))?;
        let total_source = String::from_utf8_lossy(&preamble)
            .lines()
            .find_map(|l| header::keyword(l.trim(), "total_source"))
            .and_then(|v| float::parse(v.as_bytes()));
        let parsed = (0..starts.len())
            .into_par_iter()
            .map(|s| {
                let first = starts[s];
                let last = starts.get(s + 1).map_or(marks.len(), |&n| n);
                let end = marks.get(last).map_or(opened.len(), |m| m.start);
                index_mesh(&opened, &marks[first..last], end)
            })
            .collect::<io::Result<Vec<_>>>()?;
        drop(opened);
        let mut meshes = Vec::new();
        let mut failed = Vec::new();
        for p in parsed {
            match p {
                Ok(m) => meshes.push(m),
                Err(f) => failed.push(f),
            }
        }
        Ok(Cdgs {
            source,
            chunks: index.chunks,
            total_source,
            meshes,
            failed,
        })
    }

    pub fn path(&self) -> Option<&Path> {
        self.source.path()
    }

    /// Whether the file is unchanged since it was opened.
    pub fn is_current(&self) -> bool {
        self.source.is_current()
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

    /// Mesh `id` as a meshtal tally: each voxel's source per energy bin (or
    /// line), the sum of its cells', and the Total, the sum of its bins.
    /// Errors are combined as if the cells were independent. `volume` is
    /// each voxel's, as the file gives it.
    pub fn read(&self, id: u32) -> Result<Tally, ReadError> {
        let mesh = self.mesh(id)?;
        let info = &mesh.info;
        let opened = self.source.open()?;
        let segments = self.segments(mesh, &opened)?;
        let geometry = &info.geometry;
        let (ni, nj, nk) = geometry.dims();
        let layout = Layout {
            n_bins: info.energy.n_bins(),
            voxels: ni * nj * nk,
            // Dimensions in the file's order.
            file_dims: match geometry {
                Geometry::Rectangular { .. } => [ni, nj, nk],
                Geometry::Cylindrical { .. } => [ni, nk, nj],
            },
            cylinder: matches!(geometry, Geometry::Cylindrical { .. }),
        };
        let Layout { n_bins, voxels, .. } = layout;
        let ne = info.energy.len();
        let mut value = vec![0.0; ne * voxels];
        let mut var = vec![0.0; ne * voxels];
        let mut volume: Vec<f64> = geometry.voxel_volumes();
        let mut seen = vec![false; voxels];
        let mut warnings = info.warnings.clone();
        let mut twice = Tally1::default();
        let mut mismatch = Tally1::default();
        // Segments are parsed in parallel a batch at a time, each batch put
        // in place before the next is parsed: memory beyond the result is
        // a batch's, not the mesh's.
        let batch = 4 * rayon::current_num_threads();
        for group in segments.chunks(batch) {
            let parts = group
                .par_iter()
                .map(|s| {
                    let text = opened.read(s.bytes.clone())?;
                    parse_segment(&text, s.line, &layout)
                        .map_err(|reason| ReadError::Unsupported { id, reason })
                })
                .collect::<Result<Vec<Voxels>, ReadError>>()?;
            for p in &parts {
                for (n, &v) in p.voxel.iter().enumerate() {
                    if std::mem::replace(&mut seen[v], true) {
                        twice.add(p.line[n]);
                    }
                    volume[v] = p.volume[n];
                    let bins = n * n_bins..(n + 1) * n_bins;
                    let mut sum = (0.0, 0.0);
                    for (e, (&x, &s2)) in p.value[bins.clone()].iter().zip(&p.var[bins]).enumerate()
                    {
                        value[e * voxels + v] = x;
                        var[e * voxels + v] = s2;
                        sum = (sum.0 + x, sum.1 + s2);
                    }
                    if let Some(t) = info.energy.total_index() {
                        value[t * voxels + v] = sum.0;
                        var[t * voxels + v] = sum.1;
                    }
                    // What the file says the voxel's source is, against the sum
                    // of its bins (both printed with 6 digits).
                    let off = |s: f64| (s - sum.0).abs() > 1e-4 * s.abs().max(sum.0.abs());
                    if off(p.source[n]) || off(p.cells_source[n]) {
                        mismatch.add(p.line[n]);
                    }
                }
            }
        }
        drop(opened);
        if let Some(line) = twice.first {
            warnings.push(Warning::at(
                line,
                format!(
                    "{} voxels listed more than once (the last listing is kept)",
                    twice.count
                ),
            ));
        }
        if let Some(line) = mismatch.first {
            warnings.push(Warning::at(
                line,
                format!(
                    "{} voxels' printed source (or their cells' sum) differs from the sum of \
                     their energy bins by more than 1e-4",
                    mismatch.count
                ),
            ));
        }
        // The variances become relative errors in place.
        let mut rel = var;
        for (s2, &x) in rel.iter_mut().zip(&value) {
            *s2 = if x != 0.0 { s2.sqrt() / x.abs() } else { 0.0 };
        }
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

    /// The mesh's segments, found on its first read.
    fn segments(&self, mesh: &Mesh, opened: &Opened<'_>) -> Result<Arc<Vec<Segment>>, ReadError> {
        if let Some(s) = mesh.segments.get() {
            return Ok(s.clone());
        }
        let segments = Arc::new(scan::segments(
            opened,
            &self.chunks,
            mesh.data.clone(),
            mesh.data_line,
            is_voxel_line,
        )?);
        let _ = mesh.segments.set(segments.clone());
        Ok(segments)
    }
}

/// The mesh whose marks (`mesh_id` first) are `marks`, its block ending at
/// byte `end`; or its id and why it can't be read.
fn index_mesh(
    opened: &Opened<'_>,
    marks: &[scan::Mark<Kind>],
    end: u64,
) -> io::Result<Result<Mesh, (u32, String)>> {
    let start = &marks[0];
    let find = |k: Kind| marks.iter().find(|m| m.kind == k);
    let Some(data) = find(Kind::SourceData) else {
        let line = opened.read_upto(start.start..start.end)?;
        let id = header::mesh_id(&String::from_utf8_lossy(&line)).unwrap_or(0);
        return Ok(Err((id, "no source_data line".into())));
    };
    let text = opened.read_upto(start.start..data.start)?;
    let text = String::from_utf8_lossy(&text);
    let lines: Vec<&str> = text.lines().collect();
    let info = match header::parse(&lines, start.line) {
        Ok(info) => info,
        Err(why) => {
            let id = header::mesh_id(lines.first().copied().unwrap_or("")).unwrap_or(0);
            return Ok(Err((id, why)));
        }
    };
    let mut info = info;
    let data_end = match find(Kind::EndSourceData) {
        Some(m) => m.start,
        None => {
            info.warnings.push(Warning::at(
                data.line + 1,
                "no end_source_data line: the file may be cut short",
            ));
            end
        }
    };
    Ok(Ok(Mesh {
        info,
        data: data.end..data_end,
        data_line: data.line + 1,
        segments: OnceLock::new(),
    }))
}

/// Where voxels go.
struct Layout {
    n_bins: usize,
    voxels: usize,
    file_dims: [usize; 3],
    cylinder: bool,
}

impl Layout {
    /// This crate's flat index of the voxel numbered `n` (from 0) in the file.
    fn flat(&self, n: usize) -> usize {
        let [_, d1, d2] = self.file_dims;
        let (a, b, c) = (n / (d1 * d2), (n / d2) % d1, n % d2);
        if self.cylinder {
            // File (r, theta, z); here (r, z, theta).
            (a * d2 + c) * d1 + b
        } else {
            (a * d1 + b) * d2 + c
        }
    }
}

/// A count of something, and the line it was first seen on.
#[derive(Default)]
struct Tally1 {
    count: usize,
    first: Option<u64>,
}

impl Tally1 {
    fn add(&mut self, line: u64) {
        self.count += 1;
        self.first.get_or_insert(line);
    }
}

/// The voxels of a segment: one entry each, `n_bins` values each in
/// `value` and `var` (absolute variance).
#[derive(Default)]
struct Voxels {
    voxel: Vec<usize>,
    line: Vec<u64>,
    volume: Vec<f64>,
    /// The voxel line's source.
    source: Vec<f64>,
    /// The sum of its cell lines' sources.
    cells_source: Vec<f64>,
    value: Vec<f64>,
    var: Vec<f64>,
}

/// `<voxel> <source> <volume> <cells>`: two whole numbers around two
/// numbers with exponents. A cell line's second number is a fraction with
/// none (`0  1.00000  6.08415E-04 1`), and value and error lines have no
/// whole numbers.
fn is_voxel_line(line: &[u8]) -> bool {
    let mut t = scan::tokens(line);
    let (Some(a), Some(b), Some(c), Some(d), None) =
        (t.next(), t.next(), t.next(), t.next(), t.next())
    else {
        return false;
    };
    let whole = |t: &[u8]| t.iter().all(u8::is_ascii_digit);
    let exp = |t: &[u8]| t.iter().any(|&x| x == b'E' || x == b'e');
    whole(a) && whole(d) && exp(b) && exp(c)
}

/// What the parser expects next in a voxel.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Expect {
    Cell,
    Values,
    Errors,
    Done,
}

/// The voxels of a segment, `text`, which starts after line `line0`.
fn parse_segment(text: &[u8], line0: u64, layout: &Layout) -> Result<Voxels, String> {
    let n_bins = layout.n_bins;
    let mut out = Voxels::default();
    let mut expect = Expect::Done;
    let mut cells_left = 0usize;
    let mut buf: Vec<f64> = Vec::with_capacity(n_bins);
    let mut values: Vec<f64> = Vec::with_capacity(n_bins);
    let mut voxel_line = line0 + 1;
    for (n, line) in scan::lines(text).enumerate() {
        let line_no = line0 + n as u64 + 1;
        if scan::is_blank(line) {
            continue;
        }
        let at = |m: String| format!("line {line_no}: {m}");
        let text = || String::from_utf8_lossy(line).trim().to_string();
        if is_voxel_line(line) {
            if expect != Expect::Done {
                return Err(format!("line {voxel_line}: the voxel's cells end early"));
            }
            let tok: Vec<&[u8]> = scan::tokens(line).collect();
            let int = |t: &[u8]| std::str::from_utf8(t).ok()?.parse::<usize>().ok();
            let index = int(tok[0]).ok_or_else(|| at(format!("unreadable voxel {:?}", text())))?;
            if index == 0 || index > layout.voxels {
                return Err(at(format!(
                    "voxel {index} is not in the mesh of {} voxels",
                    layout.voxels
                )));
            }
            let num =
                |t: &[u8]| float::parse(t).ok_or_else(|| at(format!("unreadable {:?}", text())));
            out.voxel.push(layout.flat(index - 1));
            out.line.push(line_no);
            out.source.push(num(tok[1])?);
            out.volume.push(num(tok[2])?);
            out.cells_source.push(0.0);
            out.value.extend(std::iter::repeat_n(0.0, n_bins));
            out.var.extend(std::iter::repeat_n(0.0, n_bins));
            cells_left =
                int(tok[3]).ok_or_else(|| at(format!("unreadable cell count {:?}", text())))?;
            expect = if cells_left == 0 {
                Expect::Done
            } else {
                Expect::Cell
            };
            voxel_line = line_no;
            continue;
        }
        match expect {
            Expect::Done if out.voxel.is_empty() => {
                return Err(at(format!("a line before the first voxel: {:?}", text())))
            }
            Expect::Done => {
                return Err(at(format!(
                    "a line after the voxel's cells (a voxel number too large for its \
                     field?): {:?}",
                    text()
                )))
            }
            Expect::Cell => {
                // `<cell> <fraction> <source> [<n>]`.
                let source = scan::tokens(line)
                    .nth(2)
                    .and_then(float::parse)
                    .ok_or_else(|| at(format!("unreadable cell line {:?}", text())))?;
                *out.cells_source.last_mut().expect("in a voxel") += source;
                buf.clear();
                expect = Expect::Values;
            }
            Expect::Values | Expect::Errors => {
                let s = std::str::from_utf8(line).map_err(|_| at("not text".into()))?;
                header::numbers_of(s, false, &mut buf)
                    .ok_or_else(|| at(format!("unreadable numbers {:?}", text())))?;
                if buf.len() > n_bins {
                    return Err(at(format!("{} numbers where {n_bins} were due", buf.len())));
                }
                if buf.len() < n_bins {
                    continue;
                }
                if expect == Expect::Values {
                    std::mem::swap(&mut values, &mut buf);
                    buf.clear();
                    expect = Expect::Errors;
                    continue;
                }
                let base = out.value.len() - n_bins;
                for e in 0..n_bins {
                    out.value[base + e] += values[e];
                    out.var[base + e] += (values[e] * buf[e]).powi(2);
                }
                cells_left -= 1;
                expect = if cells_left == 0 {
                    Expect::Done
                } else {
                    Expect::Cell
                };
            }
        }
    }
    if expect != Expect::Done {
        return Err(format!("line {voxel_line}: the voxel's cells end early"));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> Vec<u8> {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../tests/fixtures/cdgs/cdgs"
        );
        std::fs::read(path).unwrap()
    }

    #[test]
    fn voxel_lines_by_shape() {
        assert!(is_voxel_line(b"      1   6.08415E-04   4.00000E+03    1"));
        // D1SUNED 4.1's cell line, and the 3-number one F4Enix writes.
        assert!(!is_voxel_line(b" 0       1.00000   6.08415E-04 1"));
        assert!(!is_voxel_line(b"      0  1.00000   3.00000E+00"));
        assert!(!is_voxel_line(
            b" 4.51891E-04 1.55796E-04 7.28825E-07 1.0E+00"
        ));
        assert!(!is_voxel_line(b"0.008012 0.014014 0.035603 0.1"));
    }

    /// Every chunk size gives the same values: segments split anywhere.
    #[test]
    fn same_for_any_chunk() {
        let bytes = fixture();
        let whole = Cdgs::build(Source::bytes(bytes.clone()), 1 << 20).unwrap();
        for chunk in [1, 7, 64, 300] {
            let f = Cdgs::build(Source::bytes(bytes.clone()), chunk).unwrap();
            assert_eq!(f.failed(), whole.failed());
            for info in whole.tallies() {
                let (a, b) = (whole.read(info.id).unwrap(), f.read(info.id).unwrap());
                assert_eq!(a.result, b.result, "mesh {} chunk {chunk}", info.id);
                assert_eq!(a.rel_error, b.rel_error);
                assert_eq!(a.volume, b.volume);
            }
        }
    }

    fn one_mesh(data: &str) -> Cdgs {
        let text = format!(
            "num_meshes 1\ntotal_source 1.0\nmesh_id 7\nx\ncooling_time 0.0\nmesh_source 3.0\n\
             energy_type bins\nenergy_boundaries 3\n 0.0 1.0 2.0\nmesh_type rec\n\
             mesh_boundaries 3 2 2\n 0 0 0 1 0 0 0 1 0\n 0 1 2\n 0 1\n 0 1\nsource_data\n\
             {data}end_source_data\n"
        );
        Cdgs::from_bytes(text.into_bytes())
    }

    #[test]
    fn checks_what_the_file_says() {
        let voxel =
            "  1  3.00000E+00  1.00000E+00  1\n  5  1.0  3.00000E+00\n 1.0E+00 2.0E+00\n0.1 0.1\n";
        let f = one_mesh(voxel);
        let t = f.read(7).unwrap();
        assert!(t.warnings.is_empty(), "{:?}", t.warnings);
        assert_eq!(t.result[[2, 0, 0, 0, 0]], 3.0);
        // Listed twice.
        let t = one_mesh(&format!("{voxel}{voxel}")).read(7).unwrap();
        assert!(t
            .warnings
            .iter()
            .any(|w| w.message.contains("listed more than once")));
        // A source that isn't the sum of its bins.
        let off = voxel.replace("  1  3.00000E+00", "  1  4.00000E+00");
        let t = one_mesh(&off).read(7).unwrap();
        assert!(
            t.warnings.iter().any(|w| w.message.contains("differs")),
            "{:?}",
            t.warnings
        );
        // A voxel out of the mesh, and cells cut short.
        let out = voxel.replace("  1  3.00000E+00", "  9  3.00000E+00");
        assert!(one_mesh(&out).read(7).is_err());
        let short =
            "  1  3.00000E+00  1.00000E+00  2\n  5  1.0  3.00000E+00\n 1.0E+00 2.0E+00\n0.1 0.1\n";
        assert!(one_mesh(short).read(7).is_err());
    }

    #[test]
    fn an_unreadable_mesh_id_is_failed_not_zero() {
        let f = Cdgs::from_bytes(
            &b"num_meshes 1\nmesh_id x\ncomment\nsource_data\nend_source_data\n"[..],
        );
        assert_eq!(f.failed().len(), 1);
        assert!(f.failed()[0].1.contains("mesh_id"), "{:?}", f.failed());
    }
}
