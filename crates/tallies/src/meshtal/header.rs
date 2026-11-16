use crate::common::{scan, Axis, LabelKind};
use crate::meshtal::model::Geometry;

#[derive(Debug, Clone, thiserror::Error)]
pub(crate) enum HeaderError {
    #[error("could not find a \"Mesh Tally Number\" line")]
    MissingTallyNumber,
    #[error("could not parse tally id from line: {0:?}")]
    BadTallyId(String),
    #[error("unexpected end of section while looking for {0}")]
    UnexpectedEnd(&'static str),
    #[error("could not parse bin boundaries from line: {0:?}")]
    BadBinBoundaries(String),
    #[error("could not parse origin/axis/vec from line: {0:?}")]
    BadCylindricalOrigin(String),
    #[error("more than two energy/time bin lines: {0:?}")]
    TooManyBinLines(String),
}

/// The discrete bins a D1SUNED mesh prints in place of energy or time
/// bounds (`Daughter Nuclides list:`, `Decay times:`, ...). A label is
/// None where the file overflowed its field (`**********`); `overflowed`
/// says that happened, after which even how many labels there are is
/// unknown.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct HeaderLabels {
    pub kind: LabelKind,
    pub values: Vec<f64>,
    pub overflowed: bool,
}

impl HeaderLabels {
    fn parse(kind: LabelKind, line: &str) -> Option<HeaderLabels> {
        let after = line.split_once(':')?.1;
        // D1SUNED writes each label in 10 columns, so an overflowed one
        // runs into its neighbours: take what reads before the first `*`.
        let (good, overflowed) = match after.find('*') {
            Some(p) => (&after[..p], true),
            None => (after, false),
        };
        let mut values = Vec::new();
        for tok in good.split_whitespace() {
            match tok.parse::<f64>() {
                Ok(v) => values.push(v),
                // A token cut by the first `*`.
                Err(_) if overflowed => break,
                Err(_) => return None,
            }
        }
        if values.is_empty() && !overflowed {
            return None;
        }
        Some(HeaderLabels {
            kind,
            values,
            overflowed,
        })
    }

    /// The label line's kind, from its text.
    fn kind_of(line: &str) -> Option<LabelKind> {
        let t = line.trim_start();
        if t.starts_with("Daughter Nuclides list:") {
            Some(LabelKind::Nuclide)
        } else if t.starts_with("Parent Nuclides list:") {
            Some(LabelKind::Parent)
        } else if t.starts_with("Cells list:") {
            Some(LabelKind::Cell)
        } else if t.starts_with("Decay times:") {
            Some(LabelKind::DecayTime)
        } else {
            None
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct SectionHeader {
    pub id: u32,
    pub particle: String,
    pub comment: Option<String>,
    pub geometry: Geometry,
    pub energy_bounds: Vec<f64>,
    pub time_bounds: Vec<f64>,
    /// D1SUNED discrete bins in place of the energy bounds (which are then
    /// left at their default, one bin).
    pub energy_labels: Option<HeaderLabels>,
    /// The same, in place of the time bounds.
    pub time_labels: Option<HeaderLabels>,
}

impl SectionHeader {
    /// Number of real energy groups implied by `energy_bounds`. A tally
    /// with no `EMESH` card still has one implicit group.
    pub(crate) fn n_energy_from_bounds(&self) -> usize {
        match &self.energy_labels {
            Some(l) => l.values.len().max(1),
            None => self.energy_bounds.len().saturating_sub(1).max(1),
        }
    }

    /// Number of real time bins implied by `time_bounds`. A tally with no
    /// `TMESH` card still has one implicit bin.
    pub(crate) fn n_time_from_bounds(&self) -> usize {
        match &self.time_labels {
            Some(l) => l.values.len().max(1),
            None => self.time_bounds.len().saturating_sub(1).max(1),
        }
    }

    /// The labels, if any, that overflowed and are still unknown.
    pub(crate) fn overflowed_labels(&self) -> Option<&HeaderLabels> {
        [&self.energy_labels, &self.time_labels]
            .into_iter()
            .flatten()
            .find(|l| l.overflowed)
    }

    /// The energy axis, with MCNP's Total bin when there is more than one
    /// bin (`None`: decide from the data, as a matrix layout does).
    pub(crate) fn energy_axis(&self, total: Option<bool>) -> Axis {
        axis_of(&self.energy_bounds, self.energy_labels.as_ref(), total)
    }

    /// The time axis, the same way.
    pub(crate) fn time_axis(&self, total: Option<bool>) -> Axis {
        axis_of(&self.time_bounds, self.time_labels.as_ref(), total)
    }
}

fn axis_of(bounds: &[f64], labels: Option<&HeaderLabels>, total: Option<bool>) -> Axis {
    match (labels, total) {
        (Some(l), None) => Axis::labelled(l.kind, l.values.clone()),
        (Some(l), Some(t)) => Axis::labelled_with_total(l.kind, l.values.clone(), t),
        (None, None) => Axis::mcnp(bounds.to_vec()),
        (None, Some(t)) => Axis::bounded(bounds.to_vec(), t),
    }
}

/// Parse the header of one `Mesh Tally Number` section.
///
/// `lines` is every line belonging to this section (from the `Mesh Tally
/// Number` line up to, but not including, the next section's). Returns the
/// parsed header plus the index into `lines` where the data block starts
/// (the first line after the header that is not blank).
pub(crate) fn parse_header(lines: &[&str]) -> Result<(SectionHeader, usize), HeaderError> {
    let mut i = 0;

    let id_line = lines.first().ok_or(HeaderError::MissingTallyNumber)?;
    let id = parse_tally_id(id_line)?;
    i += 1;

    // Optional FC comment line: present iff the next line is NOT the
    // "<particle>  mesh tally." line.
    let mut comment = None;
    if !is_particle_line(get(lines, i)?) {
        comment = Some(lines[i].trim().to_string());
        i += 1;
    }
    let particle = parse_particle(get(lines, i)?);
    i += 1;

    i = skip_blank(lines, i);
    i = find_line_containing(lines, i, "Tally bin boundaries:")?;
    i += 1;
    i = skip_blank(lines, i);

    let geometry = if get(lines, i)?.contains("origin at") {
        let (origin, axis, vec) = parse_cylindrical_origin(lines[i])?;
        i += 1;
        let r_bounds = parse_direction_line(get(lines, i)?, "R direction:")?;
        i += 1;
        let z_bounds = parse_direction_line(get(lines, i)?, "Z direction:")?;
        i += 1;
        let theta_bounds = parse_direction_line(get(lines, i)?, "Theta direction (revolutions):")?;
        i += 1;
        Geometry::Cylindrical {
            r_bounds,
            z_bounds,
            theta_bounds,
            origin,
            axis,
            vec,
        }
    } else {
        let x_bounds = parse_direction_line(get(lines, i)?, "X direction:")?;
        i += 1;
        let y_bounds = parse_direction_line(get(lines, i)?, "Y direction:")?;
        i += 1;
        let z_bounds = parse_direction_line(get(lines, i)?, "Z direction:")?;
        i += 1;
        Geometry::Rectangular {
            x_bounds,
            y_bounds,
            z_bounds,
            transform: None,
        }
    };

    // The `Time bin boundaries:` line, when present (only when a `TMESH`
    // card was used), always precedes `Energy bin boundaries:` — confirmed
    // against `examples/meshtal` tallies 84/94. D1SUNED prints discrete
    // bins (`ETAG`/`TTAG`: nuclides, cells, decay times) in place of
    // either, in the same order: of two bin lines the first is time, the
    // second energy; a lone one is energy.
    let mut bin_lines: Vec<(BinLine, &str)> = Vec::new();
    loop {
        let line = get(lines, i)?;
        let kind = if line.contains("Time bin boundaries:") {
            BinLine::TimeBounds
        } else if line.contains("Energy bin boundaries:") {
            BinLine::EnergyBounds
        } else if let Some(k) = HeaderLabels::kind_of(line) {
            BinLine::Labels(k)
        } else {
            break;
        };
        bin_lines.push((kind, line));
        i += 1;
        if kind == BinLine::EnergyBounds {
            break;
        }
    }
    let (time_line, energy_line) = match bin_lines.as_slice() {
        [] => return Err(HeaderError::BadBinBoundaries(get(lines, i)?.to_string())),
        [e] => (None, *e),
        [t, e] => (Some(*t), *e),
        [.., last] => return Err(HeaderError::TooManyBinLines(last.1.to_string())),
    };
    let bad = |line: &str| HeaderError::BadBinBoundaries(line.to_string());
    let (mut time_bounds, mut time_labels) = (vec![0.0, 1e36], None);
    if let Some((kind, line)) = time_line {
        match kind {
            BinLine::TimeBounds => {
                time_bounds = parse_floats_after_colon(line).ok_or_else(|| bad(line))?
            }
            BinLine::Labels(k) => {
                time_labels = Some(HeaderLabels::parse(k, line).ok_or_else(|| bad(line))?)
            }
            BinLine::EnergyBounds => return Err(bad(line)),
        }
    }
    let (mut energy_bounds, mut energy_labels) = (vec![0.0, 1e36], None);
    match energy_line {
        (BinLine::EnergyBounds, line) => {
            energy_bounds = parse_floats_after_colon(line).ok_or_else(|| bad(line))?
        }
        (BinLine::Labels(k), line) => {
            energy_labels = Some(HeaderLabels::parse(k, line).ok_or_else(|| bad(line))?)
        }
        (BinLine::TimeBounds, line) => return Err(bad(line)),
    }

    i = skip_blank(lines, i);

    let header = SectionHeader {
        id,
        particle,
        comment,
        geometry,
        energy_bounds,
        time_bounds,
        energy_labels,
        time_labels,
    };
    Ok((header, i))
}

/// A line of bins after the mesh bounds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BinLine {
    TimeBounds,
    EnergyBounds,
    Labels(LabelKind),
}

fn get<'a>(lines: &[&'a str], i: usize) -> Result<&'a str, HeaderError> {
    lines
        .get(i)
        .copied()
        .ok_or(HeaderError::UnexpectedEnd("more header lines"))
}

fn skip_blank(lines: &[&str], mut i: usize) -> usize {
    while i < lines.len() && lines[i].trim().is_empty() {
        i += 1;
    }
    i
}

fn find_line_containing(lines: &[&str], start: usize, needle: &str) -> Result<usize, HeaderError> {
    for (offset, line) in lines[start..].iter().enumerate() {
        if line.contains(needle) {
            return Ok(start + offset);
        }
    }
    Err(HeaderError::UnexpectedEnd("Tally bin boundaries:"))
}

fn parse_tally_id(line: &str) -> Result<u32, HeaderError> {
    let after = line
        .find("Mesh Tally Number")
        .map(|idx| &line[idx + "Mesh Tally Number".len()..])
        .ok_or(HeaderError::MissingTallyNumber)?;
    after
        .split_whitespace()
        .next()
        .and_then(|tok| tok.parse::<u32>().ok())
        .ok_or_else(|| HeaderError::BadTallyId(line.to_string()))
}

fn is_particle_line(line: &str) -> bool {
    line.contains("mesh tally.")
}

fn parse_particle(line: &str) -> String {
    line.split("mesh tally.")
        .next()
        .unwrap_or_default()
        .trim()
        .to_string()
}

fn parse_direction_line(line: &str, label: &str) -> Result<Vec<f64>, HeaderError> {
    if !line.contains(label) {
        return Err(HeaderError::BadBinBoundaries(line.to_string()));
    }
    parse_floats_after_colon(line).ok_or_else(|| HeaderError::BadBinBoundaries(line.to_string()))
}

fn parse_floats_after_colon(line: &str) -> Option<Vec<f64>> {
    let after = line.split_once(':')?.1;
    let values: Option<Vec<f64>> = after
        .split_whitespace()
        .map(|tok| tok.parse::<f64>().ok())
        .collect();
    values.filter(|v| !v.is_empty())
}

type Triple = [f64; 3];

/// Parse the `origin at X Y Z axis in X Y Z direction, VEC direction X Y Z`
/// line that precedes cylindrical mesh bin boundaries.
fn parse_cylindrical_origin(line: &str) -> Result<(Triple, Triple, Triple), HeaderError> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    let err = || HeaderError::BadCylindricalOrigin(line.to_string());

    let at = tokens.iter().position(|&t| t == "at").ok_or_else(err)?;
    let origin = parse_triple(&tokens, at + 1).ok_or_else(err)?;

    let axis_in = tokens
        .iter()
        .position(|&t| t == "in")
        .filter(|&p| p > at)
        .ok_or_else(err)?;
    let axis = parse_triple(&tokens, axis_in + 1).ok_or_else(err)?;

    let vec_kw = tokens.iter().position(|&t| t == "VEC").ok_or_else(err)?;
    // tokens[vec_kw] == "VEC", tokens[vec_kw + 1] == "direction", then 3 floats.
    let vec = parse_triple(&tokens, vec_kw + 2).ok_or_else(err)?;

    Ok((origin, axis, vec))
}

fn parse_triple(tokens: &[&str], start: usize) -> Option<[f64; 3]> {
    let a = tokens.get(start)?.parse().ok()?;
    let b = tokens.get(start + 1)?.parse().ok()?;
    let c = tokens.get(start + 2)?.parse().ok()?;
    Some([a, b, c])
}

/// A section's header and where its data block starts (bytes from the
/// section's start; its end when there is no data).
pub(crate) fn read_header(section: &str) -> Result<(SectionHeader, usize), HeaderError> {
    // Headers take a dozen lines or so: look at the first few, and at all
    // of them only when that is not enough.
    let mut take = 64;
    loop {
        let lines: Vec<&str> = section.lines().take(take).collect();
        let complete = lines.len() < take;
        let r = parse_header(&lines);
        let more = !complete
            && match &r {
                Ok((_, i)) => *i >= lines.len(),
                Err(e) => matches!(e, HeaderError::UnexpectedEnd(_)),
            };
        if more {
            take = usize::MAX;
            continue;
        }
        return r.map(|(h, i)| (h, scan::line_start(section.as_bytes(), i)));
    }
}

/// What the lines before the first tally say about the run.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FileHeader {
    /// The code that wrote the file, as it names itself (`mcnp6`).
    pub code: String,
    pub version: String,
    /// When the problem ran, as printed after `probid =`.
    pub probid: String,
    /// The problem's title (the input's first line).
    pub title: String,
    /// The number of histories the tallies are normalized to.
    pub histories: Option<f64>,
}

/// Read the lines before the first section: a first line like `mcnp6
/// version 6 ld=01/26/23 probid = 10/15/26 14:36:50`, the title, and
/// `Number of histories used for normalizing tallies = 100.00`. Anything
/// missing is left empty.
pub(crate) fn read_file_header(text: &str) -> FileHeader {
    let mut lines = text.lines();
    let mut h = FileHeader::default();
    if let Some(first) = lines.next() {
        let mut words = first.split_whitespace();
        h.code = words.next().unwrap_or_default().to_string();
        let words: Vec<&str> = words.collect();
        if let Some(p) = words.iter().position(|w| *w == "version") {
            h.version = words.get(p + 1).copied().unwrap_or_default().to_string();
        }
        if let Some((_, probid)) = first.split_once("probid =") {
            h.probid = probid.trim().to_string();
        }
    }
    h.title = lines.next().unwrap_or_default().trim().to_string();
    for line in lines {
        if let Some((_, n)) = line.split_once("Number of histories used for normalizing tallies =")
        {
            h.histories = n.trim().parse().ok();
        }
    }
    h
}
