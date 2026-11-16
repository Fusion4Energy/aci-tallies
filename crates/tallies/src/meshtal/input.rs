//! Replacing a section's rounded header geometry with the exact FMESH
//! definition from the MCNP input.

use std::collections::HashMap;

use aci_tallies_geom::mcnp_input::{FmeshCard, FmeshError};

use super::header::{HeaderLabels, SectionHeader};
use crate::common::{Input, LabelKind, Warning};
use crate::meshtal::model::Geometry;

/// The FMESH cards of an MCNP input, by tally number.
pub(crate) struct InputMeshes {
    cards: HashMap<u32, Result<FmeshCard, FmeshError>>,
}

impl InputMeshes {
    /// Every FMESH card of `input`. A tally number given twice keeps its
    /// first card (MCNP rejects such an input anyway).
    pub fn of(input: &Input) -> InputMeshes {
        let mut cards = HashMap::new();
        for card in aci_tallies_geom::mcnp_input::fmeshes_of(input.model()) {
            let id = match &card {
                Ok(c) => c.id,
                Err(e) => e.id,
            };
            cards.entry(id).or_insert(card);
        }
        InputMeshes { cards }
    }

    pub fn get(&self, id: u32) -> Option<&Result<FmeshCard, FmeshError>> {
        self.cards.get(&id)
    }
}

/// Bounds printed with `%.2f` (theta `%.3f`) are within half a unit of the
/// last decimal; a little slack covers the float formatting itself.
const FIXED_TOL: f64 = 0.0051;
const THETA_TOL: f64 = 0.00051;
/// Origin, axes and energy/time bounds are printed with `%.2E`: three
/// significant digits.
const SIG3_TOL: f64 = 0.0051;

/// Check that `card` is the mesh the header describes, at the header's
/// print precision, and if so move its exact geometry and energy/time
/// bounds into `hdr`. On a mismatch `hdr` is left as it was.
fn apply(hdr: &mut SectionHeader, card: &FmeshCard) -> Result<(), String> {
    let fixed = |tol: f64| move |a: f64, b: f64| (a - b).abs() <= tol;
    let sig3 = |a: f64, b: f64| (a - b).abs() <= SIG3_TOL * a.abs().max(b.abs()) + 1e-12;
    let same = |name: &str, header: &[f64], exact: &[f64], close: &dyn Fn(f64, f64) -> bool| {
        if header.len() != exact.len() {
            return Err(format!(
                "{name} has {} bounds in the header, {} in the input",
                header.len(),
                exact.len()
            ));
        }
        match header.iter().zip(exact).position(|(&a, &b)| !close(a, b)) {
            Some(n) => Err(format!(
                "{name} bound {n} is {} in the header, {} in the input",
                header[n], exact[n]
            )),
            None => Ok(()),
        }
    };

    match (&hdr.geometry, &card.local) {
        (
            Geometry::Rectangular {
                x_bounds: x,
                y_bounds: y,
                z_bounds: z,
                ..
            },
            Geometry::Rectangular {
                x_bounds: ex,
                y_bounds: ey,
                z_bounds: ez,
                ..
            },
        ) => {
            same("x", x, ex, &fixed(FIXED_TOL))?;
            same("y", y, ey, &fixed(FIXED_TOL))?;
            same("z", z, ez, &fixed(FIXED_TOL))?;
        }
        (
            Geometry::Cylindrical {
                r_bounds: r,
                z_bounds: z,
                theta_bounds: t,
                origin,
                axis,
                vec,
            },
            Geometry::Cylindrical {
                r_bounds: er,
                z_bounds: ez,
                theta_bounds: et,
                origin: eo,
                axis: ea,
                vec: ev,
            },
        ) => {
            same("r", r, er, &fixed(FIXED_TOL))?;
            same("z", z, ez, &fixed(FIXED_TOL))?;
            same("theta", t, et, &fixed(THETA_TOL))?;
            same("origin", origin, eo, &sig3)?;
            same("axis", axis, ea, &sig3)?;
            same("vec", vec, ev, &sig3)?;
        }
        _ => return Err("the header and the input disagree on GEOM".into()),
    }
    // A meshinfo header has no energy bounds to compare, nor has a
    // D1SUNED one whose bins are nuclides (or cells, decay times) instead.
    let energy_bounded = hdr.energy_labels.is_none() && !hdr.energy_bounds.is_empty();
    let time_bounded = hdr.time_labels.is_none();
    if let Some(e) = card.energy_bounds.as_ref().filter(|_| energy_bounded) {
        same("energy", &hdr.energy_bounds, e, &sig3)?;
    }
    if let Some(t) = card.time_bounds.as_ref().filter(|_| time_bounded) {
        same("time", &hdr.time_bounds, t, &sig3)?;
    }
    let energy_labels = match &hdr.energy_labels {
        Some(l) if l.overflowed => Some(recovered(
            l,
            card.energy_tag.as_deref(),
            card.energy_values.as_deref(),
            "ETAG/EMESH",
        )?),
        _ => None,
    };
    let time_labels = match &hdr.time_labels {
        Some(l) if l.overflowed => Some(recovered(
            l,
            card.time_tag.as_deref(),
            card.time_values.as_deref(),
            "TTAG/TMESH",
        )?),
        _ => None,
    };

    hdr.geometry = card.geometry().map_err(|e| e.to_string())?;
    if let Some(e) = card.energy_bounds.as_ref().filter(|_| energy_bounded) {
        hdr.energy_bounds = e.clone();
    }
    if let Some(t) = card.time_bounds.as_ref().filter(|_| time_bounded) {
        hdr.time_bounds = t.clone();
    }
    if let Some(l) = energy_labels {
        hdr.energy_labels = Some(l);
    }
    if let Some(l) = time_labels {
        hdr.time_labels = Some(l);
    }
    Ok(())
}

/// The labels D1SUNED meant to print for an overflowed list, from the
/// FMESH card's tag and bins: `0.001` (every other nuclide or cell), then
/// `<bin>.<n>` for the n-th bin from 2 on. The labels the header did
/// print must agree.
fn recovered(
    header: &HeaderLabels,
    tag: Option<&str>,
    values: Option<&[f64]>,
    card_words: &str,
) -> Result<HeaderLabels, String> {
    let expected = match header.kind {
        LabelKind::Nuclide => &["nuc"][..],
        LabelKind::Parent => &["parent", "nuc"][..],
        LabelKind::Cell => &["cel", "cell"][..],
        LabelKind::DecayTime | LabelKind::Line => &[][..],
    };
    let (Some(tag), Some(values)) = (tag, values) else {
        return Err(format!(
            "the header's {} list overflowed and the card has no {card_words}",
            header.kind.name()
        ));
    };
    if !expected.contains(&tag) {
        return Err(format!(
            "the header's {} list overflowed but the card's tag is {tag:?}",
            header.kind.name()
        ));
    }
    let mut labels = vec![0.001];
    for (n, v) in values.iter().enumerate() {
        let label = format!("{}.{:03}", *v as i64, n + 2);
        labels.push(label.parse().map_err(|_| format!("bad label {label}"))?);
    }
    let printed = &header.values;
    if printed.len() > labels.len()
        || printed
            .iter()
            .zip(&labels)
            // Absolute: the labels are ids, which a relative tolerance blurs.
            .any(|(a, b)| (a - b).abs() > 1e-6)
    {
        return Err(format!(
            "the header's {} list {printed:?} does not start like the card's {labels:?}",
            header.kind.name()
        ));
    }
    Ok(HeaderLabels {
        kind: header.kind,
        values: labels,
        overflowed: false,
    })
}

/// A section header with the input's exact geometry applied (when given
/// and it matches), and what to warn about.
pub(crate) struct Prepared {
    pub header: SectionHeader,
    pub geometry_from_input: bool,
    pub energy_weighted: Option<bool>,
    pub warnings: Vec<Warning>,
}

pub(crate) fn prepare(hdr: SectionHeader, input: Option<&InputMeshes>) -> Prepared {
    prepare_from(hdr, input, "meshtal")
}

/// [`prepare`] for a section of a `file` file (named in the warnings).
pub(crate) fn prepare_from(
    mut hdr: SectionHeader,
    input: Option<&InputMeshes>,
    file: &str,
) -> Prepared {
    let mut warnings = Vec::new();
    let mut geometry_from_input = false;
    let mut energy_weighted = None;
    if let Some(input) = input {
        let id = hdr.id;
        let problem = match input.get(id) {
            None => Some(format!("the input has no FMESH{id} card")),
            Some(Err(e)) => Some(format!("the input's {e}")),
            Some(Ok(card)) => match apply(&mut hdr, card) {
                Ok(()) => {
                    geometry_from_input = true;
                    energy_weighted = Some(card.energy_weighted);
                    None
                }
                Err(why) => Some(format!(
                    "the input's FMESH{id} card (line {}) is not this mesh, so the \
                     input is probably not the one this file came from ({why})",
                    card.line
                )),
            },
        };
        if let Some(problem) = problem {
            warnings.push(Warning::new(format!(
                "{problem}; geometry is from the {file} header"
            )));
        }
    }
    // MCNP prints bounds with 2 decimals (3 for theta), so bins finer than
    // that collapse to zero width in the header. The values are still fine;
    // what is derived from the geometry (volumes, VTK) is not. (Geometry
    // from the input is checked when it's read.)
    if let Err(e) = hdr.geometry.validate() {
        warnings.push(Warning::new(format!(
            "invalid geometry from the {file} header: {e} (bins finer than the \
             printed precision?); geometry-derived values will be wrong"
        )));
    }
    Prepared {
        header: hdr,
        geometry_from_input,
        energy_weighted,
        warnings,
    }
}
