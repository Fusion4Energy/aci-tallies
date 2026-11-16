//! The columnar layouts (`out=col`, `out=cf`): one row per voxel.
//!
//! Rows are in the order energy (slowest), time, i, j, k (fastest),
//! confirmed against MCNP 6.2 and 6.3 output, Total bins included where
//! row-major order over the Total-padded shape puts them. That is the
//! memory order of [`Tally::result`](super::Tally::result): a row's index
//! within its block is its flat index, and the rows of one (energy, time)
//! bin are consecutive, so any range of them can be read on its own.

use std::ops::Range;

use crate::common::{float, scan, Warning};

use super::model::OutFormat;

/// The columns of a block, from its column header row.
///
/// `Energy` (if present) comes first and `Time` (if present) right after
/// it: `Time X Y Z ...`, `Energy Time X Y Z ...`. D1SUNED names discrete
/// bins there instead (`Nuclide bin X Y Z ...`, `Energy Nuclide X ...`,
/// `Nuclide Time X ...`, `Cell bin X ...`): one column each.
#[derive(Clone, Debug)]
pub(crate) struct Columns {
    /// Where Result is: after the Energy/Time and coordinate columns.
    result_col: usize,
    /// Columns a full row has.
    expected_cols: usize,
    /// Where Volume is, in `out=cf`.
    vol_col: Option<usize>,
}

impl Columns {
    pub fn from_header_row(line: &[u8]) -> Columns {
        // Every label before the first coordinate is one leading column,
        // but for the `bin` of `Nuclide bin` and `Cell bin`.
        let leading = scan::tokens(line)
            .take_while(|t| !matches!(*t, b"X" | b"R"))
            .filter(|t| *t != b"bin")
            .count()
            .min(2);
        let has_volume = memchr::memmem::find(line, b"Volume").is_some();
        let result_col = 3 + leading;
        Columns {
            result_col,
            expected_cols: result_col + 2 + 2 * usize::from(has_volume),
            vol_col: has_volume.then_some(result_col + 2),
        }
    }

    pub fn has_volume(&self) -> bool {
        self.vol_col.is_some()
    }

    pub fn out_format(&self) -> OutFormat {
        if self.has_volume() {
            OutFormat::Cf
        } else {
            OutFormat::Col
        }
    }
}

/// Where parsed rows go: rows `start..start + result.len()` of the block
/// (`NaN`-filled beforehand), and the volumes of rows `start..start +
/// volume.len()`.
pub(crate) struct RowsOut<'a> {
    pub start: u64,
    pub result: &'a mut [f64],
    pub rel_error: &'a mut [f64],
    pub volume: &'a mut [f64],
}

impl<'a> RowsOut<'a> {
    /// Cut into consecutive parts, one per entry of `lens` (their sum the
    /// whole length).
    pub fn split(self, lens: impl IntoIterator<Item = usize>) -> Vec<RowsOut<'a>> {
        let RowsOut {
            mut start,
            mut result,
            mut rel_error,
            mut volume,
        } = self;
        let mut out = Vec::new();
        for n in lens {
            let (r, rest_r) = result.split_at_mut(n);
            let (e, rest_e) = rel_error.split_at_mut(n);
            let (v, rest_v) = volume.split_at_mut(n.min(volume.len()));
            out.push(RowsOut {
                start,
                result: r,
                rel_error: e,
                volume: v,
            });
            start += n as u64;
            (result, rel_error, volume) = (rest_r, rest_e, rest_v);
        }
        out
    }
}

/// Parse the rows of `text` (whole lines) that land in `out`, given that
/// its first non-blank line is row `first_row` of the block (negative when
/// `text` starts before the block: those lines are only counted) and its
/// first line is line `first_line + 1` of the file. Blank lines are not
/// rows. A short row, or one whose Result or Rel Error can't be read,
/// still takes its slot (left `NaN`, with a warning), so the rows after it
/// stay where they belong. Allocation-free per row.
pub(crate) fn parse_rows(
    text: &[u8],
    first_row: i64,
    first_line: u64,
    cols: &Columns,
    out: RowsOut<'_>,
) -> Vec<Warning> {
    let mut warnings = Vec::new();
    let start = out.start as i64;
    let end = start + out.result.len() as i64;
    if end <= start {
        return warnings;
    }
    let mut row = first_row;
    for (n, line) in scan::lines(text).enumerate() {
        if scan::is_blank(line) {
            continue;
        }
        let r = row;
        row += 1;
        if r < start {
            continue;
        }
        if r >= end {
            break;
        }
        let k = (r - start) as usize;

        let (mut result, mut rel_error, mut volume) = (None, None, None);
        let mut n_tokens = 0usize;
        for (i, tok) in scan::tokens(line).enumerate() {
            n_tokens = i + 1;
            if i == cols.result_col {
                result = float::parse(tok);
            } else if i == cols.result_col + 1 {
                rel_error = float::parse(tok);
            } else if Some(i) == cols.vol_col {
                volume = float::parse(tok);
            }
        }
        let at = first_line + n as u64 + 1;
        let text = || String::from_utf8_lossy(line).trim().to_string();
        if n_tokens + 1 == cols.expected_cols {
            // A negative value can run into the coordinate before it
            // (`112.500-1.57688E-11`, when the coordinate fills its field).
            if let Some((r, e, v)) = unglued(line, cols) {
                out.result[k] = r;
                out.rel_error[k] = e;
                if let Some(slot) = out.volume.get_mut(k) {
                    *slot = v.unwrap_or(f64::NAN);
                }
                continue;
            }
        }
        if n_tokens < cols.expected_cols {
            warnings.push(Warning::at(
                at,
                format!(
                    "row {r}: expected {} columns, found {n_tokens}: {:?}",
                    cols.expected_cols,
                    text()
                ),
            ));
            continue;
        }
        let (Some(result), Some(rel_error)) = (result, rel_error) else {
            warnings.push(Warning::at(
                at,
                format!("row {r}: could not read Result/Rel Error: {:?}", text()),
            ));
            continue;
        };
        out.result[k] = result;
        out.rel_error[k] = rel_error;
        if let Some(v) = out.volume.get_mut(k) {
            *v = volume.unwrap_or(f64::NAN);
        }
    }
    warnings
}

/// Result, Rel Error and Volume of a row one token short because a
/// negative Result ran into the coordinate before it: the token where a
/// sign follows a digit and starts a number with an exponent is split
/// there. None when no token splits that way.
fn unglued(line: &[u8], cols: &Columns) -> Option<(f64, f64, Option<f64>)> {
    let mut tokens: Vec<&[u8]> = Vec::with_capacity(cols.expected_cols);
    let mut split = false;
    for tok in scan::tokens(line) {
        let cut = (!split && tokens.len() < cols.result_col)
            .then(|| glued_at(tok))
            .flatten();
        match cut {
            Some(p) => {
                tokens.push(&tok[..p]);
                tokens.push(&tok[p..]);
                split = true;
            }
            None => tokens.push(tok),
        }
    }
    if !split || tokens.len() != cols.expected_cols {
        return None;
    }
    let result = float::parse(tokens[cols.result_col])?;
    let rel_error = float::parse(tokens[cols.result_col + 1])?;
    let volume = cols.vol_col.and_then(|c| float::parse(tokens[c]));
    Some((result, rel_error, volume))
}

/// Where `tok` splits into a fixed-point number and a number with an
/// `E` exponent, `112.500-1.57688E-11` at its `-`.
fn glued_at(tok: &[u8]) -> Option<usize> {
    (1..tok.len()).find(|&p| {
        matches!(tok[p], b'-' | b'+')
            && tok[p - 1].is_ascii_digit()
            && !tok[..p].iter().any(|c| matches!(c, b'E' | b'e'))
            && tok[p..].iter().any(|c| matches!(c, b'E' | b'e'))
            && float::parse(&tok[..p]).is_some()
            && float::parse(&tok[p..]).is_some()
    })
}

/// The warnings for a block of `seen` rows where `expected` were due.
pub(crate) fn row_count_warning(seen: u64, expected: u64) -> Option<Warning> {
    if seen < expected {
        Some(Warning::new(format!(
            "expected {expected} data rows, found only {seen}; the rest are NaN"
        )))
    } else if seen > expected {
        Some(Warning::new(format!(
            "found {seen} data rows but expected {expected}; the extra ones are ignored"
        )))
    } else {
        None
    }
}

/// The rows of energy index `e` and time index `t`, in a block of `n_time`
/// time indices and `voxels` voxels.
pub(crate) fn bin_rows(e: usize, t: usize, n_time: usize, voxels: usize) -> Range<u64> {
    let start = ((e * n_time + t) * voxels) as u64;
    start..start + voxels as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Source;
    use crate::meshtal::Reader;
    use ndarray::Array5;
    use proptest::prelude::*;

    /// A one-tally file: an `n1 x n2 x n3` mesh without energy or time
    /// bins, and these data rows.
    fn file(n1: usize, n2: usize, n3: usize, rows: &[String]) -> String {
        let b = |n: usize| {
            (0..=n)
                .map(|i| format!("{i}.00"))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let mut text = format!(
            " Mesh Tally Number 4\n neutron  mesh tally.\n\n Tally bin boundaries:\n    \
             X direction: {}\n    Y direction: {}\n    Z direction: {}\n    \
             Energy bin boundaries: 0.00E+00 1.00E+36\n\n   \
             X         Y         Z     Result     Rel Error\n",
            b(n1),
            b(n2),
            b(n3)
        );
        for r in rows {
            text.push_str(r);
            text.push('\n');
        }
        text
    }

    fn next_unit(state: &mut u64) -> f64 {
        // splitmix64, so that shrinking works.
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Values as MCNP prints them, and their rows.
    fn synthetic(n: [usize; 3], seed: u64) -> (Array5<f64>, Array5<f64>, Vec<String>) {
        let shape = (1, 1, n[0], n[1], n[2]);
        let mut state = seed;
        let mut rows = Vec::new();
        let mut result = Array5::zeros(shape);
        let mut rel_error = Array5::zeros(shape);
        for (r, e) in result.iter_mut().zip(rel_error.iter_mut()) {
            let v = format!("{:.5E}", (next_unit(&mut state) - 0.5) * 2.0e5);
            let err = format!("{:.4}", next_unit(&mut state));
            *r = v.parse().unwrap();
            *e = err.parse().unwrap();
            rows.push(format!("  0.0 0.0 0.0 {v} {err}"));
        }
        (result, rel_error, rows)
    }

    fn read(text: &str, chunk: usize) -> crate::meshtal::Tally {
        let reader = Reader::build(Source::bytes(text.as_bytes().to_vec()), None, chunk).unwrap();
        reader.read(4).unwrap()
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 64, ..ProptestConfig::default() })]

        #[test]
        fn round_trips(n1 in 1usize..=4, n2 in 1usize..=4, n3 in 1usize..=4,
                       chunk in 1usize..=200, seed in any::<u64>()) {
            let (result, rel_error, rows) = synthetic([n1, n2, n3], seed);
            let t = read(&file(n1, n2, n3, &rows), chunk);
            prop_assert!(t.warnings.is_empty(), "{:?}", t.warnings);
            prop_assert!(t.volume.is_none());
            prop_assert_eq!(t.result, result);
            prop_assert_eq!(t.rel_error, rel_error);
        }

        #[test]
        fn damage_stays_local(n1 in 1usize..=4, n2 in 1usize..=4, n3 in 1usize..=4,
                              chunk in 1usize..=120, seed in any::<u64>(),
                              at in 0usize..64, kind in 0u8..3) {
            let (result, rel_error, mut rows) = synthetic([n1, n2, n3], seed);
            let at = at % rows.len();
            match kind {
                0 => {
                    let short: Vec<&str> = rows[at].split_whitespace().take(4).collect();
                    rows[at] = short.join(" ");
                }
                1 => {
                    let mut t: Vec<String> = rows[at].split_whitespace().map(str::to_string).collect();
                    t[3] = "notanumber".to_string();
                    rows[at] = t.join(" ");
                }
                // A blank line is not a row, wherever it falls.
                _ => rows.insert(at, String::new()),
            }
            let t = read(&file(n1, n2, n3, &rows), chunk);
            for (flat, ((got, want), (got_e, want_e))) in t.result.iter().zip(&result)
                .zip(t.rel_error.iter().zip(&rel_error)).enumerate() {
                if kind != 2 && flat == at {
                    prop_assert!(got.is_nan() && got_e.is_nan());
                } else {
                    prop_assert_eq!(got, want);
                    prop_assert_eq!(got_e, want_e);
                }
            }
            if kind == 2 {
                prop_assert!(t.warnings.is_empty(), "{:?}", t.warnings);
            } else {
                prop_assert_eq!(t.warnings.len(), 1, "{:?}", t.warnings);
                let prefix = format!("row {at}:");
                prop_assert!(t.warnings[0].message.starts_with(&prefix));
                // 10 header lines before the first row.
                prop_assert_eq!(t.warnings[0].line, Some(at as u64 + 11));
            }
        }
    }

    #[test]
    fn too_few_rows() {
        let (_, _, rows) = synthetic([2, 2, 2], 42);
        for chunk in [1, 50, 1 << 20] {
            let t = read(&file(2, 2, 2, &rows[..5]), chunk);
            assert_eq!(t.warnings.len(), 1, "{:?}", t.warnings);
            assert!(t.warnings[0]
                .message
                .contains("expected 8 data rows, found only 5"));
            assert!(t.result.iter().skip(5).all(|v| v.is_nan()));
        }
    }

    #[test]
    fn too_many_rows() {
        let (result, _, mut rows) = synthetic([2, 1, 2], 7);
        rows.push("0.0 0.0 0.0 notanumber 0.02".to_string());
        rows.push("0.0 0.0 0.0 1.0 0.02".to_string());
        for chunk in [1, 60, 1 << 20] {
            let t = read(&file(2, 1, 2, &rows), chunk);
            // Only the count: the extra rows are not read.
            assert_eq!(t.warnings.len(), 1, "{:?}", t.warnings);
            assert!(t.warnings[0]
                .message
                .contains("found 6 data rows but expected 4"));
            assert_eq!(t.result, result);
        }
    }

    #[test]
    fn columns() {
        let c = Columns::from_header_row(
            b"  Energy   Time   X   Y   Z  Result  Rel Error  Volume  Rslt * Vol",
        );
        assert_eq!((c.result_col, c.expected_cols, c.vol_col), (5, 9, Some(7)));
        assert_eq!(c.out_format(), OutFormat::Cf);
        let c = Columns::from_header_row(b"   Time   R   Z   Th  Result  Rel Error");
        assert_eq!((c.result_col, c.expected_cols, c.vol_col), (4, 6, None));
        // D1SUNED discrete bins.
        for (row, result_col) in [
            (
                &b"   Nuclide bin    X         Y         Z     Result     Rel Error"[..],
                4,
            ),
            (
                b"   Cell bin       X         Y         Z     Result     Rel Error",
                4,
            ),
            (
                b"   Energy      Nuclide         X    Y    Z     Result     Rel Error",
                5,
            ),
            (
                b"   Nuclide      Time         X    Y    Z     Result     Rel Error",
                5,
            ),
        ] {
            let c = Columns::from_header_row(row);
            assert_eq!(
                (c.result_col, c.expected_cols),
                (result_col, result_col + 2)
            );
        }
    }

    #[test]
    fn negative_value_run_into_its_coordinate() {
        let mut rows = vec!["0.0 0.0 0.0 1.0 0.5".to_string(); 4];
        rows[2] = "1.500   1.500    12.500-1.57688E-11 2.61718E+00".to_string();
        let t = read(&file(2, 2, 1, &rows), 1 << 20);
        assert!(t.warnings.is_empty(), "{:?}", t.warnings);
        assert_eq!(t.result[[0, 0, 1, 0, 0]], -1.57688e-11);
        assert_eq!(t.rel_error[[0, 0, 1, 0, 0]], 2.61718);
        assert_eq!(glued_at(b"12.500-1.57688E-11"), Some(6));
        assert_eq!(glued_at(b"1.57688E-11"), None);
        assert_eq!(glued_at(b"1.5-104"), None);
    }
}
