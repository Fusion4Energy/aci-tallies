//! D1SUNED's discrete bins in meshtal files (`tests/fixtures/d1suned_bins`,
//! see its README): values encode their bins, `1000 e + 100 t + voxel + 1`.

use aci_tallies::common::{Bins, Input, LabelKind};
use aci_tallies::meshtal::Reader;

mod common;

use common::fixture;

const TOTAL: f64 = 9.0;

fn value(e: f64, t: f64, voxel: usize) -> f64 {
    1000.0 * e + 100.0 * t + voxel as f64 + 1.0
}

#[test]
fn nuclides_cells_and_decay_times() {
    let r = Reader::open(fixture("d1suned_bins/bins.msht")).unwrap();
    let failed: Vec<u32> = r.failed().map(|(id, _)| id).collect();
    assert_eq!(failed, [94], "only the overflowed list needs the input");

    let t = r.read(44).unwrap();
    assert!(t.warnings.is_empty(), "{:?}", t.warnings);
    assert_eq!(
        t.energy.bins,
        Bins::Labels {
            kind: LabelKind::Nuclide,
            values: vec![0.001, 26059.002, 27060.003]
        }
    );
    assert!(t.energy.total);
    for (e, label) in [1.0, 2.0, 3.0, TOTAL].into_iter().enumerate() {
        for v in 0..2 {
            assert_eq!(t.result[[e, 0, v, 0, 0]], value(label, 1.0, v));
        }
    }

    let t = r.read(64).unwrap();
    assert_eq!(t.time.label_kind(), Some(LabelKind::DecayTime));
    assert!(!t.time.total, "decay times have no Total");
    assert_eq!(t.energy.label_kind(), Some(LabelKind::Nuclide));
    assert_eq!(t.result.shape(), [3, 2, 2, 1, 1]);
    assert_eq!(t.result[[2, 1, 1, 0, 0]], value(TOTAL, 2.0, 1));

    for (id, kind) in [(74, LabelKind::Parent), (84, LabelKind::Cell)] {
        assert_eq!(r.read(id).unwrap().energy.label_kind(), Some(kind));
    }
}

#[test]
fn an_overflowed_list_from_the_input() {
    let input = Input::parse(&std::fs::read_to_string(fixture("d1suned_bins/input")).unwrap());
    let r = Reader::open_with_input(fixture("d1suned_bins/bins.msht"), &input).unwrap();
    let t = r.read(94).unwrap();
    let labels = t.energy.labels().unwrap();
    assert_eq!(labels.len(), 3);
    assert!((labels[1] - 41092900.002).abs() < 1e-6, "{labels:?}");
    assert_eq!(t.result[[3, 0, 1, 0, 0]], value(TOTAL, 1.0, 1));
}

#[test]
fn a_lone_energy_bin_printed_as_total() {
    let r = Reader::open(fixture("d1suned_bins/bins.msht")).unwrap();
    let t = r.read(104).unwrap();
    assert!(t.warnings.is_empty(), "{:?}", t.warnings);
    assert!(!t.energy.total);
    assert_eq!(t.result.iter().copied().collect::<Vec<_>>(), [1.0, 2.0]);
}
