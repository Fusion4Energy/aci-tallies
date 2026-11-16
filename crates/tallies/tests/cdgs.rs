//! CDGS files (`tests/fixtures/cdgs`, see its README).

use aci_tallies::cdgs::Cdgs;
use aci_tallies::meshtal::Geometry;

mod common;

use common::fixture;

#[test]
fn cylinder_voxels_in_r_z_theta_order() {
    let f = Cdgs::open(fixture("cdgs/cdgs")).unwrap();
    assert!(f.failed().is_empty(), "{:?}", f.failed());
    let t = f.read(1).unwrap();
    assert!(matches!(t.geometry, Geometry::Cylindrical { .. }));
    assert_eq!(t.result.shape(), [3, 1, 2, 3, 2]);
    for a in 0..2 {
        for b in 0..2 {
            for c in 0..3 {
                let s = (100 * a + 10 * b + c + 1) as f64;
                assert_eq!(t.result[[0, 0, a, c, b]], s);
                assert_eq!(t.result[[2, 0, a, c, b]], 3.0 * s, "the Total");
            }
        }
    }
}

#[test]
fn cells_add_up() {
    let f = Cdgs::open(fixture("cdgs/cdgs")).unwrap();
    let t = f.read(2).unwrap();
    assert_eq!(t.result[[0, 0, 0, 0, 0]], 4.0);
    assert_eq!(t.result[[2, 0, 0, 0, 0]], 10.0);
    assert_eq!(t.result[[2, 0, 1, 0, 0]], 0.0, "a voxel not listed");
}
