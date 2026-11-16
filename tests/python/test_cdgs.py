"""D1SUNED CDGS files (`tests/fixtures/cdgs`, written by its make_cdgs.py in
D1SUNED 4.1's layout)."""

import os
from pathlib import Path

import aci_tallies
import numpy as np
import pytest
from aci_tallies import cdgs

CDGS = Path(__file__).parents[1] / "fixtures" / "cdgs" / "cdgs"


@pytest.fixture(scope="module")
def f() -> cdgs.CdgsFile:
    f = aci_tallies.open(CDGS)
    assert isinstance(f, cdgs.CdgsFile)
    return f


def test_file(f: cdgs.CdgsFile) -> None:
    assert aci_tallies.detect(CDGS) == "cdgs"
    assert f.ids == [1, 2, 3, 4, 5] and f.failed == []
    assert f.total_source == 1.0
    info = f.info(1)
    assert info.cooling_time == 1e5 and info.strength == pytest.approx(2052.0)
    assert info.comment.startswith("CDGS produced by D1S-UNED")


def test_cylinder_in_r_z_theta_order(f: cdgs.CdgsFile) -> None:
    t = f[1]
    g = t.geometry
    assert g.kind == "cylindrical" and tuple(g.shape) == (2, 3, 2)
    np.testing.assert_array_equal(g.bounds[1], [0, 1, 2, 3])  # z
    np.testing.assert_array_equal(g.bounds[2], [0, 0.5, 1])  # theta
    assert g.axis == [0, 1, 0] and g.vec == [-1, 0, 0]
    assert t.out_format == "cdgs" and t.energy.total
    for a in range(2):
        for b in range(2):
            for c in range(3):
                s = 100 * a + 10 * b + c + 1
                # result: (energy, time, r, z, theta)
                np.testing.assert_allclose(t.result[:, 0, a, c, b], [s, 2 * s, 3 * s])
                np.testing.assert_allclose(t.rel_error[:2, 0, a, c, b], [0.1, 0.2])
    # the Total's error: that of the sum of the bins
    np.testing.assert_allclose(t.rel_error[2, 0, 0, 0, 0], np.hypot(0.1, 0.4) / 3)


def test_voxels_of_several_cells(f: cdgs.CdgsFile) -> None:
    t = f[2]
    # the voxel's source is the sum of its cells'; voxel 2 is not listed
    np.testing.assert_allclose(t.result[:, 0, :, 0, 0], [[4, 0], [6, 0], [10, 0]])
    np.testing.assert_allclose(t.rel_error[0, 0, 0, 0, 0], np.hypot(0.1, 0.6) / 4)
    assert t.volume is not None and t.volume[0, 0, 0] == 8.0


def test_lines_and_overflowed_bounds(f: cdgs.CdgsFile) -> None:
    lines = f[3]
    assert lines.energy.kind == "line" and lines.energy.total
    np.testing.assert_array_equal(lines.energy.labels, [0.5, 1.0, 1.5])
    np.testing.assert_allclose(lines.result[:, 0, 0, 0, 0], [1, 2, 3, 6])
    over = f[4]
    bounds = over.energy.bounds
    assert bounds is not None and np.isnan(bounds[1]) and not over.energy.total
    assert over.result[0, 0, 0, 0, 0] == 5.0
    assert any("overflowed" in w.message for w in over.warnings)


def test_frame_of_a_tr_card(f: cdgs.CdgsFile) -> None:
    t = f[5]
    m = np.asarray(t.geometry.transform)
    # local x along global y, local y along -x, origin (1, 2, 3)
    np.testing.assert_allclose(m @ [1, 0, 0, 1], [1, 3, 3, 1], atol=1e-12)
    np.testing.assert_allclose(m @ [0, 1, 0, 1], [0, 2, 3, 1], atol=1e-12)
    assert any("TR card" in w.message for w in t.warnings)
    assert f.info(5).strength == 7.0  # labelled total_source, as F4Enix writes it


D1SUNED = os.environ.get("D1SUNED_TESTING")


@pytest.mark.skipif(not D1SUNED, reason="D1SUNED_TESTING (its Testing/D1SUNED) not set")
def test_d1suned_files() -> None:
    """D1SUNED's own CDGS files: every voxel's Total adds up to the mesh's
    source."""
    assert D1SUNED is not None
    files = [p for p in Path(D1SUNED).iterdir() if p.is_file()]
    read = 0
    for path in files:
        if aci_tallies.detect(path) != "cdgs":
            continue
        f = cdgs.open(path)
        assert f.failed == []
        for id in f.ids:
            t = f[id]
            np.testing.assert_allclose(
                t.result[-1].sum(), f.info(id).strength, rtol=1e-4, atol=1e-30
            )
            read += 1
    assert read > 0
