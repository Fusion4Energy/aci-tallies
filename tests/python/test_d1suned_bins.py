"""D1SUNED meshtal sections with discrete bins (`tests/fixtures/d1suned_bins`,
see its README): nuclides, parents, cells and decay times in place of energy
or time bins, an overflowed nuclide list, a lone energy bin printed as its
Total, and a negative value run into its coordinate."""

from pathlib import Path

import numpy as np
import pytest
from aci_tallies import Input, meshtal

BINS = Path(__file__).parents[2] / "tests" / "fixtures" / "d1suned_bins"
MSHT, INPUT = str(BINS / "bins.msht"), str(BINS / "input")
TOTAL = 9


def expected(e: list[int], t: list[int]) -> np.ndarray:
    """`result` as the fixture writes it: 1000 e + 100 t + voxel + 1."""
    out = np.empty((len(e), len(t), 2, 1, 1))
    for i, ei in enumerate(e):
        for j, tj in enumerate(t):
            out[i, j, :, 0, 0] = [1000 * ei + 100 * tj + 1, 1000 * ei + 100 * tj + 2]
    return out


def test_every_section_reads():
    f = meshtal.open(MSHT)
    assert f.ids == [44, 54, 64, 74, 84, 104, 114]
    assert [i for i, _ in f.failed] == [94]
    assert "overflowed" in f.failed[0][1] and "input=" in f.failed[0][1]
    for t in f.read_all():
        assert not t.warnings, (t.id, t.warnings)


def test_nuclides_in_place_of_energy():
    t = meshtal.open(MSHT)[44]
    assert t.energy.kind == "nuclide" and t.energy.bounds is None
    np.testing.assert_array_equal(t.energy.labels, [0.001, 26059.002, 27060.003])
    assert t.energy.total and t.energy.n_bins == 3 and len(t.energy) == 4
    assert t.time.kind == "bounds" and len(t.time) == 1
    np.testing.assert_array_equal(t.result, expected([1, 2, 3, TOTAL], [1]))
    np.testing.assert_array_equal(t.rel_error[0, 0, :, 0, 0], [0.01, 0.02])


def test_nuclides_in_place_of_time():
    t = meshtal.open(MSHT)[54]
    assert t.time.kind == "nuclide" and t.time.total
    assert t.energy.kind == "bounds" and t.energy.total
    np.testing.assert_array_equal(t.energy.bounds, [0.0, 1.0, 2.0])
    np.testing.assert_array_equal(t.result, expected([1, 2, TOTAL], [1, 2, 3, TOTAL]))


def test_decay_times_have_no_total():
    t = meshtal.open(MSHT)[64]
    assert t.time.kind == "decay_time" and not t.time.total
    np.testing.assert_array_equal(t.time.labels, [1.0, 2.0])
    assert t.energy.kind == "nuclide" and t.energy.total
    np.testing.assert_array_equal(t.result, expected([1, 2, TOTAL], [1, 2]))


def test_parents_and_cells():
    f = meshtal.open(MSHT)
    assert f[74].energy.kind == "parent"
    np.testing.assert_array_equal(f[74].result, expected([1, 2, 3, TOTAL], [1]))
    assert f[84].energy.kind == "cell"
    np.testing.assert_array_equal(f[84].energy.labels, [0.001, 10.002])
    np.testing.assert_array_equal(f[84].result, expected([1, 2, TOTAL], [1]))


def test_an_overflowed_list_from_the_input():
    f = meshtal.open(MSHT, input=Input(INPUT))
    assert not f.failed
    t = f[94]
    np.testing.assert_allclose(t.energy.labels, [0.001, 41092900.002, 99973182.003])
    np.testing.assert_array_equal(t.result, expected([1, 2, 3, TOTAL], [1]))


def test_a_lone_energy_bin_printed_as_total():
    t = meshtal.open(MSHT)[104]
    assert t.out_format == "ij" and not t.energy.total
    np.testing.assert_array_equal(t.result[0, 0, :, 0, 0], [1.0, 2.0])
    np.testing.assert_array_equal(t.rel_error[0, 0, :, 0, 0], [0.1, 0.2])


def test_a_negative_value_run_into_its_coordinate():
    t = meshtal.open(MSHT)[114]
    assert t.result[0, 0, 1, 0, 0] == -1.57688e-11
    assert t.rel_error[0, 0, 1, 0, 0] == 2.61718


def test_exports_name_the_labels():
    pl = pytest.importorskip("polars")
    df = meshtal.open(MSHT)[44].to_polars()
    assert "energy_nuclide" in df.columns and "energy_low" not in df.columns
    assert df["energy_nuclide"].is_nan().sum() == 2  # the Total's
    assert isinstance(df, pl.DataFrame)
