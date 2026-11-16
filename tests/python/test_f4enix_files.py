"""F4Enix's meshtal and CuV test files (`F4ENIX_RESOURCES`, its
`tests/resources`), when present: every section reads, with no warning,
but for the ones only F4Enix reads (D1SUNED source-importance meshes)."""

import os
import re
from pathlib import Path

import aci_tallies
import numpy as np
import pytest
from aci_tallies import cuv, meshtal

RESOURCES = Path(
    os.environ.get("F4ENIX_RESOURCES", Path.home() / "Coding/F4Enix/tests/resources")
)
MESHTAL = RESOURCES / "meshtal"
# D1SUNED source importance (`GEOMETRY IMESH`, target tallies).
NOT_READ = {"test_srcimp": {114, 5114, 8114}}

pytestmark = pytest.mark.skipif(
    not MESHTAL.is_dir(), reason="F4Enix's test files are not here"
)


def files() -> list[Path]:
    if not MESHTAL.is_dir():
        return []
    return sorted(
        p
        for p in MESHTAL.rglob("*")
        if p.is_file() and aci_tallies.detect(p) in ("meshtal", "cuvmsh", "meshinfo")
    )


@pytest.mark.parametrize("path", files(), ids=lambda p: str(p.relative_to(MESHTAL)))
def test_reads(path: Path) -> None:
    f = aci_tallies.open(path)
    assert {i for i, _ in f.failed} == NOT_READ.get(path.name, set())
    for id in f.ids:
        t = f[id]
        if isinstance(f, meshtal.MeshtalFile):
            assert not t.warnings, (id, t.warnings)


def _vtr_total(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """`Value - Total` and `Error - Total` of an ASCII VTR, VTK's cell order."""
    text = path.read_text()

    def array(name: str) -> np.ndarray:
        m = re.search(rf'Name="{name}"[^>]*>(.*?)</DataArray>', text, re.DOTALL)
        assert m is not None, name
        return np.array(m.group(1).split(), dtype=float)

    return array("Value - Total"), array("Error - Total")


@pytest.mark.parametrize(
    ("expected", "kwargs"),
    [
        ("cuvmsh_44_CuV_NONE.vtr", {}),
        ("cuvmsh_44_CuV_VTOT_NONE.vtr", {"norm": "ctot"}),
        ("cuvmsh_44_CuV_CELF_NONE.vtr", {"norm": "celf"}),
        ("cuvmsh_44_CuV_CELF10.vtr", {"cells": [10], "norm": "celf"}),
    ],
)
def test_cuv_mesh_as_f4enix_makes_it(expected: str, kwargs: dict) -> None:
    """`to_mesh(cells=, norm=)` against F4Enix's own VTK files of the same
    cuvmsh (printed with 6 digits)."""
    t = cuv.open(MESHTAL / "test_CuV" / "cuvmsh").to_mesh(44, **kwargs)
    value, error = _vtr_total(MESHTAL / "test_CuV" / expected)
    # VTK cells run x fastest: (i, j, k) -> (k, j, i).
    np.testing.assert_allclose(
        t.result[-1, 0].transpose(2, 1, 0).ravel(), value, rtol=1e-5
    )
    np.testing.assert_allclose(
        t.rel_error[-1, 0].transpose(2, 1, 0).ravel(), error, rtol=1e-4
    )
