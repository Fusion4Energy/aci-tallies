"""F4Enix's meshtal and CuV test files (`F4ENIX_RESOURCES`, its
`tests/resources`), when present: every section reads, with no warning,
but for the ones only F4Enix reads (D1SUNED source-importance meshes) and an
overflowed nuclide list, which needs the input."""

import os
from pathlib import Path

import aci_tallies
import pytest
from aci_tallies import Input, meshtal

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
    input = None
    if path.name == "meshtal_d1s_overflow":
        input = Input(MESHTAL / "input_d1s_overflow")
    f = aci_tallies.open(path, input=input)
    assert {i for i, _ in f.failed} == NOT_READ.get(path.name, set())
    for id in f.ids:
        t = f[id]
        if isinstance(f, meshtal.MeshtalFile):
            assert not t.warnings, (id, t.warnings)
