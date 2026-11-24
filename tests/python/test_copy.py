"""Every object is immutable: `copy` and `deepcopy` give back the object itself."""

import copy
import inspect
from pathlib import Path

import aci_tallies
from aci_tallies import _native, cdgs, cuv, mctal, meshtal

FIXTURES = Path(__file__).parents[1] / "fixtures"
MESHTAL = Path(__file__).parents[2] / "examples" / "meshtal"
D1S = FIXTURES / "d1suned"


def _objects() -> list[object]:
    m = meshtal.open(MESHTAL)
    t = m[14]
    c = mctal.open(FIXTURES / "mctal" / "kcode.m")
    levels = FIXTURES / "mctal" / "levels"
    bins = mctal.open(levels.with_suffix(".m"), input=levels.with_suffix(".mcnp"))[
        84
    ].bins
    assert bins is not None
    # no FMESH card of this file's tallies in this input: a warning
    warned = meshtal.open(
        MESHTAL, input=aci_tallies.Input.from_text("t\n1 0 -1\n\n1 so 1\n")
    )
    ct = c.read_all()[0]
    u = cuv.open(D1S / "cuvmsh")
    g = cdgs.open(FIXTURES / "cdgs" / "cdgs")
    return [
        aci_tallies.Input(D1S / "input"),
        m,
        m.info(14),
        t,
        t.geometry,
        t.energy,
        c,
        c.info(ct.id),
        ct,
        c.kcode,
        u,
        u.info(4),
        u[4],
        g,
        g.info(g.ids[0]),
        ct.tfc,
        bins[0],
        next(w for t in warned.tallies for w in t.warnings),
    ]


def test_copy_and_deepcopy_are_the_object_itself():
    objects = _objects()
    classes = {
        c
        for c in vars(_native).values()
        if inspect.isclass(c) and not issubclass(c, BaseException)
    }
    assert {type(obj) for obj in objects} == classes
    for obj in objects:
        assert copy.copy(obj) is obj, type(obj)
        assert copy.deepcopy(obj) is obj, type(obj)


def test_deepcopy_of_a_container_shares_them():
    t = meshtal.open(MESHTAL)[14]
    box = [t, {"k": t}]
    new = copy.deepcopy(box)
    assert new is not box and new[1] is not box[1]
    assert new[0] is t and new[1]["k"] is t
    assert new[0].result.sum() == t.result.sum()
