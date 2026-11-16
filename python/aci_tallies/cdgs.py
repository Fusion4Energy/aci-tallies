"""D1SUNED Common Decay Gamma Source (CDGS) files: decay gamma sources on
meshes, as D1SUNED (and F4Enix) write them.

>>> from aci_tallies import cdgs
>>> f = cdgs.open("cdgs")
>>> f.ids                          # mesh ids
>>> f.info(1).cooling_time         # s; also strength (photons/s), comment
>>> src = f[1]                     # MeshTally: source per voxel and energy bin
>>> src.result.shape               # (energy, 1, i, j, k), Total energy bin last

Each voxel's source is the sum of its cells' (photons/s); voxels the file
leaves out have none. `energy` has the bins (`kind == "bounds"`) or the
gamma lines (`kind == "line"`, MeV) of `energy_type`. Cylindrical meshes
are given in this package's (r, z, theta) order.
"""

from __future__ import annotations

import os

from ._native import CdgsFile, CdgsInfo, open_cdgs

__all__ = ["CdgsFile", "CdgsInfo", "open"]


def open(path: str | os.PathLike[str]) -> CdgsFile:
    """Open a CDGS file."""
    return open_cdgs(path)
