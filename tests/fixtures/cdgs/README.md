# CDGS

`cdgs` is synthetic, written by `make_cdgs.py` in D1SUNED 4.1's layout
(`cdgsV2_print` and `cuv_cdgs_print` in its `fmesh_mod.F90`):

| Mesh | What it has |
|---|---|
| 1 | a cylinder, its bounds printed r, theta, z and its voxels numbered r, theta, z (z fastest); each voxel's source encodes its indices |
| 2 | voxels of several cells (`out=cuvcdgs`), one voxel not listed |
| 3 | gamma lines (`energy_type lines`) |
| 4 | an energy bound too large for its field (`******`) |
| 5 | a mesh placed by a TR card (origin, x axis, y axis), its source labelled `total_source` as F4Enix writes it |

`tests/python/test_cdgs.py` also reads D1SUNED's own CDGS files when
`D1SUNED_TESTING` names its `Testing/D1SUNED` directory.
