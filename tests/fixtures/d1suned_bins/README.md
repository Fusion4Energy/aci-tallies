# D1SUNED discrete bins

`bins.msht` is synthetic, written by `make_bins.py` in the layout of
D1SUNED 4.1.x meshtal files (and, for tally 104, 3.1.4's): 2 x 1 x 1 voxels
whose `result` encodes its bins, `1000 e + 100 t + voxel + 1` (1-based
indices; 9 for a Total), with errors `0.01 (voxel + 1)`.

| Tally | What it has |
|---|---|
| 44 | `Daughter Nuclides list:` in place of energy bins (`ETAG=nuc`) |
| 54 | nuclides in place of time bins (`TTAG=nuc`), energy bins as well |
| 64 | `Decay times:` in place of time bins (no Total), nuclides of energy |
| 74 | `Parent Nuclides list:` |
| 84 | `Cells list:` |
| 94 | a nuclide list that overflowed its fields (`**********`); `input`'s FMESH94 card has the nuclides |
| 104 | `out=ij` with one energy bin printed as `Total Energy Bin` (D1SUNED 3.1.4) |
| 114 | a negative value run into its coordinate (`112.500-1.57688E-11`) |
