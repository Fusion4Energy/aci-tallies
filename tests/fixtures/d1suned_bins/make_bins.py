"""Write tests/fixtures/d1suned_bins/{bins.msht,input}: small D1SUNED-style
meshtal sections (2x1x1 voxels) whose values encode their bins:
result = 1000*e + 100*t + voxel + 1 (e, t 1-based; Total = 9), error 0.01*(voxel+1)."""

import sys

out = sys.argv[1]
H = (
    " d1sune version 411   ld=05/19/22  probid =  08/29/23 14:45:17\n"
    " D1SUNED discrete bins\n"
    " Number of histories used for normalizing tallies =          1000.00\n"
)
GEOM = (
    "\n Tally bin boundaries:\n"
    "    X direction:      0.00      1.00      2.00\n"
    "    Y direction:      0.00      1.00\n"
    "    Z direction:      0.00      1.00\n"
)
XS = [0.5, 1.5]


def val(e, t, v):
    return f"{1000 * e + 100 * t + v + 1:.5E} {0.01 * (v + 1):.5E}"


def lab(i, n, total):
    return "Total" if (total and i == n) else str(i + 1)


sections = []


def section(id, flags, lines, colhead, e_items, t_items):
    s = (
        f"\n Mesh Tally Number {id:9d}\ndecay photon   mesh tally.\n"
        + "".join(f" {f}\n" for f in flags)
        + GEOM
    )
    s += (
        "".join(f"    {l}\n" for l in lines)
        + "\n"
        + f"   {colhead}    X         Y         Z     Result     Rel Error\n"
    )
    for ei, e in enumerate(e_items):
        for ti, t in enumerate(t_items):
            for v, x in enumerate(XS):
                lead = " ".join(f"{c:>10}" for c in (e[1], t[1]) if c is not None)
                s += f" {lead} {x:9.3f} {0.5:9.3f} {0.5:9.3f} {val(e[0], t[0], v)}\n"
    sections.append(s)


nuc3 = [(1, "1"), (2, "2"), (3, "3"), (9, "Total")]
none = [(1, None)]
section(
    44,
    ["Energy binning is used as user bin."],
    ["Daughter Nuclides list:     0.001 26059.002 27060.003"],
    "Nuclide bin",
    nuc3,
    none,
)
section(
    54,
    ["Time binning is used as user bin."],
    [
        "Daughter Nuclides list:     0.001 26059.002 27060.003",
        "Energy bin boundaries: 0.00E+00 1.00E+00 2.00E+00",
    ],
    "Energy      Nuclide",
    [(1, "1.000E+00"), (2, "2.000E+00"), (9, "Total")],
    nuc3,
)
section(
    64,
    ["Time binning is used as user bin.", "Energy binning is used as user bin."],
    ["Decay times: 1.00E+00 2.00E+00", "Daughter Nuclides list:     0.001 26059.002"],
    "Nuclide      Time",
    [(1, "1"), (2, "2"), (9, "Total")],
    [(1, "1.000E+00"), (2, "2.000E+00")],
)
section(
    74,
    ["Energy binning is used as user bin."],
    ["Parent Nuclides list:     0.001 26059.002 27060.003"],
    "Nuclide bin",
    nuc3,
    none,
)
section(
    84,
    ["Energy binning is used as user bin."],
    ["Cells list:      0.001     10.002"],
    "Cell bin",
    [(1, "1"), (2, "2"), (9, "Total")],
    none,
)
section(
    94,
    ["Energy binning is used as user bin."],
    ["Daughter Nuclides list:     0.001********************"],
    "Nuclide bin",
    nuc3,
    none,
)
# A lone energy bin printed as D1SUNED 3.1.4 does, as its "Total".
sections.append(
    "\n Mesh Tally Number       104\n This is a neutron mesh tally.\n"
    + GEOM
    + "    Energy bin boundaries: 1.00E-11 5.00E+01\n\nTotal Energy Bin\n\n  Z bin:      0.00   -      1.00\n\n"
    "     Tally Results:  X (across) by Y (down)\n                 0.50        1.50\n         0.50 1.00000E+00 2.00000E+00\n\n"
    "     Relative Errors\n                 0.50        1.50\n         0.50 0.10000     0.20000\n\n"
)
# A negative value run into its coordinate.
sections.append(
    "\n Mesh Tally Number       114\n neutron   mesh tally.\n"
    + GEOM.replace("0.00      1.00      2.00", "0.00    112.00    225.00")
    + "    Energy bin boundaries: 0.00E+00 1.00E+36\n\n   X         Y         Z     Result     Rel Error\n"
    "    56.000     0.500     0.500 1.00000E+00 1.00000E-01\n"
    "   168.500     0.500   112.500-1.57688E-11 2.61718E+00\n"
)
with open(out + "/bins.msht", "w") as f:
    f.write(H + "".join(sections))
with open(out + "/input", "w") as f:
    f.write(
        "D1SUNED discrete bins\n1 0 -1\n2 0 1\n\n1 so 100\n\nmode p\nnps 1000\n"
        "FMESH94:p origin 0 0 0 imesh 2 iints 2 jmesh 1 kmesh 1\n"
        "     etag=nuc emesh 41092900 99973182\n"
    )
