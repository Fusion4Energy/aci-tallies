"""Write tests/fixtures/cdgs/cdgs: a synthetic CDGS file in D1SUNED 4.1's
layout (fmesh_mod.F90, cdgsV2_print and cuv_cdgs_print), whose values say
where they belong. Usage: python make_cdgs.py <output directory>."""

import sys


def nums(values, fmt="{:12.5E}", per_line=7):
    out = ""
    for i in range(0, len(values), per_line):
        out += "".join(fmt.format(v) for v in values[i : i + per_line]) + "\n"
    return out


def header(
    mesh_id, source, energy_type, energies, mesh_type, bounds, frame, key="mesh_source"
):
    text = f"mesh_id {mesh_id:6d}\n"
    text += "CDGS produced by D1S-UNED simulation.                   synthetic\n"
    text += "cooling_time 1.000E+05\n"
    text += f"{key} {source:10.4E}\n"
    text += f"energy_type {energy_type}\n"
    text += f"energy_boundaries {len(energies):5d}\n"
    text += "".join(e if isinstance(e, str) else f"{e:6.3f}  " for e in energies) + "\n"
    text += f"mesh_type {mesh_type}\n"
    text += "mesh_boundaries " + "".join(f"  {len(b):5d}" for b in bounds) + "\n"
    text += "".join(f"   {v:12.5E}" for v in frame) + "\n"
    for b in bounds:
        text += nums(b, "{:12.5E} ")
    return text + "source_data\n"


def voxel(index, cells, volume):
    """cells: (cell, fraction, values, errors)."""
    total = sum(sum(c[2]) for c in cells)
    text = f"{index:7d}  {total:12.5E}  {volume:12.5E}  {len(cells):3d}\n"
    for cell, fraction, values, errors in cells:
        text += f"{cell:7d}  {fraction:7.5f}  {sum(values):12.5E}\n"
        text += nums(values, "{:12.5E}", 8)
        text += nums(errors, "{:8.6f} ", 8)
    return text, total


IDENTITY = [0, 0, 0, 1, 0, 0, 0, 1, 0]
meshes = []

# 1: a cylinder (bounds printed r, theta, z; voxels r slowest, z fastest).
# A voxel's source is s = 100 a + 10 b + c + 1 (a, b, c its r, theta, z
# indices), s times (bin + 1) in each of 2 energy bins.
r, theta, z = [0, 1, 2], [0, 0.5, 1], [0, 1, 2, 3]
body, total, index = "", 0.0, 0
for a in range(2):
    for b in range(2):
        for c in range(3):
            index += 1
            s = 100 * a + 10 * b + c + 1
            text, t = voxel(index, [(0, 1.0, [s, 2 * s], [0.1, 0.2])], 1.0)
            body, total = body + text, total + t
cyl = [0, 0, 0, 0, 1, 0, -1, 0, 0]  # origin, AXS, VEC
meshes.append(header(1, total, "bins", [0, 1, 2], "cyl", [r, theta, z], cyl) + body)

# 2: voxels of several cells (cuvcdgs): the voxel's source is the sum of its
# cells'; voxel 2 has none and is not listed.
cells = [(5, 0.3, [1.0, 2.0], [0.1, 0.1]), (6, 0.7, [3.0, 4.0], [0.2, 0.2])]
text, total = voxel(1, cells, 8.0)
meshes.append(
    header(2, total, "bins", [0, 1, 2], "rec", [[0, 2, 4], [0, 2], [0, 2]], IDENTITY)
    + text
)

# 3: gamma lines, and an energy bound too large for its field.
text, total = voxel(1, [(0, 1.0, [1.0, 2.0, 3.0], [0, 0, 0])], 1.0)
meshes.append(
    header(3, total, "lines", [0.5, 1.0, 1.5], "rec", [[0, 1]] * 3, IDENTITY) + text
)
text, total = voxel(1, [(0, 1.0, [5.0], [0])], 1.0)
meshes.append(
    header(4, total, "bins", [0.001, "  ******"], "rec", [[0, 1]] * 3, IDENTITY) + text
)

# 5: placed by a TR card: origin, x axis, y axis; its source labelled as F4Enix
# writes it (total_source).
text, total = voxel(1, [(0, 1.0, [7.0], [0])], 1.0)
frame = [1, 2, 3, 0, 1, 0, -1, 0, 0]
meshes.append(
    header(5, total, "bins", [0, 1], "rec", [[0, 1]] * 3, frame, "total_source") + text
)

out = f"num_meshes {len(meshes):5d}\ntotal_source 1.0000E+00\n"
out += "".join(m + "end_source_data\n" for m in meshes)
with open(sys.argv[1] + "/cdgs", "w") as f:
    f.write(out)
