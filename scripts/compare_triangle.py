"""Mesh every example with this mesher and with Shewchuk's Triangle, and compare counts and angles.

Triangle is not part of this repository (its licence does not allow redistribution). Point this script at
your own build with --triangle PATH, the TRIANGLE environment variable, or a `triangle` on PATH.

Both meshers get the same boundary: the points this mesher triangulates (each side already divided into pieces
no longer than max_edge, if there is one) joined into the same loops. Triangle runs with -p (keep the boundary),
-q<angle> (the same minimum angle) and, when there is a size limit h, -a<area> with the area of an equilateral
triangle of side h, sqrt(3)/4 h^2, because Triangle limits area, not edge length.

usage: python scripts/compare_triangle.py [--triangle PATH] [--out docs/data/triangle.json]
"""

import argparse
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ["airfoil", "slot", "sharp_corner", "square_hole", "slanted_hole"]


def find_triangle(given):
    path = given or os.environ.get("TRIANGLE") or shutil.which("triangle")
    if not path or not Path(path).is_file():
        raise SystemExit("Triangle not found: pass --triangle PATH, set TRIANGLE, or put `triangle` on PATH")
    return str(Path(path).resolve())


def oriented(loop, outer):
    area = sum(loop[i][0] * loop[(i + 1) % len(loop)][1] - loop[(i + 1) % len(loop)][0] * loop[i][1] for i in range(len(loop)))
    return loop if (area > 0) == outer else loop[::-1]


def ring_sizes(spec):
    """Number of boundary points per loop, as src/refine.cpp:build_boundary samples them."""
    h = spec.get("max_edge") or 0
    sizes = []
    for loop in [spec["outer"], *spec.get("holes", [])]:
        n = 0
        for k, a in enumerate(loop):
            b = loop[(k + 1) % len(loop)]
            n += max(1, math.ceil(math.dist(a, b) / h)) if h > 0 else 1
        sizes.append(n)
    return sizes


def inside(loop, p):
    result = False
    j = len(loop) - 1
    for i in range(len(loop)):
        a, b = loop[i], loop[j]
        if (a[1] > p[1]) != (b[1] > p[1]) and p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]:
            result = not result
        j = i
    return result


def point_in_hole(hole):
    """Any point strictly inside the hole: the vertex average, or else a point near a convex corner."""
    c = [sum(p[0] for p in hole) / len(hole), sum(p[1] for p in hole) / len(hole)]
    if inside(hole, c):
        return c
    for k, p in enumerate(hole):
        q, r = hole[k - 1], hole[(k + 1) % len(hole)]
        m = [(q[0] + p[0] + r[0]) / 3, (q[1] + p[1] + r[1]) / 3]
        if inside(hole, m):
            return m
    raise ValueError("could not find a point inside a hole")


def angles(p, q, r):
    """The three angles in degrees, with the formula of src/quality.cpp:triangle_angles."""
    out = []
    for c, n, m in ((p, q, r), (q, r, p), (r, p, q)):
        ux, uy, vx, vy = n[0] - c[0], n[1] - c[1], m[0] - c[0], m[1] - c[1]
        out.append(math.degrees(math.atan2(abs(ux * vy - uy * vx), ux * vx + uy * vy)))
    return out


def stats(points, triangles, target):
    smallest = [min(angles(*(points[i] for i in t))) for t in triangles]
    largest = [max(angles(*(points[i] for i in t))) for t in triangles]
    edges = {tuple(sorted((t[k], t[(k + 1) % 3]))) for t in triangles for k in range(3)}
    return {
        "points": len(points),
        "triangles": len(triangles),
        "min_angle_deg": round(min(smallest), 4),
        "max_angle_deg": round(max(largest), 4),
        "below_target": sum(a < target for a in smallest),
        "longest_edge": round(max(math.dist(points[a], points[b]) for a, b in edges), 6),
    }


def read_table(path, width):
    rows = [line.split() for line in Path(path).read_text().splitlines() if line.strip() and not line.startswith("#")]
    count = int(rows[0][0])
    return [[float(x) for x in row[1 : 1 + width]] for row in rows[1 : 1 + count]]


def run_triangle(triangle, work, boundary, sizes, spec):
    lines = [f"{len(boundary)} 2 0 0"]
    lines += [f"{i + 1} {p[0]!r} {p[1]!r}" for i, p in enumerate(boundary)]
    segments, start = [], 0
    for n in sizes:
        segments += [(start + k, start + (k + 1) % n) for k in range(n)]
        start += n
    lines.append(f"{len(segments)} 0")
    lines += [f"{i + 1} {a + 1} {b + 1}" for i, (a, b) in enumerate(segments)]
    holes = [point_in_hole(oriented(h, False)) for h in spec.get("holes", [])]
    lines.append(str(len(holes)))
    lines += [f"{i + 1} {p[0]!r} {p[1]!r}" for i, p in enumerate(holes)]
    poly = work / "input.poly"
    poly.write_text("\n".join(lines) + "\n")
    switches = f"-pq{spec['min_angle_deg']:g}"
    if spec.get("max_edge"):
        switches += f"a{math.sqrt(3) / 4 * spec['max_edge'] ** 2:.10g}"
    switches += "Q"
    subprocess.run([triangle, switches, str(poly)], check=True)
    nodes = read_table(work / "input.1.node", 2)
    elements = read_table(work / "input.1.ele", 3)
    first = int(Path(work / "input.1.node").read_text().split("\n")[1].split()[0])  # 1 unless built with -z
    triangles = [[int(v) - first for v in t] for t in elements]
    for f in work.glob("input.1.*"):  # Triangle's own output files carry its full command line; keep only the numbers
        f.unlink()
    return switches, nodes, triangles


def read_vtk(path):
    tokens = Path(path).read_text().split()
    i = tokens.index("POINTS")
    n = int(tokens[i + 1])
    points = [[float(tokens[i + 3 + 3 * k]), float(tokens[i + 4 + 3 * k])] for k in range(n)]
    j = tokens.index("CELLS")
    m = int(tokens[j + 1])
    triangles = [[int(tokens[j + 4 + 4 * k + c]) for c in range(3)] for k in range(m)]
    return points, triangles


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--triangle", help="path to the triangle executable")
    parser.add_argument("--mesher", default=str(ROOT / "build" / "mesher"))
    parser.add_argument("--out", default=str(ROOT / "docs" / "data" / "triangle.json"))
    args = parser.parse_args()
    triangle = find_triangle(args.triangle)
    version = "1.6"  # Triangle prints no version; 1.6 is the current release

    rows = []
    for name in EXAMPLES:
        work = ROOT / "build" / "triangle" / name
        work.mkdir(parents=True, exist_ok=True)
        spec = json.loads((ROOT / "examples" / f"{name}.json").read_text())
        subprocess.run([args.mesher, str(ROOT / "examples" / f"{name}.json"), str(work / "mesh.vtk"),
                        "--steps", str(work / "steps.json")], check=True, capture_output=True)
        boundary = json.loads((work / "steps.json").read_text())["points"]
        sizes = ring_sizes(spec)
        assert sum(sizes) == len(boundary), f"{name}: boundary sampling does not match the mesher's"
        mine = stats(*read_vtk(work / "mesh.vtk"), spec["min_angle_deg"])
        switches, nodes, triangles = run_triangle(triangle, work, boundary, sizes, spec)
        theirs = stats(nodes, triangles, spec["min_angle_deg"])
        row = {"example": name, "target_deg": spec["min_angle_deg"], "max_edge": spec.get("max_edge"),
               "boundary_points": len(boundary), "triangle_switches": switches, "mine": mine, "triangle": theirs,
               "triangle_ratio": round(mine["triangles"] / theirs["triangles"], 3)}
        rows.append(row)
        print(f"{name:13s} mine {mine['points']:5d} pts {mine['triangles']:5d} tris min {mine['min_angle_deg']:7.3f}  "
              f"Triangle {theirs['points']:5d} pts {theirs['triangles']:5d} tris min {theirs['min_angle_deg']:7.3f}  "
              f"ratio {row['triangle_ratio']}")

    Path(args.out).write_text(json.dumps({
        "triangle_version": version,
        "note": "Same boundary points for both; Triangle -p -q<target>, plus -a<sqrt(3)/4 h^2> when max_edge = h.",
        "examples": rows,
    }, indent=1) + "\n")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
