"""Plot a mesher result: the mesh coloured by triangle minimum angle, and an angle histogram.

With --stages it instead reads a steps.json file and draws the five steps (the points on the shape's
edges, points connected, edges put back, outside removed, thin triangles fixed), one small PNG each.

usage: python scripts/plot.py mesh.vtk [report.json] [-o mesh.png]
       python scripts/plot.py --stages steps.json [-o folder] [--dpi 100]
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.tri import Triangulation  # noqa: E402

# The web page's palette (docs/deck.css).
SURFACE = "#FDFDFB"
INK = "#1A1C19"
INK_SECONDARY = "#40443D"
MUTED = "#5B6057"
GRID = "#D3D8CD"
AXIS = "#82897E"
BAR = "#1D5FA6"
# One hue, light (small angle) to dark (large angle).
BLUES = LinearSegmentedColormap.from_list("blues", ["#E3ECF6", "#9DB9DA", "#5786BE", "#1D5FA6", "#0F355E"])
# Thin triangles (an angle under the target) and highlighted edges in the stage pictures.
THIN = "#EBC9A5"
THIN_EDGE = "#A85600"
# Bin width of the all-angles histogram when no report is given; the report states its own.
ANGLE_BIN_DEG = 10


def read_vtk(path):
    """Parse the legacy ASCII files written by `mesher` (points, triangles, min_angle)."""
    tokens = Path(path).read_text().split()
    points = triangles = min_angle = None
    i = 0
    while i < len(tokens):
        if tokens[i] == "POINTS":
            n = int(tokens[i + 1])
            points = np.array(tokens[i + 3 : i + 3 + 3 * n], float).reshape(n, 3)[:, :2]
            i += 3 + 3 * n
        elif tokens[i] == "CELLS":
            n = int(tokens[i + 1])
            triangles = np.array(tokens[i + 3 : i + 3 + 4 * n], int).reshape(n, 4)[:, 1:]
            i += 3 + 4 * n
        elif tokens[i] == "LOOKUP_TABLE" and triangles is not None:
            n = len(triangles)
            min_angle = np.array(tokens[i + 2 : i + 2 + n], float)
            i += 2 + n
        else:
            i += 1
    return points, triangles, min_angle


def all_angles(points, triangles):
    p = points[triangles]
    angles = []
    for k in range(3):
        u = p[:, (k + 1) % 3] - p[:, k]
        v = p[:, (k + 2) % 3] - p[:, k]
        cross = np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
        angles.append(np.degrees(np.arctan2(cross, (u * v).sum(axis=1))))
    return np.concatenate(angles)


def replay_stages(run):
    """Replay steps.json and return the triangles after each of the five steps, plus what to highlight."""
    n = len(run["points"])
    build = run["build"]
    live = {(n, n + 1, n + 2)}  # the big enclosing triangle the mesh starts from

    def apply(change):
        live.difference_update(tuple(t) for t in change["removed"])
        live.update(tuple(t) for t in change["added"])

    for change in build["insertions"]:
        apply(change)
    connected = set(live)
    for change in build["recoveries"]:
        apply(change)
    forced = set(live)
    apply(build["exterior"])
    apply(build["legalize"])
    trimmed = set(live)
    for step in run["steps"]:
        apply(step)
    build_points = run["points"] + build["enclosing"]
    final_points = run["points"] + [s["point"] for s in run["steps"]]
    recovered = [r["segment"] for r in build["recoveries"]]
    return [
        ("0-points", build_points, set(), []),
        ("1-connect", build_points, connected, []),
        ("2-force-edges", build_points, forced, recovered),
        ("3-remove-outside", build_points, trimmed, []),
        ("4-refine", final_points, live, []),
    ]


def draw_stages(steps_path, folder, dpi=100):
    """One small picture per step: thin triangles (under the target) in orange, the rest by smallest angle."""
    run = json.loads(Path(steps_path).read_text())
    target = run["input"]["min_angle_deg"]
    loops = [run["input"]["outer"], *run["input"]["holes"]]
    corners = np.array([p for loop in loops for p in loop], float)
    low, high = corners.min(axis=0), corners.max(axis=0)
    margin = 0.08 * (high - low).max()
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for name, points, tris, highlight in replay_stages(run):
        points = np.array(points, float)
        tris = np.array(sorted(tris), int).reshape(-1, 3)
        fig, ax = plt.subplots(figsize=(3.2, 3.2), facecolor=SURFACE)
        # Before the outside is removed the mesh still reaches the big enclosing triangle, so angles mean
        # nothing yet: draw edges only. After that, fill each triangle by its smallest angle.
        filled = len(tris) and (tris.max() < len(points) - 3 or name == "4-refine")
        smallest = all_angles(points, tris).reshape(3, -1).min(axis=0) if len(tris) else []
        for t, a in zip(tris, smallest):
            if not filled:
                colour = "none"
            elif a < target:
                colour = THIN
            else:
                colour = BLUES(0.1 + 0.5 * (a - target) / (60 - target))
            ax.fill(*points[t].T, facecolor=colour, edgecolor=INK_SECONDARY, linewidth=0.6)
        for loop in loops:
            ax.add_patch(plt.Polygon(loop, closed=True, fill=False, edgecolor=INK, linewidth=1.6))
        for a, b in highlight:
            ax.plot(*points[[a, b]].T, color=THIN_EDGE, linewidth=3.5, solid_capstyle="round")
        ax.plot(*points[: len(run["points"])].T, "o", color=INK, markersize=2.5)
        ax.set_xlim(low[0] - margin, high[0] + margin)
        ax.set_ylim(low[1] - margin, high[1] + margin)
        ax.set_aspect("equal")
        ax.axis("off")
        output = folder / f"step{name}.png"
        fig.savefig(output, dpi=dpi, facecolor=SURFACE, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        print(f"wrote {output}")


def style(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=9)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("vtk", nargs="?")
    parser.add_argument("report", nargs="?")
    parser.add_argument("-o", "--output", help="PNG path (default: next to the .vtk); with --stages, a folder")
    parser.add_argument("--stages", metavar="STEPS_JSON", help="draw the five step pictures from a steps.json file")
    parser.add_argument("--dpi", type=int, default=100, help="resolution of the step pictures (default 100)")
    args = parser.parse_args()
    if args.stages:
        draw_stages(args.stages, args.output or ".", args.dpi)
        return
    if not args.vtk:
        parser.error("give a mesh.vtk file, or --stages steps.json")

    points, triangles, min_angle = read_vtk(args.vtk)
    report = json.loads(Path(args.report).read_text()) if args.report else None
    target = report["targets"]["min_angle_deg"] if report else None

    fig, (mesh_ax, hist_ax) = plt.subplots(
        1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.35, 1]}, facecolor=SURFACE
    )

    tri = Triangulation(points[:, 0], points[:, 1], triangles)
    shading = mesh_ax.tripcolor(tri, facecolors=min_angle, cmap=BLUES, vmin=0, vmax=60)
    mesh_ax.triplot(tri, color=SURFACE, linewidth=0.5)
    mesh_ax.set_aspect("equal")
    mesh_ax.axis("off")
    title = f"{len(points)} points, {len(triangles)} triangles"
    mesh_ax.set_title(title, color=INK, fontsize=12, loc="left")
    bar = fig.colorbar(shading, ax=mesh_ax, fraction=0.04, pad=0.02)
    bar.set_label("triangle minimum angle (deg)", color=INK_SECONDARY, fontsize=9)
    bar.ax.tick_params(colors=MUTED, labelcolor=INK_SECONDARY, labelsize=8)
    bar.outline.set_visible(False)
    if target is not None:
        bar.ax.axhline(target, color=INK, linewidth=1.5)

    style(hist_ax)
    if report:
        width = report["angle_histogram"]["bin_width_deg"]
        counts = np.array(report["angle_histogram"]["counts"])
    else:
        width = ANGLE_BIN_DEG
        counts, _ = np.histogram(all_angles(points, triangles), bins=np.arange(0, 180 + width, width))
    edges = np.arange(len(counts)) * width
    hist_ax.bar(edges, counts, width=width, align="edge", color=BAR, edgecolor=SURFACE, linewidth=2)
    hist_ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    hist_ax.set_axisbelow(True)
    hist_ax.set_xlim(0, 180)
    hist_ax.set_xticks(range(0, 181, 30))
    hist_ax.set_xlabel("triangle angle (deg)", color=INK_SECONDARY, fontsize=9)
    hist_ax.set_ylabel("number of angles", color=INK_SECONDARY, fontsize=9)
    hist_ax.set_title(f"All triangle angles, {width:g}-degree bins", color=INK, fontsize=12, loc="left")
    if report:
        hist_ax.axvline(target, color=INK, linewidth=1.5, linestyle="--")
        hist_ax.text(
            0.98,
            0.95,
            f"dashed line: target minimum {target:g}°\nsmallest angle in mesh: {report['min_angle_deg']:.2f}°",
            transform=hist_ax.transAxes,
            ha="right",
            va="top",
            color=INK_SECONDARY,
            fontsize=9,
        )

    output = args.output or str(Path(args.vtk).with_suffix(".png"))
    fig.tight_layout()
    fig.savefig(output, dpi=150, facecolor=SURFACE)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
