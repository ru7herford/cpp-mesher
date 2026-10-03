"""Pictures for page 3 ("Where this fits") and the README: the airfoil as a curve, as points, as my triangle
mesh, and as a gmsh mesh with prism (boundary) layers; plus a 3D wing surface mesh. gmsh draws only these
pictures; the mesher itself never uses it.

Needs numpy, matplotlib, pillow and gmsh. Run from the repository root after building the mesher:
    ./build/mesher examples/airfoil.json build/airfoil.vtk
    python scripts/bigpicture.py
"""
import json
import sys
from pathlib import Path

import gmsh
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402
from PIL import Image, ImageChops, ImageDraw  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
from plot import read_vtk  # noqa: E402

PANEL, INK, TRI, LAYER, POINT = "#FDFDFB", "#1A1C19", "#C9D9EC", "#D8CDBB", "#1A1C19"
VIEW = (-0.15, 1.15, -0.37, 0.37)  # x0, x1, y0, y1 around the airfoil (chord 1)
FOIL = np.array(json.loads(Path("examples/airfoil.json").read_text())["holes"][0], float)
OUT = Path("docs/img")


def naca0012(n=400):
    """The true NACA 0012 curve (closed trailing edge), for panel 1."""
    x = (1 - np.cos(np.linspace(0, np.pi, n))) / 2
    t = 0.12 / 0.2 * (0.2969 * np.sqrt(x) - 0.126 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    return np.concatenate([np.c_[x, t][::-1], np.c_[x, -t][1:]])


def panel():
    fig, ax = plt.subplots(figsize=(4.4, 4.4 * (VIEW[3] - VIEW[2]) / (VIEW[1] - VIEW[0])), dpi=200)
    ax.set_xlim(VIEW[:2]); ax.set_ylim(VIEW[2:]); ax.set_aspect("equal"); ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)
    return fig, ax


def save(fig, name, trim=False):
    """Write a webp; with trim, cut the empty margin around the drawing (keeping a small border)."""
    png = OUT / f"{name}.png"
    fig.savefig(png, facecolor=PANEL)
    plt.close(fig)
    image = Image.open(png).convert("RGB")
    if trim:
        ink = ImageChops.difference(image, Image.new("RGB", image.size, PANEL)).convert("L").point(lambda v: v > 8 and 255)
        x0, y0, x1, y1 = ink.getbbox()
        image = image.crop((max(x0 - 20, 0), max(y0 - 20, 0), min(x1 + 20, image.width), min(y1 + 20, image.height)))
    image.save(OUT / f"{name}.webp", quality=88)
    png.unlink()
    print("wrote", OUT / f"{name}.webp")


def airfoil_on_top(ax, outline=FOIL):
    ax.fill(*outline.T, color=PANEL, zorder=3)
    ax.plot(*np.vstack([outline, outline[:1]]).T, color=INK, lw=1.4, zorder=4)


def gmsh_prism_layers():
    """Airfoil in the same box as the example, with 12 thin layers against the wall and triangles outside."""
    box = json.loads(Path("examples/airfoil.json").read_text())["outer"]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    def loop(points, size):
        tags = [gmsh.model.geo.addPoint(x, y, 0, size) for x, y in points]
        lines = [gmsh.model.geo.addLine(tags[i], tags[(i + 1) % len(tags)]) for i in range(len(tags))]
        return gmsh.model.geo.addCurveLoop(lines), lines

    outer, _ = loop(box, 0.25)
    inner, wall = loop(FOIL, 0.03)
    gmsh.model.geo.addPlaneSurface([outer, inner])
    gmsh.model.geo.synchronize()
    layers = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(layers, "CurvesList", wall)
    gmsh.model.mesh.field.setNumber(layers, "Size", 0.004)
    gmsh.model.mesh.field.setNumber(layers, "Ratio", 1.25)
    gmsh.model.mesh.field.setNumber(layers, "Thickness", 0.06)
    gmsh.model.mesh.field.setNumber(layers, "Quads", 1)
    gmsh.model.mesh.field.setAsBoundaryLayer(layers)
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xy = {t: coords[3 * i : 3 * i + 2] for i, t in enumerate(tags)}
    cells = []
    for kind, nodes in zip(*[gmsh.model.mesh.getElements(2)[i] for i in (0, 2)]):
        k = 3 if kind == 2 else 4
        cells += [([xy[n] for n in nodes[i : i + k]], k) for i in range(0, len(nodes), k)]
    gmsh.finalize()
    return cells


def gmsh_wing():
    """Surface triangles of a straight wing: the airfoil pushed 2 chords sideways."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    points = [gmsh.model.occ.addPoint(x, y, 0) for x, y in FOIL]
    lines = [gmsh.model.occ.addLine(points[i], points[(i + 1) % len(points)]) for i in range(len(points))]
    face = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(lines)])
    gmsh.model.occ.extrude([(2, face)], 0, 0, 2.0)
    gmsh.model.occ.synchronize()
    gmsh.option.setNumber("Mesh.MeshSizeMax", 0.1)
    gmsh.model.mesh.generate(2)
    tags, coords, _ = gmsh.model.mesh.getNodes()
    xyz = {t: coords[3 * i : 3 * i + 3] for i, t in enumerate(tags)}
    triangles = []
    for kind, nodes in zip(*[gmsh.model.mesh.getElements(2)[i] for i in (0, 2)]):
        if kind == 2:
            triangles += [[xyz[n] for n in nodes[i : i + 3]] for i in range(0, len(nodes), 3)]
    gmsh.finalize()
    return np.array(triangles)


def main():
    OUT.mkdir(exist_ok=True)
    # 1. The shape: the true curve.
    fig, ax = panel()
    airfoil_on_top(ax, naca0012())
    save(fig, "where-1-curve")
    # 2. The curve cut into points.
    fig, ax = panel()
    airfoil_on_top(ax)
    ax.plot(*FOIL.T, "o", color=POINT, ms=3.2, zorder=5)
    save(fig, "where-2-points")
    # 3. My triangles.
    points, triangles, _ = read_vtk("build/airfoil.vtk")
    fig, ax = panel()
    ax.add_collection(PolyCollection(points[triangles], facecolors=TRI, edgecolors=INK, linewidths=0.35))
    airfoil_on_top(ax)
    save(fig, "where-3-triangles")
    # 4. Prism layers (gmsh).
    cells = gmsh_prism_layers()
    fig, ax = panel()
    ax.add_collection(PolyCollection([c for c, _ in cells], facecolors=[LAYER if k == 4 else TRI for _, k in cells],
                                     edgecolors=INK, linewidths=0.3))
    airfoil_on_top(ax)
    save(fig, "where-4-layers")
    # 3D: wing surface mesh, shaded by how much each triangle faces the light.
    triangles = gmsh_wing()[:, :, [0, 2, 1]]  # x along the chord, y along the span, z up
    normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    normals /= np.linalg.norm(normals, axis=1, keepdims=True)
    light = np.array([-0.4, -0.6, 0.7]) / np.linalg.norm([-0.4, -0.6, 0.7])
    shade = 0.55 + 0.45 * np.abs(normals @ light)
    base = np.array([int(TRI[i : i + 2], 16) / 255 for i in (1, 3, 5)])
    fig = plt.figure(figsize=(6, 4), dpi=200)
    ax = fig.add_subplot(projection="3d")
    ax.set_facecolor(PANEL)
    ax.add_collection3d(Poly3DCollection(triangles, facecolors=np.clip(base * shade[:, None], 0, 1),
                                         edgecolors=INK, linewidths=0.15))
    ax.set_xlim(0, 1); ax.set_ylim(0, 2); ax.set_zlim(-0.3, 0.3); ax.set_box_aspect((1, 2, 0.6))
    ax.view_init(elev=22, azim=-128)
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)
    save(fig, "where-3d-wing", trim=True)
    # The README shows the four 2D pictures side by side, the third (this project) outlined in the page accent.
    panels = [Image.open(OUT / f"where-{n}.webp") for n in ("1-curve", "2-points", "3-triangles", "4-layers")]
    w, h, gap = *panels[0].size, 24
    strip = Image.new("RGB", (4 * w + 3 * gap, h), "white")
    for k, panel_image in enumerate(panels):
        strip.paste(panel_image, (k * (w + gap), 0))
    ImageDraw.Draw(strip).rectangle((2 * (w + gap), 0, 2 * (w + gap) + w - 1, h - 1), outline="#A3213B", width=6)
    strip.resize((1800, round(1800 * h / strip.width)), Image.LANCZOS).save("media/where.png", optimize=True)
    print("wrote media/where.png")


if __name__ == "__main__":
    main()
