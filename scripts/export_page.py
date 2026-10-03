"""Write the web page's data (docs/data/) from the mesher's own output.

Not part of the five steps: it runs the C++ program and the tests, and keeps what they print.

For every example it runs build/mesher and keeps its report.json and steps.json, times the run (median of
--runs runs, recording off), and replays steps.json to get the mesh quality after every refinement step
(progress.json). It also runs the three test programs and keeps their summaries (checks.json), and runs every
example at higher angle targets to find where refinement stops converging (limits.json).

usage: python scripts/export_page.py [--build build] [--runs 15]
"""

import argparse
import datetime
import json
import math
import re
import statistics
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "docs" / "data"

# Shown on the page in this order. `focus` is the zoomed view (x0, y0, x1, y1).
EXAMPLES = [
    {"id": "airfoil", "title": "Airfoil", "short": "airfoil", "focus": [-0.15, -0.35, 1.15, 0.35], "focus_label": "Zoom on airfoil",
     "note": "A NACA 0012 airfoil (48 points, cosine spacing) in a channel, with no size limit: the triangles grade "
             "themselves outwards from the airfoil."},
    {"id": "slot", "title": "Slot", "short": "slot",
     "note": "A narrow slot with a slanted end. It is the one example where every construction stage has work to do."},
    {"id": "sharp_corner", "title": "Sharp corner", "short": "sharp corner", "focus": [-0.1, -0.12, 1.5, 0.45], "focus_label": "Zoom on corner",
     "note": "A polygon with a 10° corner, sharper than the 30° target, and a size limit of 0.3."},
    {"id": "square_hole", "title": "Square with hole", "short": "square",
     "note": "A 4 × 4 square with a 1 × 1 hole and a size limit of 0.5."},
    {"id": "slanted_hole", "title": "Slanted hole", "short": "slanted plate",
     "note": "A slanted four-sided plate with a slanted four-sided hole and a size limit of 0.3: no side is "
             "parallel to an axis."},
]
SMALLEST_ANGLE_BIN_DEG = 5  # bins of the page's "smallest angle per triangle" histogram, 0 to 60 degrees
SWEEP_TARGETS = [31, 32, 33, 34, 35, 36]


def angles(p, q, r):
    """The three angles in degrees, computed exactly as src/quality.cpp:triangle_angles does."""
    out = []
    for c, n, m in ((p, q, r), (q, r, p), (r, p, q)):
        ux, uy, vx, vy = n[0] - c[0], n[1] - c[1], m[0] - c[0], m[1] - c[1]
        out.append(math.atan2(abs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / math.pi)
    return out


def progress(run):
    """Quality of the mesh at refinement step 0 (the first full mesh) and after every step."""
    points = run["points"] + [s["point"] for s in run["steps"]]
    target = run["input"]["min_angle_deg"]
    bins = 60 // SMALLEST_ANGLE_BIN_DEG
    live = {tuple(t): min(angles(*(points[i] for i in t))) for t in run["triangles"]}
    series = {"below": [], "min_angle": [], "triangles": [], "hist": []}

    def record():
        smallest = list(live.values())
        hist = [0] * bins
        for a in smallest:
            hist[min(bins - 1, int(a // SMALLEST_ANGLE_BIN_DEG))] += 1
        series["below"].append(sum(a < target for a in smallest))
        series["min_angle"].append(round(min(smallest), 3))
        series["triangles"].append(len(smallest))
        series["hist"].append(hist)

    record()
    for step in run["steps"]:
        for t in step["removed"]:
            del live[tuple(t)]
        for t in step["added"]:
            live[tuple(t)] = min(angles(*(points[i] for i in t)))
        record()
    return {"target_deg": target, "bin_deg": SMALLEST_ANGLE_BIN_DEG, **series}


def mesher(build, spec_path, out_dir, steps=True):
    out_dir.mkdir(parents=True, exist_ok=True)
    command = [str(build / "mesher"), str(spec_path), str(out_dir / "mesh.vtk"), "--report", str(out_dir / "report.json")]
    if steps:
        command += ["--steps", str(out_dir / "steps.json")]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 2:
        raise SystemExit(f"mesher failed on {spec_path}: {result.stderr.strip()}")
    return result.returncode, json.loads((out_dir / "report.json").read_text())


def summary_line(output, program):
    match = re.search(r"^SUMMARY (\{.*\})$", output, re.MULTILINE)
    if not match:
        raise SystemExit(f"{program} printed no SUMMARY line")
    return json.loads(match.group(1))


def run_tests(build, scratch):
    results = {}
    for name, args in (("examples", [str(ROOT / "examples"), str(scratch)]), ("predicates", []), ("fuzz", [])):
        program = build / f"test_{name}"
        result = subprocess.run([str(program), *args], capture_output=True, text=True)
        results[name] = summary_line(result.stdout, program.name)
        if result.returncode != 0 or results[name]["failures"]:
            raise SystemExit(f"{program.name} failed:\n{result.stdout}")
        print(f"  test_{name}: passed")
    return results


def toolchain(build):
    """Compiler and CMake versions, as CMake recorded them for this build."""
    info = next(build.glob("CMakeFiles/*/CMakeCXXCompiler.cmake")).read_text()
    field = lambda name: re.search(rf'set\(CMAKE_CXX_COMPILER_{name} "([^"]*)"\)', info).group(1)
    names = {"GNU": "GCC", "Clang": "Clang", "AppleClang": "Apple Clang", "MSVC": "MSVC"}
    cache = (build / "CMakeCache.txt").read_text()
    part = lambda name: re.search(rf"^CMAKE_CACHE_{name}_VERSION:INTERNAL=(\d+)$", cache, re.MULTILINE).group(1)
    version = ".".join(part(p) for p in ("MAJOR", "MINOR", "PATCH"))
    return {"compiler": f"{names.get(field('ID'), field('ID'))} {field('VERSION')}", "cmake": version}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--build", default=str(ROOT / "build"), help="CMake build folder with mesher and the tests")
    parser.add_argument("--runs", type=int, default=15, help="timed runs per example (the median is kept)")
    args = parser.parse_args()
    build = Path(args.build).resolve()

    with tempfile.TemporaryDirectory() as tmp:
        scratch = Path(tmp)
        print("tests")
        checks = run_tests(build, scratch)

        listing = []
        for example in EXAMPLES:
            name = example["id"]
            spec = ROOT / "examples" / f"{name}.json"
            out = DATA / name
            code, report = mesher(build, spec, scratch / name)
            if code != 0:
                raise SystemExit(f"{name}: targets not met")
            out.mkdir(parents=True, exist_ok=True)
            for f in ("report.json", "steps.json"):
                (out / f).write_text((scratch / name / f).read_text())
            run = json.loads((out / "steps.json").read_text())
            (out / "progress.json").write_text(json.dumps(progress(run), separators=(",", ":")) + "\n")

            # Timing: recording off, so this is the cost of meshing alone (triangulation + refinement).
            times = []
            for _ in range(args.runs):
                _, timed = mesher(build, spec, scratch / f"{name}_timed", steps=False)
                times.append(timed["time_ms"]["triangulate"] + timed["time_ms"]["refine"])
            listing.append({**example, "time_ms": round(statistics.median(times), 2), "timed_runs": args.runs,
                            "points": report["points"], "triangles": report["triangles"],
                            "min_angle_deg": report["min_angle_deg"], "max_angle_deg": report["max_angle_deg"],
                            "steps": len(run["steps"]), "boundary_points": len(run["points"])})
            print(f"  {name}: {report['points']} points, {report['triangles']} triangles, "
                  f"min angle {report['min_angle_deg']:.3f}, {len(run['steps'])} steps, median {listing[-1]['time_ms']} ms")

        print("angle sweep")
        sweep = []
        for example in EXAMPLES:
            spec = json.loads((ROOT / "examples" / f"{example['id']}.json").read_text())
            for target in SWEEP_TARGETS:
                path = scratch / f"sweep_{example['id']}_{target}.json"
                path.write_text(json.dumps({**spec, "min_angle_deg": target}))
                code, report = mesher(build, path, scratch / "sweep", steps=False)
                r = report["refinement"]
                sweep.append({"example": example["id"], "target_deg": target, "met": code == 0,
                              "points": report["points"], "min_angle_deg": round(report["min_angle_deg"], 3),
                              "hit_point_limit": r["hit_point_limit"],
                              "below_target": r["triangles_below_min_angle"] - r["of_which_at_sharp_input_corners"]})
                print(f"  {example['id']} at {target}: {'met' if code == 0 else 'not met'}, {report['points']} points")

    recorded = datetime.date.today().isoformat()
    (DATA / "examples.json").write_text(json.dumps({"recorded": recorded, **toolchain(build), "examples": listing}, indent=1, ensure_ascii=False) + "\n")
    (DATA / "checks.json").write_text(json.dumps(checks, indent=1) + "\n")
    (DATA / "limits.json").write_text(json.dumps({"max_points": 100000, "runs": sweep}, indent=1) + "\n")
    print(f"wrote {DATA}")


if __name__ == "__main__":
    main()
