// Reading the input file and writing the three output files. Not part of the
// five steps; main.cpp calls these before and after.
#pragma once

#include "cdt2d.hpp"
#include "quality.hpp"
#include "refine.hpp"

#include <string>
#include <vector>

namespace mesher {

/// Everything the input file says: the polygon and the targets.
struct Input {
    Domain domain;
    Settings settings;
};

/// Read {"outer": [[x,y],...], "holes": [[[x,y],...],...], "min_angle_deg": a, "max_edge": h}.
/// "holes" and "max_edge" are optional; any other key is an error.
Input read_input(const std::string& path);

/// Legacy ASCII VTK unstructured grid with a per-triangle `min_angle` cell field.
void write_vtk(const std::string& path, const std::vector<Point>& points, const std::vector<Tri>& triangles);

/// What goes into report.json.
struct Report {
    const Input& input;
    std::size_t points;
    std::size_t triangles;
    std::size_t boundary_segments;
    double domain_area;
    const Quality& quality;
    const RefineStats& stats;
    double triangulate_ms;
    double refine_ms;
};

/// Counts, angles, areas, timings: one small JSON file.
void write_report(const std::string& path, const Report& report);

/// What goes into steps.json.
struct Replay {
    const Input& input;
    const std::vector<Point>& boundary_points;  // the triangulator's input points
    const BuildTrace& build;
    const std::vector<Tri>& initial_triangles;  // the triangulation refinement starts from
    const std::vector<Point>& final_points;
    const std::vector<Step>& steps;
};

/// Every change of the run, in order, so it can be replayed exactly (the web
/// page does this): the boundary points, the triangles each change of steps 2
/// to 3 removed and added, then every point of step 5 with its reason.
void write_steps(const std::string& path, const Replay& replay);

}  // namespace mesher
