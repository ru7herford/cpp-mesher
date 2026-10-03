// The `mesher` program: reads a polygon, runs all five steps, writes the results.
//
// usage: mesher <input.json> <out.vtk> [--report report.json] [--steps steps.json]
//
// Start reading here. The five steps are:
//   1. place points along the polygon's edges          (build_boundary, refine.cpp)
//   2. connect the boundary points into triangles      (Cdt constructor, cdt2d.cpp)
//   3. force the boundary edges into the mesh          (Cdt constructor, cdt2d.cpp)
//   4. remove the triangles outside the polygon        (Cdt constructor, cdt2d.cpp)
//   5. fix thin triangles until none are left          (refine, refine.cpp)
// This file only calls them in order, times them and writes the files.

#include "cdt2d.hpp"
#include "io.hpp"
#include "quality.hpp"
#include "refine.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>
#include <vector>

namespace {

// Milliseconds since `since`, for the timings in the report.
double elapsed_ms(std::chrono::steady_clock::time_point since) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - since).count();
}

// Wrong arguments: print how to call the program and exit with code 2.
int usage() {
    std::fprintf(stderr, "usage: mesher <input.json> <out.vtk> [--report report.json] [--steps steps.json]\n");
    return 2;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) return usage();
    const std::string input_path = argv[1], vtk_path = argv[2];
    std::string report_path, steps_path;
    for (int i = 3; i < argc; ++i) {
        const std::string flag = argv[i];
        if (i + 1 >= argc) return usage();
        if (flag == "--report") report_path = argv[++i];
        else if (flag == "--steps") steps_path = argv[++i];
        else return usage();
    }

    try {
        const mesher::Input input = mesher::read_input(input_path);
        const auto& settings = input.settings;

        // Steps 1 to 4: place points along the boundary, then build the first mesh.
        auto start = std::chrono::steady_clock::now();
        const mesher::Boundary boundary = mesher::build_boundary(input.domain, settings.max_edge);
        mesher::BuildTrace build;  // filled only when a steps file is requested
        mesher::Cdt cdt(boundary.points, boundary.segments, steps_path.empty() ? nullptr : &build);
        const double triangulate_ms = elapsed_ms(start);
        const std::vector<mesher::Tri> initial_triangles = cdt.triangles();

        // Step 5: add points until no triangle is thin (or too long, with a size limit).
        std::vector<mesher::Step> steps;
        start = std::chrono::steady_clock::now();
        const mesher::RefineStats stats =
            mesher::refine(cdt, boundary, settings, steps_path.empty() ? nullptr : &steps);
        const double refine_ms = elapsed_ms(start);

        // Measure the finished mesh and write the files.
        const auto& points = cdt.points();
        const std::vector<mesher::Tri> triangles = cdt.triangles();
        const mesher::Quality quality = mesher::measure(points, triangles);
        double domain_area = std::abs(mesher::signed_area(input.domain.outer));
        for (const auto& hole : input.domain.holes) domain_area -= std::abs(mesher::signed_area(hole));

        mesher::write_vtk(vtk_path, points, triangles);
        if (!report_path.empty()) {
            mesher::write_report(report_path, {input, points.size(), triangles.size(), cdt.segments().size(),
                                               domain_area, quality, stats, triangulate_ms, refine_ms});
        }
        if (!steps_path.empty()) {
            mesher::write_steps(steps_path, {input, boundary.points, build, initial_triangles, points, steps});
        }

        // Short summary on screen. Exit code 0 means the targets were met, 1 means they were not.
        std::printf("%zu points, %zu triangles, %zu boundary segments\n", points.size(), triangles.size(),
                    cdt.segments().size());
        std::printf("min angle %.2f deg (target %.1f), max angle %.2f deg\n", quality.min_angle,
                    settings.min_angle_deg, quality.max_angle);
        if (settings.max_edge > 0.0) {
            std::printf("longest edge %.4g (target %.4g)\n", quality.max_edge, settings.max_edge);
        }
        if (stats.targets_met()) {
            std::printf("targets met");
            if (stats.at_sharp_corner > 0) {
                std::printf(", except %d triangles next to input corners sharper than 60 deg, which refinement leaves alone",
                            stats.at_sharp_corner);
            }
            std::printf("\n");
        } else {
            std::printf("targets NOT met: %d triangles below the angle target, %d over the edge target%s%s\n",
                        stats.below_min_angle - stats.at_sharp_corner, stats.over_max_edge,
                        stats.refused_splits > 0 ? "; some segments became too short to split" : "",
                        stats.hit_point_limit ? "; point limit reached" : "");
        }
        std::printf("time: triangulation %.1f ms, refinement %.1f ms\n", triangulate_ms, refine_ms);
        return stats.targets_met() ? 0 : 1;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "mesher: %s\n", error.what());
        return 2;
    }
}
