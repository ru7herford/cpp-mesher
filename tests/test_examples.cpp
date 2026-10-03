// Checks all five steps, end to end, on the five example shapes.
//
// Runs the example inputs through the full pipeline and checks each mesh
// against its input geometry (tests/mesh_checks.hpp), including the exact
// constrained Delaunay check of every interior edge. Then checks the recorded
// history: replaying it must rebuild every mesh exactly, and recording it must
// not change the result. Last, bad inputs must be refused with a clear message.

#include "cdt2d.hpp"
#include "io.hpp"
#include "mesh_checks.hpp"
#include "quality.hpp"
#include "refine.hpp"

#include <cstdio>
#include <fstream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

using mesher::Point;
using mesher::Tri;

namespace {

int failures = 0;
int recovered_segments = 0;  // across all examples, so the recovery stage stays tested

void check(bool condition, const std::string& name, const std::string& what) {
    if (!condition) {
        std::printf("FAIL %s: %s\n", name.c_str(), what.c_str());
        ++failures;
    }
}

/// Applies a recorded change to a set of live triangles; every removed triangle
/// must be live and every added one new.
void apply(std::set<Tri>& live, const mesher::Delta& change, const std::string& name) {
    for (const Tri& t : change.removed) check(live.erase(t) == 1, name, "replay removes a triangle that is not there");
    for (const Tri& t : change.added) check(live.insert(t).second, name, "replay adds a triangle twice");
}

void check_replay(const std::string& name, const mesher::Input& input, const std::vector<Point>& points,
                  const std::vector<Tri>& triangles) {
    const mesher::Boundary boundary = mesher::build_boundary(input.domain, input.settings.max_edge);
    mesher::BuildTrace build;
    mesher::Cdt cdt(boundary.points, boundary.segments, &build);
    const std::vector<Tri> initial = cdt.triangles();
    std::vector<mesher::Step> steps;
    mesher::refine(cdt, boundary, input.settings, &steps);
    const std::set<Tri> final_set(triangles.begin(), triangles.end());
    const auto recorded = cdt.triangles();
    check(cdt.points() == points && std::set<Tri>(recorded.begin(), recorded.end()) == final_set, name,
          "recording the history changed the mesh");

    // Construction starts from the enclosing triangle, whose corners follow the input points.
    const int n = static_cast<int>(boundary.points.size());
    std::set<Tri> live{{n, n + 1, n + 2}};
    for (const mesher::Delta& insertion : build.insertions) apply(live, insertion, name);
    for (const auto& recovery : build.recoveries) apply(live, recovery.flips, name);
    apply(live, build.exterior, name);
    apply(live, build.legalize, name);
    check(live == std::set<Tri>(initial.begin(), initial.end()), name, "construction replay misses the initial mesh");
    for (const mesher::Step& step : steps) apply(live, step, name);
    check(live == final_set, name, "refinement replay misses the final mesh");
    recovered_segments += static_cast<int>(build.recoveries.size());
}

struct Totals {
    int examples = 0;
    std::size_t triangles = 0;
    std::size_t interior_edges = 0;
    std::size_t not_delaunay = 0;
};

void run(const std::string& path, Totals& totals) {
    const std::string name = path.substr(path.find_last_of('/') + 1);
    const mesher::Input input = mesher::read_input(path);
    const mesher::Boundary boundary = mesher::build_boundary(input.domain, input.settings.max_edge);
    mesher::Cdt cdt(boundary.points, boundary.segments);
    const mesher::RefineStats stats = mesher::refine(cdt, boundary, input.settings, nullptr);

    const checks::Result r = checks::check_mesh(input, cdt, stats);
    for (const std::string& what : r.failures) check(false, name, what);
    check_replay(name, input, cdt.points(), cdt.triangles());

    ++totals.examples;
    totals.triangles += r.triangles;
    totals.interior_edges += r.interior_edges;
    totals.not_delaunay += r.not_delaunay;
    std::printf("%-18s %5zu points %5zu triangles  min angle %6.3f  interior edges %5zu, not Delaunay %zu\n",
                name.c_str(), cdt.points().size(), r.triangles, r.min_angle, r.interior_edges, r.not_delaunay);
}

/// Bad inputs must stop with an error whose message names the problem.
int check_refusals(const std::string& scratch) {
    struct Case {
        const char* json;
        const char* message;  // a piece of the expected error message
    };
    const Case cases[] = {
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "min_angle_deg": 30, "max_egde": 0.1})", "unknown input key \"max_egde\""},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "min_angle_deg": 30, "min_angle_deg": 20})", "appears twice"},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1],[0,0]], "min_angle_deg": 30})", "ends with its first point"},
        {R"({"outer": [[0,0],[1,0],[1,0],[1,1],[0,1]], "min_angle_deg": 30})", "same point twice in a row"},
        {R"({"outer": [[0,0],[3,0],[0,1],[1,2]], "min_angle_deg": 30})", "polygon sides intersect"},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "holes": [[[2,2],[3,2],[3,3]]], "min_angle_deg": 30})", "every hole must lie inside"},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "min_angle_deg": 60})", "between 0 and 60"},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "min_angle_deg": +30})", "expected a value"},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "min_angle_deg": 30, "max_edge": .5})", "expected a value"},
        {R"({"outer": [[0,0],[1,0],[1,1],[0,1]], "min_angle_deg": 30, "max_edge": -1})", "must be positive"},
    };
    int refused = 0;
    for (const Case& c : cases) {
        const std::string path = scratch + "/refusal_input.json";
        std::ofstream(path) << c.json;
        std::string message = "(accepted)";
        try {
            const mesher::Input input = mesher::read_input(path);
            mesher::build_boundary(input.domain, input.settings.max_edge);
        } catch (const std::exception& error) {
            message = error.what();
        }
        const bool ok = message.find(c.message) != std::string::npos;
        check(ok, "refusals", std::string(c.json) + " gave: " + message);
        refused += ok ? 1 : 0;
    }
    std::remove((scratch + "/refusal_input.json").c_str());
    return refused;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "examples";
    const std::string scratch = argc > 2 ? argv[2] : ".";
    Totals totals;
    for (const char* example :
         {"airfoil.json", "slot.json", "sharp_corner.json", "square_hole.json", "slanted_hole.json"}) {
        try {
            run(dir + "/" + example, totals);
        } catch (const std::exception& error) {
            std::printf("FAIL %s: %s\n", example, error.what());
            ++failures;
        }
    }
    check(recovered_segments > 0, "examples", "no example exercises segment recovery");
    const int refused = check_refusals(scratch);
    std::printf("refused %d bad inputs with the expected message\n", refused);
    std::printf("SUMMARY {\"test\": \"examples\", \"examples\": %d, \"triangles\": %zu, \"interior_edges\": %zu, "
                "\"not_delaunay\": %zu, \"recovered_segments\": %d, \"bad_inputs_refused\": %d, \"failures\": %d}\n",
                totals.examples, totals.triangles, totals.interior_edges, totals.not_delaunay, recovered_segments,
                refused, failures);
    std::printf(failures == 0 ? "all checks passed\n" : "%d checks failed\n", failures);
    return failures == 0 ? 0 : 1;
}
