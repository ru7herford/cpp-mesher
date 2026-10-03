// Checks all five steps on many random shapes.
//
// Random polygons, with and without a hole, at angle targets of 20, 25 and 30
// degrees, with and without a size limit; then wedges whose tip is sharper
// than the target. Every mesh goes through tests/mesh_checks.hpp: valid,
// boundary kept, exact constrained Delaunay property, area, topology, and no
// triangle below the target except at an input corner sharper than the target.
//
// The polygons are star-shaped around the origin (corners sorted by angle,
// every gap under 180 degrees), which makes them simple by construction while
// still producing sharp and nearly straight corners.

#include "cdt2d.hpp"
#include "io.hpp"
#include "mesh_checks.hpp"
#include "refine.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <numbers>
#include <string>
#include <vector>

using mesher::Point;

namespace {

int failures = 0;

/// Deterministic and identical on every platform (std:: distributions are not).
struct Random {
    std::uint64_t state;
    std::uint64_t next() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uniform(double low, double high) { return low + (high - low) * static_cast<double>(next() >> 11) * 0x1.0p-53; }
    int integer(int low, int high) { return low + static_cast<int>(next() % static_cast<std::uint64_t>(high - low + 1)); }
};

/// Corners at sorted random angles, no gap under `min_gap` or over pi - min_gap, random radii.
std::vector<Point> star(Random& random, int corners, double radius_low, double radius_high, Point centre, double scale) {
    constexpr double kMinGap = 0.15;  // radians
    const double tau = 2.0 * std::numbers::pi;
    std::vector<double> angle;
    while (true) {
        angle.clear();
        for (int k = 0; k < corners; ++k) angle.push_back(random.uniform(0.0, tau));
        std::sort(angle.begin(), angle.end());
        bool ok = true;
        for (int k = 0; k < corners; ++k) {
            const double gap = k + 1 < corners ? angle[k + 1] - angle[k] : angle[0] + tau - angle[k];
            ok = ok && gap > kMinGap && gap < std::numbers::pi - kMinGap;
        }
        if (ok) break;
    }
    std::vector<Point> loop;
    for (const double theta : angle) {
        const double r = scale * random.uniform(radius_low, radius_high);
        loop.push_back({centre[0] + r * std::cos(theta), centre[1] + r * std::sin(theta)});
    }
    return loop;
}

/// Distance from c to the nearest side of the loop.
double clearance(const std::vector<Point>& loop, const Point& c) {
    double best = INFINITY;
    for (std::size_t k = 0; k < loop.size(); ++k) {
        const Point& a = loop[k];
        const Point& b = loop[(k + 1) % loop.size()];
        const double dx = b[0] - a[0], dy = b[1] - a[1];
        const double t = std::clamp(((c[0] - a[0]) * dx + (c[1] - a[1]) * dy) / (dx * dx + dy * dy), 0.0, 1.0);
        best = std::min(best, std::hypot(a[0] + t * dx - c[0], a[1] + t * dy - c[1]));
    }
    return best;
}

struct Totals {
    int meshes = 0, with_holes = 0, with_size_limit = 0, with_small_angles = 0;
    long points = 0, triangles = 0, interior_edges = 0, not_delaunay = 0;
    double worst_ratio = INFINITY;  // smallest (mesh min angle / smaller of target and sharpest input corner)
    std::string worst_case;
    int max_points = 0;
};

/// The input as mesher JSON, printed with a failure so the case can be rerun on its own.
std::string to_json(const mesher::Input& input) {
    const auto loop = [](const std::vector<Point>& points) {
        std::string text = "[";
        for (std::size_t i = 0; i < points.size(); ++i) {
            char pair[64];
            std::snprintf(pair, sizeof pair, "%s[%.17g, %.17g]", i ? ", " : "", points[i][0], points[i][1]);
            text += pair;
        }
        return text + "]";
    };
    std::string text = "{\"outer\": " + loop(input.domain.outer);
    if (!input.domain.holes.empty()) text += ", \"holes\": [" + loop(input.domain.holes[0]) + "]";
    char settings[96];
    std::snprintf(settings, sizeof settings, ", \"min_angle_deg\": %.17g", input.settings.min_angle_deg);
    text += settings;
    if (input.settings.max_edge > 0.0) {
        std::snprintf(settings, sizeof settings, ", \"max_edge\": %.17g", input.settings.max_edge);
        text += settings;
    }
    return text + "}";
}

/// Returns the mesh's point count and smallest angle (0, 0 when it could not be meshed).
std::pair<std::size_t, double> mesh_and_check(const std::string& name, const mesher::Input& input, Totals& totals) {
    try {
        const mesher::Boundary boundary = mesher::build_boundary(input.domain, input.settings.max_edge);
        mesher::Cdt cdt(boundary.points, boundary.segments);
        const mesher::RefineStats stats = mesher::refine(cdt, boundary, input.settings, nullptr);
        const checks::Result r = checks::check_mesh(input, cdt, stats);
        for (const std::string& what : r.failures) {
            std::printf("FAIL %s: %s\n", name.c_str(), what.c_str());
            ++failures;
        }
        if (!r.failures.empty()) std::printf("input: %s\n", to_json(input).c_str());
        ++totals.meshes;
        totals.with_holes += input.domain.holes.empty() ? 0 : 1;
        totals.with_size_limit += input.settings.max_edge > 0.0 ? 1 : 0;
        totals.points += static_cast<long>(cdt.points().size());
        totals.max_points = std::max(totals.max_points, static_cast<int>(cdt.points().size()));
        totals.triangles += static_cast<long>(r.triangles);
        totals.interior_edges += static_cast<long>(r.interior_edges);
        totals.not_delaunay += static_cast<long>(r.not_delaunay);
        if (r.below_target > 0) ++totals.with_small_angles;
        // How far below what the input allows (the target, or a sharper input corner) the mesh went.
        const double ratio = r.min_angle / r.sharpest_allowed;
        if (ratio < totals.worst_ratio) {
            totals.worst_ratio = ratio;
            char text[200];
            std::snprintf(text, sizeof text, "%s: smallest angle %.2f deg, target %.0f deg, sharpest corner %.2f deg",
                          name.c_str(), r.min_angle, input.settings.min_angle_deg, r.smallest_input_angle);
            totals.worst_case = text;
        }
        return {cdt.points().size(), r.min_angle};
    } catch (const std::exception& error) {
        std::printf("FAIL %s: %s\n", name.c_str(), error.what());
        ++failures;
    }
    return {0, 0.0};
}

}  // namespace

int main(int argc, char** argv) {
    const int polygons = argc > 1 ? std::stoi(argv[1]) : 300;
    const auto start = std::chrono::steady_clock::now();
    Random random{20261002};
    Totals totals;

    for (int n = 0; n < polygons; ++n) {
        mesher::Input input;
        const double scale = std::ldexp(1.0, random.integer(-3, 6));
        const Point centre{random.uniform(-10, 10) * scale, random.uniform(-10, 10) * scale};
        input.domain.outer = star(random, random.integer(3, 12), 0.35, 1.0, centre, scale);
        if (n % 2 == 1) {
            const double room = clearance(input.domain.outer, centre);
            input.domain.holes.push_back(star(random, random.integer(3, 6), 0.3, 0.7, centre, room));
        }
        input.settings.min_angle_deg = 20.0 + 5.0 * (n % 3);
        if (n % 4 >= 2) input.settings.max_edge = scale * random.uniform(0.1, 0.4);
        mesh_and_check("polygon " + std::to_string(n), input, totals);
    }

    std::string wedges;  // [tip, max_edge, points, smallest angle] for the summary
    // Wedges: a tip of 0.5 to 15 degrees, with and without a size limit, at the 30 degree target.
    for (const double tip : {0.5, 1.0, 2.0, 5.0, 10.0, 15.0}) {
        for (const double max_edge : {0.0, 0.3}) {
            mesher::Input input;
            const double a = tip * std::numbers::pi / 180.0;
            const Point b{4.0 * std::cos(a), 4.0 * std::sin(a)};
            input.domain.outer = {{0, 0}, {4, 0}, {4, b[1] + 2}, {b[0], b[1] + 2}, b};
            input.settings.min_angle_deg = 30.0;
            input.settings.max_edge = max_edge;
            char name[64];
            std::snprintf(name, sizeof name, "wedge %.1f deg%s", tip, max_edge > 0.0 ? ", max_edge 0.3" : "");
            const auto [points, min_angle] = mesh_and_check(name, input, totals);
            std::printf("%-24s %6zu points, smallest angle %.4f deg\n", name, points, min_angle);
            char row[96];
            std::snprintf(row, sizeof row, "%s[%g, %g, %zu, %.4f]", wedges.empty() ? "" : ", ", tip, max_edge, points, min_angle);
            wedges += row;
        }
    }

    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::printf("%d meshes (%d with a hole, %d with a size limit, %d keep triangles under the target)\n", totals.meshes,
                totals.with_holes, totals.with_size_limit, totals.with_small_angles);
    std::printf("%ld points, %ld triangles, largest mesh %d points; %ld interior edges, %ld not Delaunay\n",
                totals.points, totals.triangles, totals.max_points, totals.interior_edges, totals.not_delaunay);
    std::printf("worst smallest angle: %.4f of the target or sharpest corner (%s)\n", totals.worst_ratio,
                totals.worst_case.c_str());
    std::printf("SUMMARY {\"test\": \"fuzz\", \"polygons\": %d, \"wedges\": %d, \"meshes\": %d, \"with_holes\": %d, "
                "\"with_size_limit\": %d, \"with_small_angles\": %d, \"points\": %ld, \"triangles\": %ld, "
                "\"largest_mesh_points\": %d, \"interior_edges\": %ld, \"not_delaunay\": %ld, \"worst_ratio\": %.6f, "
                "\"worst_case\": \"%s\", \"wedges_tip_maxedge_points_minangle\": [%s], \"seconds\": %.2f, \"failures\": %d}\n",
                polygons, totals.meshes - polygons, totals.meshes, totals.with_holes, totals.with_size_limit,
                totals.with_small_angles, totals.points, totals.triangles, totals.max_points, totals.interior_edges,
                totals.not_delaunay, totals.worst_ratio, totals.worst_case.c_str(), wedges.c_str(), seconds, failures);
    std::printf(failures == 0 ? "all checks passed\n" : "%d checks failed\n", failures);
    return failures == 0 ? 0 : 1;
}
