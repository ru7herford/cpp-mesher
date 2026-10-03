// Shared by test_examples.cpp and test_fuzz.cpp: the checks run on every finished mesh.
//
// Checks a finished mesh against its input geometry, without trusting the
// mesher's own bookkeeping: orientation, edge sharing, boundary, area,
// topology, the constrained Delaunay property (with exact arithmetic) and
// where triangles below the angle target are allowed to stay.
//
// The polygon helpers below repeat small pieces of src/ on purpose, so that a
// bug there cannot hide itself here.

#pragma once

#include "cdt2d.hpp"
#include "exact.hpp"
#include "io.hpp"
#include "quality.hpp"
#include "refine.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <string>
#include <vector>

namespace checks {

using mesher::Edge;
using mesher::Point;
using mesher::Tri;

inline double signed_area(const std::vector<Point>& loop) {
    double twice = 0.0;
    for (std::size_t i = 1; i + 1 < loop.size(); ++i) {
        twice += (loop[i][0] - loop[0][0]) * (loop[i + 1][1] - loop[0][1]) -
                 (loop[i][1] - loop[0][1]) * (loop[i + 1][0] - loop[0][0]);
    }
    return 0.5 * twice;
}

inline bool inside(const std::vector<Point>& loop, const Point& p) {
    bool result = false;
    for (std::size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++) {
        const Point& a = loop[i];
        const Point& b = loop[j];
        if ((a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) {
            result = !result;
        }
    }
    return result;
}

/// Parameter of p along a->b, or -1 when p is not on the segment (within tol).
inline double on_side(const Point& a, const Point& b, const Point& p, double tol) {
    const double dx = b[0] - a[0], dy = b[1] - a[1];
    const double length2 = dx * dx + dy * dy;
    const double t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / length2;
    const double off = std::abs((p[0] - a[0]) * dy - (p[1] - a[1]) * dx) / std::sqrt(length2);
    return (t >= -1e-12 && t <= 1 + 1e-12 && off <= tol) ? t : -1.0;
}

struct Corner {
    Point apex, next, previous;
    double angle;  // degrees, measured through the domain
};

/// Every input corner, with the domain on the left of each side.
inline std::vector<Corner> corners(const mesher::Domain& domain) {
    std::vector<std::vector<Point>> loops{domain.outer};
    loops.insert(loops.end(), domain.holes.begin(), domain.holes.end());
    std::vector<Corner> result;
    for (std::size_t l = 0; l < loops.size(); ++l) {
        auto loop = loops[l];
        if ((l == 0) != (signed_area(loop) > 0.0)) std::reverse(loop.begin(), loop.end());
        for (std::size_t k = 0; k < loop.size(); ++k) {
            const Point& c = loop[k];
            const Point& n = loop[(k + 1) % loop.size()];
            const Point& p = loop[(k + loop.size() - 1) % loop.size()];
            double angle = std::atan2((n[0] - c[0]) * (p[1] - c[1]) - (n[1] - c[1]) * (p[0] - c[0]),
                                      (n[0] - c[0]) * (p[0] - c[0]) + (n[1] - c[1]) * (p[1] - c[1]));
            angle = angle * 180.0 / std::numbers::pi;
            if (angle < 0.0) angle += 360.0;
            result.push_back({c, n, p, angle});
        }
    }
    return result;
}

struct Result {
    std::vector<std::string> failures;
    std::size_t triangles = 0;
    std::size_t interior_edges = 0;  // edges that are not boundary segments: each one gets the exact in-circle test
    std::size_t not_delaunay = 0;
    double min_angle = 180.0;            // smallest angle in the mesh
    double smallest_input_angle = 360.0;  // smallest corner of the input polygon
    int below_target = 0;                 // triangles whose smallest angle is under the target
    double sharpest_allowed = 0.0;        // the smaller of the target and the smallest input corner
};

inline Result check_mesh(const mesher::Input& input, const mesher::Cdt& cdt, const mesher::RefineStats& stats) {
    Result r;
    const auto fail = [&r](const std::string& what) { r.failures.push_back(what); };
    const auto& points = cdt.points();
    const std::vector<Tri> triangles = cdt.triangles();
    const double target = input.settings.min_angle_deg;
    r.triangles = triangles.size();

    // Orientation (exact) and edge sharing: each directed edge once, each undirected edge in one or two triangles.
    std::map<std::pair<int, int>, int> apex;  // directed edge -> opposite vertex
    for (const Tri& t : triangles) {
        if (exact::orient(points[t[0]].data(), points[t[1]].data(), points[t[2]].data()) <= 0) {
            fail("a triangle is not strictly counter-clockwise");
        }
        for (int i = 0; i < 3; ++i) {
            if (!apex.emplace(std::pair{t[i], t[(i + 1) % 3]}, t[(i + 2) % 3]).second) {
                fail("a directed edge is used by two triangles");
            }
        }
    }
    std::set<Edge> edges, boundary_edges;
    for (const auto& [e, w] : apex) {
        const Edge key = mesher::edge_key(e.first, e.second);
        edges.insert(key);
        const auto twin = apex.find({e.second, e.first});
        if (twin == apex.end()) {
            boundary_edges.insert(key);
            continue;
        }
        if (e.first > e.second || cdt.is_segment(e.first, e.second)) continue;
        // The constrained Delaunay property, edge by edge: the vertex across a non-segment edge
        // must not lie strictly inside the circumcircle of the triangle on this side.
        ++r.interior_edges;
        if (exact::incircle(points[e.first].data(), points[e.second].data(), points[w].data(),
                            points[twin->second].data()) > 0) {
            ++r.not_delaunay;
        }
    }
    if (r.not_delaunay > 0) fail(std::to_string(r.not_delaunay) + " edges are not Delaunay");
    if (boundary_edges != cdt.segments()) fail("the mesh boundary is not exactly the set of boundary segments");

    // Topology: a polygon with h holes has Euler characteristic V - E + F = 1 - h.
    const long euler = static_cast<long>(points.size()) - static_cast<long>(edges.size()) +
                       static_cast<long>(triangles.size());
    if (euler != 1 - static_cast<long>(input.domain.holes.size())) fail("V - E + F does not match the number of holes");

    // Every input side is covered by a chain of mesh edges.
    std::vector<std::vector<Point>> loops{input.domain.outer};
    loops.insert(loops.end(), input.domain.holes.begin(), input.domain.holes.end());
    double span = 0.0;
    for (const Point& p : input.domain.outer) span = std::max({span, std::abs(p[0]), std::abs(p[1])});
    const double tol = 1e-9 * span;
    for (const auto& loop : loops) {
        for (std::size_t k = 0; k < loop.size(); ++k) {
            const Point& a = loop[k];
            const Point& b = loop[(k + 1) % loop.size()];
            std::vector<std::pair<double, int>> chain;
            for (std::size_t v = 0; v < points.size(); ++v) {
                const double t = on_side(a, b, points[v], tol);
                if (t >= 0.0) chain.emplace_back(t, static_cast<int>(v));
            }
            std::sort(chain.begin(), chain.end());
            bool covered = chain.size() >= 2 && points[chain.front().second] == a && points[chain.back().second] == b;
            for (std::size_t i = 0; covered && i + 1 < chain.size(); ++i) {
                covered = edges.count(mesher::edge_key(chain[i].second, chain[i + 1].second)) > 0;
            }
            if (!covered) fail("an input side is not a chain of mesh edges");
        }
    }

    // No triangle in a hole or outside, and the areas add up.
    double area = 0.0;
    for (const Tri& t : triangles) {
        const Point c{(points[t[0]][0] + points[t[1]][0] + points[t[2]][0]) / 3.0,
                      (points[t[0]][1] + points[t[1]][1] + points[t[2]][1]) / 3.0};
        bool in_hole = false;
        for (const auto& hole : input.domain.holes) in_hole = in_hole || inside(hole, c);
        if (!inside(input.domain.outer, c) || in_hole) fail("a triangle lies outside the domain");
        area += mesher::triangle_area(points, t);
    }
    double expected = std::abs(signed_area(input.domain.outer));
    for (const auto& hole : input.domain.holes) expected -= std::abs(signed_area(hole));
    if (std::abs(area - expected) > 1e-9 * expected) fail("mesh area differs from the domain area");

    // Triangles under the target may stay only next to an input corner sharper than 60 degrees:
    // touching a corner sharper than the target, or with their shortest edge joining the two
    // sides of a corner sharper than 60 degrees (seditious edges). Nowhere else.
    constexpr double kSeditiousCornerDeg = 60.0;
    const std::vector<Corner> all = corners(input.domain);
    std::vector<Corner> sharp;
    for (const Corner& c : all) {
        r.smallest_input_angle = std::min(r.smallest_input_angle, c.angle);
        if (c.angle < kSeditiousCornerDeg) sharp.push_back(c);
    }
    int unexplained = 0;
    for (const Tri& t : triangles) {
        const auto angles = mesher::triangle_angles(points, t);
        const int k = static_cast<int>(std::min_element(angles.begin(), angles.end()) - angles.begin());
        r.min_angle = std::min(r.min_angle, angles[k]);
        if (angles[k] >= target) continue;
        ++r.below_target;
        const Point& p = points[t[(k + 1) % 3]];
        const Point& q = points[t[(k + 2) % 3]];
        bool explained = false;
        for (const Corner& c : sharp) {
            const bool touches = c.angle < target &&
                                 (points[t[0]] == c.apex || points[t[1]] == c.apex || points[t[2]] == c.apex);
            const bool spans = (on_side(c.apex, c.next, p, tol) >= 0 && on_side(c.apex, c.previous, q, tol) >= 0) ||
                               (on_side(c.apex, c.previous, p, tol) >= 0 && on_side(c.apex, c.next, q, tol) >= 0);
            explained = explained || touches || spans;
        }
        unexplained += explained ? 0 : 1;
    }
    if (unexplained > 0) fail(std::to_string(unexplained) + " triangles below the angle target away from sharp corners");
    r.sharpest_allowed = std::min(target, r.smallest_input_angle);
    if (stats.hit_point_limit) fail("refinement hit the point limit");
    if (!stats.targets_met()) fail("refinement reports unmet targets");
    if (input.settings.max_edge > 0.0 && mesher::measure(points, triangles).max_edge > input.settings.max_edge) {
        fail("an edge is longer than max_edge");
    }
    return r;
}

}  // namespace checks
