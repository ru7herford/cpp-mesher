// Angles, lengths and areas of triangles. Used by step 5 (is a triangle
// thin?) and by the report. Plain doubles: these are measurements, not
// yes/no decisions that could tangle the mesh, so rounding does no harm.

#include "quality.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace mesher {

double distance(const Point& a, const Point& b) { return std::hypot(b[0] - a[0], b[1] - a[1]); }

std::array<double, 3> triangle_angles(const std::vector<Point>& points, const Tri& t) {
    std::array<double, 3> angles{};
    for (int i = 0; i < 3; ++i) {
        const Point& centre = points[t[i]];
        const Point& next = points[t[(i + 1) % 3]];
        const Point& previous = points[t[(i + 2) % 3]];
        const double ux = next[0] - centre[0], uy = next[1] - centre[1];
        const double vx = previous[0] - centre[0], vy = previous[1] - centre[1];
        // Angle between the two edges u, v that leave this corner:
        // atan2(|u x v|, u . v). It stays accurate for angles near 0 and 180
        // degrees, where acos(u . v / (|u| |v|)) does not.
        angles[i] = std::atan2(std::abs(ux * vy - uy * vx), ux * vx + uy * vy) * 180.0 / std::numbers::pi;
    }
    return angles;
}

double triangle_area(const std::vector<Point>& points, const Tri& t) {
    const Point& a = points[t[0]];
    const Point& b = points[t[1]];
    const Point& c = points[t[2]];
    return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
}

double longest_edge(const std::vector<Point>& points, const Tri& t) {
    return std::max({distance(points[t[0]], points[t[1]]), distance(points[t[1]], points[t[2]]),
                     distance(points[t[2]], points[t[0]])});
}

Quality measure(const std::vector<Point>& points, const std::vector<Tri>& triangles) {
    Quality q;
    if (triangles.empty()) return q;
    q.min_angle = q.min_area = std::numeric_limits<double>::infinity();
    for (const Tri& t : triangles) {
        const auto angles = triangle_angles(points, t);
        for (double angle : angles) {
            q.min_angle = std::min(q.min_angle, angle);
            q.max_angle = std::max(q.max_angle, angle);
            const auto bin = static_cast<std::size_t>(angle / kAngleHistogramBinDeg);
            q.angle_histogram[std::min(kAngleHistogramBins - 1, bin)] += 1;
        }
        const double a = distance(points[t[1]], points[t[2]]);
        const double b = distance(points[t[2]], points[t[0]]);
        const double c = distance(points[t[0]], points[t[1]]);
        const double area = triangle_area(points, t);
        // Circumradius R = abc / (4 area); quality number = R / shortest edge.
        q.max_radius_edge = std::max(q.max_radius_edge, a * b * c / (4.0 * area) / std::min({a, b, c}));
        q.min_area = std::min(q.min_area, area);
        q.max_area = std::max(q.max_area, area);
        q.total_area += area;
        q.max_edge = std::max(q.max_edge, std::max({a, b, c}));
    }
    return q;
}

}  // namespace mesher
