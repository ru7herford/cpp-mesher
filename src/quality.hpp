// How good the triangles are: angles, edge lengths, areas. Step 5 uses the
// angles to find thin triangles; the report uses all of it.
#pragma once

#include "cdt2d.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace mesher {

/// Width of one bin of the report's angle histogram, and the number of bins that cover 0 to 180 degrees.
inline constexpr double kAngleHistogramBinDeg = 10.0;
inline constexpr std::size_t kAngleHistogramBins = 18;
static_assert(kAngleHistogramBins * kAngleHistogramBinDeg == 180.0);

/// Straight-line distance between two points.
double distance(const Point& a, const Point& b);

/// Interior angles in degrees, at t[0], t[1], t[2].
std::array<double, 3> triangle_angles(const std::vector<Point>& points, const Tri& t);
/// Area; positive for a counter-clockwise triangle.
double triangle_area(const std::vector<Point>& points, const Tri& t);
/// Length of the longest of the three edges.
double longest_edge(const std::vector<Point>& points, const Tri& t);

/// Summary of a whole mesh, for the report.
struct Quality {
    double min_angle = 0.0;  // degrees
    double max_angle = 0.0;
    // The quality number: circumradius / shortest edge = 1 / (2 sin(smallest
    // angle)). 0.577 for an equilateral triangle, 1 at 30 degrees, bigger is thinner.
    // This is the largest over the mesh.
    double max_radius_edge = 0.0;
    double min_area = 0.0;
    double max_area = 0.0;
    double total_area = 0.0;
    double max_edge = 0.0;
    std::array<int, kAngleHistogramBins> angle_histogram{};  // every angle, in kAngleHistogramBinDeg bins
};

/// Measure every triangle of a finished mesh.
Quality measure(const std::vector<Point>& points, const std::vector<Tri>& triangles);

}  // namespace mesher
