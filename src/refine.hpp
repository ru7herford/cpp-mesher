// Step 5: fix thin triangles until none are left (Ruppert's algorithm).
// Also step 1: checking the input polygon and placing points along its sides.
#pragma once

#include "cdt2d.hpp"

#include <vector>

namespace mesher {

/// A polygon with holes. Loop orientation does not matter.
struct Domain {
    std::vector<Point> outer;
    std::vector<std::vector<Point>> holes;
};

/// The targets from the input file.
struct Settings {
    double min_angle_deg = 30.0;
    double max_edge = 0.0;     // 0 means no size limit
    int max_points = 100000;  // stop here, so a target that cannot be reached still ends
};

/// The boundary points and edges for steps 2 to 4, plus what step 5 needs to
/// know about them. A "side" is one side of the input polygon, from corner to
/// corner; it may be cut into several segments.
struct Boundary {
    std::vector<Point> points;
    std::vector<Edge> segments;
    std::vector<Edge> sides;            // the input polygon edges, between corner points
    std::vector<int> side_of_point;     // which side each point lies on; -1 for corners
    std::vector<double> corner_angle;   // angle inside the domain at each corner, in degrees; 180 elsewhere
};

/// Area inside a loop: positive when its points go counter-clockwise.
double signed_area(const std::vector<Point>& loop);

/// Check the input polygon (sides must not cross, holes must be inside), turn
/// every loop the same way round, and cut each side into equal pieces no
/// longer than `max_edge`.
Boundary build_boundary(const Domain& domain, double max_edge);

/// Counts for the report.
struct RefineStats {
    int inserted = 0;          // circumcentres added
    int splits = 0;            // segments split
    int refused_splits = 0;    // segments that became too short to split again
    int rejected_points = 0;   // circumcentres that could not be added
    bool hit_point_limit = false;
    // Final state, counted after refinement stops.
    int below_min_angle = 0;   // triangles under the angle target
    int at_sharp_corner = 0;   // of those, ones left alone at a sharp input corner (see the sharp-corner rule in assess)
    int over_max_edge = 0;
    bool targets_met() const { return below_min_angle == at_sharp_corner && over_max_edge == 0; }
};

/// Step 5: add points until every triangle meets the targets, or no more
/// progress is possible. When `steps` is given, every step is added to it.
RefineStats refine(Cdt& cdt, const Boundary& boundary, const Settings& settings, std::vector<Step>* steps);

}  // namespace mesher
