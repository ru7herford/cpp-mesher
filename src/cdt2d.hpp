// The mesh itself: triangles, and the boundary edges they must keep.
// Steps 2, 3 and 4 happen in the constructor; step 5 (refine.cpp) then calls
// insert() and split() to add points one at a time.
//
// Words used below:
//   segment  an edge of the shape that the mesh must keep. It starts as a side
//            (or a piece of a side) of the input polygon and can be split.
//   domain   the inside of the polygon, minus its holes.
//   cavity   the hole a new point makes: the triangles whose circle contains it.
//   apex     the corner of a triangle opposite a given edge.
//
// Constructor, in order:
//   1. connect all boundary points into triangles (Bowyer-Watson);
//   2. force in each segment that is missing, by edge flips (Sloan);
//   3. delete every triangle outside the domain, then flip edges until the
//      mesh is a constrained Delaunay triangulation again (Lawson flips).
// Every change reports exactly which triangles it removed and added, so a
// run can be recorded and replayed on the web page.

#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mesher {

using Point = std::array<double, 2>;
using Tri = std::array<int, 3>;   // three point indices, counter-clockwise
using Edge = std::pair<int, int>;  // two point indices, smaller one first

/// The same edge whichever way round it is given: smaller index first.
Edge edge_key(int a, int b);

/// Lower-left and upper-right corners of the points' bounding box.
std::pair<Point, Point> bounding_box(const std::vector<Point>& points);

/// What one change did to the mesh: the triangles it took away and the ones
/// it put in. Triangles that were made and removed inside the same change
/// are left out.
struct Delta {
    std::vector<Tri> removed;  // triangles that existed before the action
    std::vector<Tri> added;    // triangles that exist after it
};

/// Why step 5 added a point. Only written to the steps file for the web page;
/// the mesher never reads it back.
struct Cause {
    enum class Reason {
        encroached,         // a point sits inside the circle that has a segment as its diameter
        small_angle,        // a triangle is thin: its smallest angle is under the target
        too_long,           // a triangle has an edge longer than max_edge
        centre_encroaches,  // a triangle's circumcentre would encroach on a segment
        centre_outside,     // ... would land outside the domain or in a hole
        centre_on_segment,  // ... lands exactly on a segment
    };
    Reason reason = Reason::encroached;
    std::optional<Tri> triangle;    // the thin (or too long) triangle behind the step
    std::optional<Edge> segment;    // the segment that was split
    std::optional<int> encroacher;  // the point inside the segment's circle
    std::optional<Point> centre;    // the circumcentre that was not added
    bool shell = false;             // split at a power-of-two distance from an input corner
};

/// One step of step 5: a new point and the triangles it changed.
struct Step : Delta {
    enum class Kind { insert, split };
    Kind kind = Kind::insert;
    int point = -1;  // index of the new point
    // insert: point added inside the mesh. split: point added on a segment.
    Cause cause;
};

/// Steps 2 to 4, change by change, for the steps file. Only filled when the
/// constructor is given one; otherwise nothing is recorded.
struct BuildTrace {
    struct Recovery {
        Edge segment;
        Delta flips;
    };
    std::array<Point, 3> enclosing{};  // corners of the big starting triangle (point indices n, n+1, n+2)
    std::vector<Delta> insertions;     // step 2: one change per boundary point, in index order
    std::vector<Recovery> recoveries;  // step 3: only the segments that were missing
    Delta exterior;                    // step 4: triangles outside the domain or in holes, all removed
    Delta legalize;                    // step 4: Lawson flips that make the mesh Delaunay again
};

class Cdt {
public:
    /// Runs steps 2 to 4. `segments` must form closed loops (the outer
    /// boundary and the holes). A triangle is kept when you cross an odd
    /// number of segments to reach it from far outside. `trace`, when given,
    /// records every change.
    Cdt(std::vector<Point> points, const std::vector<Edge>& segments, BuildTrace* trace = nullptr);

    struct InsertResult {
        // inserted: done. on_segment: p is exactly on a segment, split that instead.
        // rejected: p lies outside the domain or on an existing point.
        enum class Status { inserted, on_segment, rejected };
        Status status;
        Step step;     // valid when inserted
        Edge segment;  // valid when on_segment: split this segment instead
    };

    /// Add point `p` inside the mesh (Bowyer-Watson, see cavity()). `near` is
    /// any current triangle close to `p`; it only makes finding `p` faster.
    InsertResult insert(const Point& p, const Tri& near);

    /// Add point `p` on a segment, which becomes two segments. `p` must lie
    /// between the segment's ends. Edge flips then make the mesh Delaunay again.
    Step split(const Edge& segment, const Point& p);

    const std::vector<Point>& points() const { return points_; }
    const std::set<Edge>& segments() const { return segments_; }
    std::vector<Tri> triangles() const;     // every current triangle
    bool contains(const Tri& t) const;      // whether t is still in the mesh
    bool is_segment(int a, int b) const { return segments_.count(edge_key(a, b)) > 0; }
    /// For each triangle (one or two) that has `segment` as a side, its third point.
    std::vector<int> opposite_vertices(const Edge& segment) const;

private:
    // The hole a new point makes in Bowyer-Watson: the triangles to remove,
    // and the edges around them, which get joined to the new point.
    struct Cavity {
        std::vector<int> triangles;
        std::vector<Edge> boundary;  // each edge in counter-clockwise order around the hole
    };

    // A triangle lives in a numbered slot of tris_; removed slots are reused.
    int add(const Tri& t);
    void remove(int slot);
    int triangle_at(int a, int b) const;  // the triangle with edge a->b (counter-clockwise), or -1
    int apex(int slot, int a) const;      // its third point, across from the edge that starts at a
    int locate(const Point& p, int seed) const;
    bool triangle_contains(int slot, const Point& p) const;
    Cavity cavity(int point, int located) const;
    int fill_cavity(const Cavity& cavity, int point);
    Edge flip(int a, int b);
    bool proper_crossing(int a, int b, int c, int d) const;
    void recover_segment(const Edge& segment);
    void remove_exterior(int real_point_count);
    void legalize(std::vector<Edge> stack);

    std::vector<Point> points_;
    std::vector<Tri> tris_;
    std::vector<char> alive_;
    std::vector<int> free_slots_;
    std::unordered_map<std::uint64_t, int> half_edges_;  // edge a->b -> the triangle that has it
    std::set<Edge> segments_;
    Delta* recording_ = nullptr;  // where the current change is written down, if anywhere
};

}  // namespace mesher
