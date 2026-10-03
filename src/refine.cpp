// Step 5 (Ruppert's algorithm), and step 1: placing the boundary points.
//
// Step 5 keeps two to-do lists:
//   check_  segments that a nearby point may encroach on;
//   bad_    triangles that may be thin (or too long).
// It always empties check_ first, then takes one triangle from bad_:
//   - an encroached segment is split in two;
//   - a thin triangle gets a new point at its circumcentre, unless that point
//     would encroach on a segment or land outside; then the segment is split.
// Every change puts the new triangles and segments on the lists. It stops when
// both lists are empty, or at max_points.
// The triangles are taken in the order they were found (first in, first
// out), not worst first.

#include "refine.hpp"

#include "predicates.hpp"
#include "quality.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <numbers>
#include <optional>
#include <set>
#include <stdexcept>

namespace mesher {

double signed_area(const std::vector<Point>& loop) {
    // Shoelace formula, measured from the first point: it loses less to
    // rounding when the coordinates are large.
    const Point& o = loop[0];
    double twice = 0.0;
    for (std::size_t i = 1; i + 1 < loop.size(); ++i) {
        const Point& a = loop[i];
        const Point& b = loop[i + 1];
        twice += (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    }
    return 0.5 * twice;
}

namespace {

// Input corners sharper than this get the sharp-corner rule (see assess).
// Ruppert's proof that step 5 ends needs every input corner to be at least 60 degrees.
constexpr double kSharpCornerDeg = 60.0;
// Segments are never split into pieces shorter than this fraction of the domain size.
constexpr double kMinPieceFraction = 1e-9;
// Two points count as the same distance from a corner when their distances agree
// to this fraction. Split distances are powers of two, so real differences are much larger.
constexpr double kShellTolerance = 1e-3;

// Whether p is inside a loop: count how many sides a ray from p to the right
// crosses; odd means inside. Plain doubles are enough for this input check.
bool inside_loop(const std::vector<Point>& loop, const Point& p) {
    bool inside = false;
    for (std::size_t i = 0, j = loop.size() - 1; i < loop.size(); j = i++) {
        const Point& a = loop[i];
        const Point& b = loop[j];
        if ((a[1] > p[1]) != (b[1] > p[1]) && p[0] < (b[0] - a[0]) * (p[1] - a[1]) / (b[1] - a[1]) + a[0]) {
            inside = !inside;
        }
    }
    return inside;
}

/// Angle at corner loop[k], measured on the domain's side (the left of every
/// side), from 0 to 360 degrees. Uses atan2(cross, dot) like quality.cpp.
double corner_angle(const std::vector<Point>& loop, std::size_t k) {
    const std::size_t m = loop.size();
    const Point& c = loop[k];
    const Point& next = loop[(k + 1) % m];
    const Point& previous = loop[(k + m - 1) % m];
    const double nx = next[0] - c[0], ny = next[1] - c[1];
    const double px = previous[0] - c[0], py = previous[1] - c[1];
    double angle = std::atan2(nx * py - ny * px, nx * px + ny * py) * 180.0 / std::numbers::pi;
    if (angle < 0.0) angle += 360.0;
    return angle;
}

// Refuse a polygon whose sides cross or touch, or that repeats a point.
// Checks every pair of sides, which is fine for input polygons.
void require_simple(const std::vector<std::vector<Point>>& loops) {
    std::vector<std::pair<const Point*, const Point*>> sides;
    for (const auto& loop : loops) {
        // A side of zero length has no direction, so check this before measuring any angle.
        for (std::size_t k = 0; k < loop.size(); ++k) {
            if (loop[k] == loop[(k + 1) % loop.size()]) {
                throw std::invalid_argument(k + 1 == loop.size()
                    ? "a loop ends with its first point; list every corner once, without closing the loop"
                    : "a loop lists the same point twice in a row");
            }
        }
        for (std::size_t k = 0; k < loop.size(); ++k) {
            sides.emplace_back(&loop[k], &loop[(k + 1) % loop.size()]);
            if (corner_angle(loop, k) <= 0.0) throw std::invalid_argument("a loop doubles back on itself");
        }
    }
    for (std::size_t i = 0; i < sides.size(); ++i) {
        for (std::size_t j = i + 1; j < sides.size(); ++j) {
            const auto [a, b] = sides[i];
            const auto [c, d] = sides[j];
            if (*a == *c || *a == *d || *b == *c || *b == *d) continue;  // neighbours meet at a corner
            if (segments_intersect_2d(a->data(), b->data(), c->data(), d->data())) {
                throw std::invalid_argument("polygon sides intersect; the domain must be a simple polygon with disjoint holes");
            }
        }
    }
}

// Centre of the circle through a, b and c. With b' = b - a and c' = c - a:
//   D = 2 (b'x c'y - b'y c'x),  centre = a + (|b'|^2 c'y - |c'|^2 b'y, |c'|^2 b'x - |b'|^2 c'x) / D.
// For
// three points on one line it falls back to their average, which the caller
// then refuses or handles like any other point.
Point circumcentre(const Point& a, const Point& b, const Point& c) {
    const double bx = b[0] - a[0], by = b[1] - a[1];
    const double cx = c[0] - a[0], cy = c[1] - a[1];
    const double denominator = 2.0 * (bx * cy - by * cx);
    if (denominator == 0.0) {
        return {(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};
    }
    const double b2 = bx * bx + by * by;
    const double c2 = cx * cx + cy * cy;
    return {a[0] + (b2 * cy - c2 * by) / denominator, a[1] + (bx * c2 - cx * b2) / denominator};
}

// Step 5. One Refiner does one run.
class Refiner {
public:
    Refiner(Cdt& cdt, const Boundary& boundary, const Settings& settings, std::vector<Step>* steps)
        : cdt_(cdt), boundary_(boundary), settings_(settings), steps_(steps),
          side_of_(boundary.side_of_point) {
        const auto [low, high] = bounding_box(boundary.points);
        min_piece_ = kMinPieceFraction * std::max(high[0] - low[0], high[1] - low[1]);
    }

    // The main loop. Stopping rule: stop when both to-do lists are empty (every
    // segment checked, every triangle good or left alone), or at max_points.
    RefineStats run() {
        for (const Edge& s : cdt_.segments()) check_.push_back(s);
        for (const Tri& t : cdt_.triangles()) bad_.push_back(t);
        // Encroached segments go first, as in Ruppert's algorithm: a circumcentre
        // is only safe to add once no segment near it is encroached.
        while (true) {
            if (static_cast<int>(cdt_.points().size()) >= settings_.max_points) {
                stats_.hit_point_limit = true;
                break;
            }
            if (!check_.empty()) {
                const Edge s = check_.front();
                check_.pop_front();
                if (!cdt_.is_segment(s.first, s.second)) continue;  // already split earlier
                if (const auto v = encroacher(s)) {
                    Cause cause;
                    cause.reason = Cause::Reason::encroached;
                    cause.encroacher = *v;
                    split(s, cause);
                }
                continue;
            }
            if (bad_.empty()) break;
            const Tri t = bad_.front();
            bad_.pop_front();
            if (!cdt_.contains(t)) continue;  // already removed by an earlier step
            if (const Assessment a = assess(t); needs_refinement(a)) refine_triangle(t, a);
        }
        // Count what is left, for the report and the exit code.
        for (const Tri& t : cdt_.triangles()) {
            const Assessment a = assess(t);
            stats_.below_min_angle += a.small_angle ? 1 : 0;
            stats_.at_sharp_corner += a.small_angle && a.unfixable ? 1 : 0;
            stats_.over_max_edge += a.too_long ? 1 : 0;
        }
        return stats_;
    }

private:
    // What is wrong with one triangle.
    struct Assessment {
        bool small_angle = false;  // thin: smallest angle under the target
        bool unfixable = false;  // thin, but left alone because of a sharp input corner
        bool too_long = false;
    };

    // Whether step 5 should work on this triangle.
    static bool needs_refinement(const Assessment& a) {
        return a.too_long || (a.small_angle && !a.unfixable);
    }

    // Whether point v is a corner of the input polygon (not a point added on a side or inside).
    bool is_corner(int v) const {
        return v < static_cast<int>(boundary_.points.size()) && side_of_[v] < 0;
    }

    // Is triangle t thin or too long, and if thin, can it be fixed?
    Assessment assess(const Tri& t) const {
        const auto& points = cdt_.points();
        Assessment a;
        a.too_long = settings_.max_edge > 0.0 && longest_edge(points, t) > settings_.max_edge;
        const auto angles = triangle_angles(points, t);
        const int k = static_cast<int>(std::min_element(angles.begin(), angles.end()) - angles.begin());
        a.small_angle = angles[k] < settings_.min_angle_deg;
        if (!a.small_angle) return a;

        // Rule 1: if the smallest angle sits in an input corner that is itself
        // sharper than the target, no point can fix it. Leave it.
        const int v = t[k];
        if (is_corner(v) && boundary_.corner_angle[v] < settings_.min_angle_deg) {
            a.unfixable = true;
            return a;
        }
        // Rule 2, the sharp-corner rule (Shewchuk calls these "seditious edges"):
        //   1. p and q are the two ends of the edge across from the small angle
        //      (the triangle's shortest edge).
        //   2. p lies on one side of the input polygon, q on another, and the
        //      two sides meet at a corner c sharper than 60 degrees.
        //   3. p and q are the same distance from c (splits near corners use
        //      powers of two, see split_point, so this happens on purpose).
        // Then the triangle is thin only because the corner is sharp. Fixing it
        // would add a point closer to c, making a smaller copy of the same
        // triangle, again and again forever. So leave it, even when the corner
        // is wider than the target.
        const int p = t[(k + 1) % 3], q = t[(k + 2) % 3];
        const int sp = side_of_[p], sq = side_of_[q];
        if (sp >= 0 && sq >= 0 && sp != sq) {
            const Edge& first = boundary_.sides[sp];
            const Edge& second = boundary_.sides[sq];
            for (const int c : {first.first, first.second}) {
                if ((c == second.first || c == second.second) && boundary_.corner_angle[c] < kSharpCornerDeg) {
                    const double dp = distance(points[c], points[p]), dq = distance(points[c], points[q]);
                    if (std::abs(dp - dq) <= kShellTolerance * std::max(dp, dq)) a.unfixable = true;
                }
            }
        }
        return a;
    }

    // Fix one thin (or too long) triangle: add its circumcentre, or split a
    // segment when the circumcentre is not allowed.
    void refine_triangle(const Tri& t, const Assessment& assessment) {
        Cause cause;
        cause.triangle = t;
        cause.reason = assessment.small_angle && !assessment.unfixable ? Cause::Reason::small_angle
                                                                       : Cause::Reason::too_long;
        const auto& points = cdt_.points();
        const Point& a = points[t[0]];
        const Point& b = points[t[1]];
        const Point& c = points[t[2]];
        const Point centre = circumcentre(a, b, c);
        const Point centroid{(a[0] + b[0] + c[0]) / 3.0, (a[1] + b[1] + c[1]) / 3.0};

        // The circumcentre is not added when:
        //   - the path from the triangle's middle to it crosses a segment, so it
        //     is outside the domain or in a hole: split the first segment crossed;
        //   - it would encroach on segments: split all of them.
        // Splitting usually removes the triangle; if not, it goes back on the list.
        std::vector<Edge> to_split;
        if (const auto blocking = blocking_segment(centroid, centre)) {
            to_split.push_back(*blocking);
            cause.reason = Cause::Reason::centre_outside;
        } else {
            to_split = encroached_by(centre);
            if (!to_split.empty()) cause.reason = Cause::Reason::centre_encroaches;
        }
        if (to_split.empty()) {
            auto result = cdt_.insert(centre, t);
            switch (result.status) {
            case Cdt::InsertResult::Status::inserted:
                side_of_.push_back(-1);
                ++stats_.inserted;
                record(std::move(result.step), cause);
                return;
            case Cdt::InsertResult::Status::on_segment:
                // Rare: blocking_segment usually catches a centre on a segment first. If one
                // slips past, split the segment there rather than add a point that makes flat triangles.
                to_split.push_back(result.segment);
                cause.reason = Cause::Reason::centre_on_segment;
                break;
            case Cdt::InsertResult::Status::rejected:  // outside, or on an existing point
                ++stats_.rejected_points;
                return;
            }
        }
        cause.centre = centre;
        bool progressed = false;
        for (const Edge& s : to_split) progressed = split(s, cause) || progressed;
        if (progressed) bad_.push_back(t);
    }

    /// A point that encroaches on segment `s`, if any. As in Ruppert's paper,
    /// only the third points of the one or two triangles on `s` are checked:
    /// in a Delaunay mesh, if a point that can see `s` is inside its circle,
    /// one of these is too.
    std::optional<int> encroacher(const Edge& s) const {
        const auto& points = cdt_.points();
        for (const int v : cdt_.opposite_vertices(s)) {
            if (inside_diametral_circle(points[s.first], points[s.second], points[v])) return v;
        }
        return std::nullopt;
    }

    // Encroachment test: p is inside the circle with diameter ab exactly when
    // the angle a-p-b is over 90 degrees, that is when (a - p) . (b - p) < 0.
    static bool inside_diametral_circle(const Point& a, const Point& b, const Point& p) {
        return (a[0] - p[0]) * (b[0] - p[0]) + (a[1] - p[1]) * (b[1] - p[1]) < 0.0;
    }

    // Every segment that point p would encroach on. Checks all segments.
    std::vector<Edge> encroached_by(const Point& p) const {
        std::vector<Edge> result;
        const auto& points = cdt_.points();
        for (const Edge& s : cdt_.segments()) {
            if (inside_diametral_circle(points[s.first], points[s.second], p)) result.push_back(s);
        }
        return result;
    }

    /// The segment crossed first on the straight path from `from` to `to`, if any.
    /// A quick box test skips far-away segments; the crossing test itself is exact.
    std::optional<Edge> blocking_segment(const Point& from, const Point& to) const {
        const auto& points = cdt_.points();
        std::optional<Edge> nearest;
        double nearest_t = 2.0;
        for (const Edge& s : cdt_.segments()) {
            const Point& a = points[s.first];
            const Point& b = points[s.second];
            if (std::max(a[0], b[0]) < std::min(from[0], to[0]) || std::min(a[0], b[0]) > std::max(from[0], to[0]) ||
                std::max(a[1], b[1]) < std::min(from[1], to[1]) || std::min(a[1], b[1]) > std::max(from[1], to[1])) {
                continue;
            }
            if (!segments_intersect_2d(from.data(), to.data(), a.data(), b.data())) continue;
            const double dx = to[0] - from[0], dy = to[1] - from[1];
            const double ex = b[0] - a[0], ey = b[1] - a[1];
            const double denominator = dx * ey - dy * ex;
            const double t = denominator == 0.0 ? 0.0 : ((a[0] - from[0]) * ey - (a[1] - from[1]) * ex) / denominator;
            if (t < nearest_t) {
                nearest_t = t;
                nearest = s;
            }
        }
        return nearest;
    }

    // Where to split segment s: its midpoint, except next to an input corner.
    Point split_point(const Edge& s) const {
        const auto& points = cdt_.points();
        const bool first_corner = is_corner(s.first), second_corner = is_corner(s.second);
        const Point& a = points[s.first];
        const Point& b = points[s.second];
        if (first_corner == second_corner) return {0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])};
        // Next to a corner (Ruppert's "concentric shells"): put the new point at
        // a power of two from the corner (..., 1/4, 1/2, 1, 2, ...), the one
        // closest to half the length. Then splits on the two sides of a sharp
        // corner land at equal distances from it, instead of creeping closer
        // and closer to the corner one after the other.
        const Point& apex = first_corner ? a : b;
        const Point& far = first_corner ? b : a;
        const double length = distance(apex, far);
        const double f = std::exp2(std::round(std::log2(0.5 * length))) / length;
        return {apex[0] + f * (far[0] - apex[0]), apex[1] + f * (far[1] - apex[1])};
    }

    // Which input side segment s lies on.
    int side_of_segment(const Edge& s) const {
        if (!is_corner(s.first)) return side_of_[s.first];
        if (!is_corner(s.second)) return side_of_[s.second];
        for (std::size_t i = 0; i < boundary_.sides.size(); ++i) {
            if (edge_key(boundary_.sides[i].first, boundary_.sides[i].second) == s) return static_cast<int>(i);
        }
        throw std::logic_error("segment does not belong to any input side");
    }

    // Split segment s. Returns false when the pieces would be too small to be
    // useful (a guard so the run always ends).
    bool split(const Edge& s, Cause cause) {
        if (refused_.count(s) > 0) return false;
        const auto& points = cdt_.points();
        const Point p = split_point(s);
        if (std::min(distance(points[s.first], p), distance(p, points[s.second])) < min_piece_) {
            refused_.insert(s);
            ++stats_.refused_splits;
            return false;
        }
        const int side = side_of_segment(s);
        cause.segment = s;
        cause.shell = is_corner(s.first) != is_corner(s.second);  // see split_point
        Step step = cdt_.split(s, p);
        side_of_.push_back(side);
        ++stats_.splits;
        record(std::move(step), std::move(cause));  // puts the two halves on the check list
        return true;
    }

    // After every change: put the new triangles on the bad list, and their
    // segments on the check list. Keep the step for the steps file if asked.
    void record(Step step, Cause cause) {
        for (const Tri& t : step.added) {
            bad_.push_back(t);
            // The new point may encroach on segments next to it.
            for (int i = 0; i < 3; ++i) {
                if (cdt_.is_segment(t[i], t[(i + 1) % 3])) check_.push_back(edge_key(t[i], t[(i + 1) % 3]));
            }
        }
        if (steps_ != nullptr) {
            step.cause = std::move(cause);
            steps_->push_back(std::move(step));
        }
    }

    Cdt& cdt_;
    const Boundary& boundary_;
    const Settings& settings_;
    std::vector<Step>* steps_;
    std::vector<int> side_of_;  // for every point: the input side it lies on, or -1
    double min_piece_ = 0.0;
    std::deque<Tri> bad_;      // triangles to look at
    std::deque<Edge> check_;   // segments to check for encroachment
    std::set<Edge> refused_;   // segments too short to split again
    RefineStats stats_;
};

}  // namespace

// Step 1: check the polygon and place the boundary points.
Boundary build_boundary(const Domain& domain, double max_edge) {
    std::vector<std::vector<Point>> loops{domain.outer};
    loops.insert(loops.end(), domain.holes.begin(), domain.holes.end());
    for (std::size_t i = 0; i < loops.size(); ++i) {
        if (loops[i].size() < 3) throw std::invalid_argument("every loop needs at least three points");
        const double area = signed_area(loops[i]);
        if (area == 0.0) throw std::invalid_argument("a loop has zero area");
        // Turn the outer loop counter-clockwise and the holes clockwise. Then the
        // domain is on the left of every side, which corner_angle() relies on.
        if ((i == 0) != (area > 0.0)) std::reverse(loops[i].begin(), loops[i].end());
    }
    require_simple(loops);
    // No sides cross, so one point of a hole tells where the whole hole is.
    for (std::size_t i = 1; i < loops.size(); ++i) {
        bool nested = !inside_loop(loops[0], loops[i][0]);
        for (std::size_t j = 1; j < loops.size(); ++j) nested = nested || (j != i && inside_loop(loops[j], loops[i][0]));
        if (nested) throw std::invalid_argument("every hole must lie inside the outer loop and outside the other holes");
    }

    Boundary boundary;
    const auto add_point = [&boundary](const Point& p, int side, double angle) {
        boundary.points.push_back(p);
        boundary.side_of_point.push_back(side);
        boundary.corner_angle.push_back(angle);
        return static_cast<int>(boundary.points.size()) - 1;
    };
    // Walk each loop: add its corners, and with a size limit, equally spaced
    // points along each side. Then join neighbours into segments.
    for (const auto& loop : loops) {
        const int first_side = static_cast<int>(boundary.sides.size());
        const int m = static_cast<int>(loop.size());
        std::vector<int> ring;
        for (int k = 0; k < m; ++k) {
            const Point& a = loop[k];
            const Point& b = loop[(k + 1) % m];
            ring.push_back(add_point(a, -1, corner_angle(loop, static_cast<std::size_t>(k))));
            const int pieces = max_edge > 0.0 ? std::max(1, static_cast<int>(std::ceil(distance(a, b) / max_edge))) : 1;
            for (int j = 1; j < pieces; ++j) {
                const double f = static_cast<double>(j) / pieces;
                ring.push_back(add_point({a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])}, first_side + k, 180.0));
            }
        }
        std::vector<int> corners;
        for (const int v : ring) {
            if (boundary.side_of_point[v] < 0) corners.push_back(v);
        }
        for (int k = 0; k < m; ++k) boundary.sides.push_back({corners[k], corners[(k + 1) % m]});
        for (std::size_t r = 0; r < ring.size(); ++r) {
            boundary.segments.push_back({ring[r], ring[(r + 1) % ring.size()]});
        }
    }
    return boundary;
}

RefineStats refine(Cdt& cdt, const Boundary& boundary, const Settings& settings, std::vector<Step>* steps) {
    return Refiner(cdt, boundary, settings, steps).run();
}

}  // namespace mesher
