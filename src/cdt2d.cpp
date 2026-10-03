// Steps 2, 3 and 4: build the first mesh. Also the two ways step 5 changes
// it: insert() (a point inside) and split() (a point on a segment).
//
// Every yes/no geometric question here goes through the exact predicates
// (predicates.cpp): orient2d (which side of a line is a point on?) and
// incircle (is a point inside a triangle's circumcircle?). They never give a
// wrong answer because of rounding, so the mesh never gets tangled.

#include "cdt2d.hpp"

#include "predicates.hpp"

#include <algorithm>
#include <deque>
#include <stdexcept>

namespace mesher {
namespace {

// How far the big starting triangle reaches past the domain, in domain sizes.
constexpr double kEnclosingMargin = 16.0;

// One number for the edge a->b, used as a lookup key. a->b and b->a differ:
// they belong to the triangles on the two sides of the edge.
std::uint64_t half_edge_key(int a, int b) {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32) |
           static_cast<std::uint32_t>(b);
}

}  // namespace

Edge edge_key(int a, int b) { return a < b ? Edge{a, b} : Edge{b, a}; }

std::pair<Point, Point> bounding_box(const std::vector<Point>& points) {
    Point low = points.at(0), high = points.at(0);
    for (const Point& p : points) {
        low = {std::min(low[0], p[0]), std::min(low[1], p[1])};
        high = {std::max(high[0], p[0]), std::max(high[1], p[1])};
    }
    return {low, high};
}

// Steps 2, 3 and 4.
Cdt::Cdt(std::vector<Point> points, const std::vector<Edge>& segments, BuildTrace* trace)
    : points_(std::move(points)) {
    const int n = static_cast<int>(points_.size());
    if (n < 3 || segments.size() < 3) {
        throw std::invalid_argument("a domain needs at least three points and three segments");
    }
    std::vector<Point> sorted = points_;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
        throw std::invalid_argument("the input contains duplicate points");
    }

    // Step 2 needs a mesh to put points into. Start with one big triangle
    // around everything, so every point falls inside some triangle. Its three
    // corners are removed with the rest of the outside in step 4.
    const auto [low, high] = bounding_box(points_);
    const double margin = kEnclosingMargin * std::max(high[0] - low[0], high[1] - low[1]);
    const double cx = 0.5 * (low[0] + high[0]), cy = 0.5 * (low[1] + high[1]);
    points_.push_back({cx - 2.0 * margin, cy - margin});
    points_.push_back({cx + 2.0 * margin, cy - margin});
    points_.push_back({cx, cy + 2.0 * margin});
    int seed = add({n, n + 1, n + 2});
    if (trace) std::copy(points_.begin() + n, points_.end(), trace->enclosing.begin());

    // Step 2: add the points one at a time (Bowyer-Watson).
    for (int i = 0; i < n; ++i) {
        recording_ = trace ? &trace->insertions.emplace_back() : nullptr;
        seed = fill_cavity(cavity(i, locate(points_[i], seed)), i);
    }
    for (const auto& [a, b] : segments) {
        if (a < 0 || b < 0 || a >= n || b >= n || a == b) {
            throw std::invalid_argument("segment endpoint is invalid");
        }
        if (!segments_.insert(edge_key(a, b)).second) {
            throw std::invalid_argument("the input contains duplicate segments");
        }
    }
    // Step 3: force in every segment that step 2 did not produce.
    for (const Edge& segment : segments_) {
        Delta flips;
        recording_ = trace ? &flips : nullptr;
        recover_segment(segment);
        if (trace && !flips.added.empty()) trace->recoveries.push_back({segment, std::move(flips)});
    }
    // Step 4: remove the outside, and the big triangle's corners with it.
    recording_ = trace ? &trace->exterior : nullptr;
    remove_exterior(n);
    points_.resize(static_cast<std::size_t>(n));

    // The flips of step 3 ignore the empty-circle rule, so some edges may now
    // break it. Check every edge and flip the ones that do (Lawson flips).
    // Step 5 relies on this: Bowyer-Watson only works on a Delaunay mesh.
    std::vector<Edge> all_edges;
    for (const Tri& t : triangles()) {
        for (int i = 0; i < 3; ++i) all_edges.push_back({t[i], t[(i + 1) % 3]});
    }
    recording_ = trace ? &trace->legalize : nullptr;
    legalize(std::move(all_edges));
    recording_ = nullptr;
}

// Step 5, adding a circumcentre: Bowyer-Watson for one point. Refuses points
// outside the domain, and reports a point that lands exactly on a segment.
Cdt::InsertResult Cdt::insert(const Point& p, const Tri& near) {
    InsertResult result{InsertResult::Status::rejected, {}, {}};
    const int located = locate(p, triangle_at(near[0], near[1]));
    if (located < 0) return result;
    for (int v : tris_[static_cast<std::size_t>(located)]) {
        if (points_[static_cast<std::size_t>(v)] == p) return result;
    }

    const int index = static_cast<int>(points_.size());
    points_.push_back(p);
    const Cavity region = cavity(index, located);
    // Every edge around the hole must have p strictly on its left, or the new
    // triangles (edge + p) would be flat or flipped over. Check before changing anything.
    for (const auto& [a, b] : region.boundary) {
        if (orient2d(points_[a].data(), points_[b].data(), p.data()) == PredicateSign::positive) {
            continue;
        }
        points_.pop_back();
        // p is exactly on a segment. It cannot be joined to points across the
        // segment, so the caller splits the segment there instead.
        if (is_segment(a, b)) {
            result.status = InsertResult::Status::on_segment;
            result.segment = edge_key(a, b);
            return result;
        }
        throw std::logic_error("Bowyer-Watson cavity is not star-shaped");
    }

    result.status = InsertResult::Status::inserted;
    result.step.kind = Step::Kind::insert;
    result.step.point = index;
    recording_ = &result.step;
    fill_cavity(region, index);
    recording_ = nullptr;
    return result;
}

// Step 5, splitting a segment: put p on segment ab. Each triangle a-b-w on
// either side becomes two, a-p-w and p-b-w, and ab becomes the two segments ap and pb.
Step Cdt::split(const Edge& segment, const Point& p) {
    const auto [a, b] = segment;
    if (!is_segment(a, b)) throw std::logic_error("split target is not a segment");
    const int index = static_cast<int>(points_.size());
    points_.push_back(p);

    // Check that every new triangle will be counter-clockwise (not flat or
    // flipped) before changing anything.
    for (const auto& [u, v] : {Edge{a, b}, Edge{b, a}}) {
        const int slot = triangle_at(u, v);
        if (slot < 0) continue;
        const double* w = points_[static_cast<std::size_t>(apex(slot, u))].data();
        if (orient2d(points_[u].data(), p.data(), w) != PredicateSign::positive ||
            orient2d(p.data(), points_[v].data(), w) != PredicateSign::positive) {
            points_.pop_back();
            throw std::invalid_argument("split point is not strictly inside its segment");
        }
    }

    Step step;
    step.kind = Step::Kind::split;
    step.point = index;
    recording_ = &step;
    std::vector<Edge> link;
    for (const auto& [u, v] : {Edge{a, b}, Edge{b, a}}) {
        const int slot = triangle_at(u, v);
        if (slot < 0) continue;
        const int w = apex(slot, u);
        remove(slot);
        add({u, index, w});
        add({index, v, w});
        link.push_back({w, u});
        link.push_back({v, w});
    }
    segments_.erase(segment);
    segments_.insert(edge_key(a, index));
    segments_.insert(edge_key(index, b));
    // The new triangles may break the empty-circle rule. Flip the edges around
    // them until it holds again, so the next Bowyer-Watson step works.
    legalize(std::move(link));
    recording_ = nullptr;
    return step;
}

std::vector<Tri> Cdt::triangles() const {
    std::vector<Tri> result;
    for (std::size_t slot = 0; slot < tris_.size(); ++slot) {
        if (alive_[slot]) result.push_back(tris_[slot]);
    }
    return result;
}

bool Cdt::contains(const Tri& t) const {
    const int slot = triangle_at(t[0], t[1]);
    return slot >= 0 && apex(slot, t[0]) == t[2];
}

std::vector<int> Cdt::opposite_vertices(const Edge& segment) const {
    std::vector<int> result;
    for (const auto& [u, v] : {segment, Edge{segment.second, segment.first}}) {
        const int slot = triangle_at(u, v);
        if (slot >= 0) result.push_back(apex(slot, u));
    }
    return result;
}

// Put a triangle into the mesh and return its slot. Refuses a triangle that is
// not strictly counter-clockwise: that would be a bug, so it stops loudly.
int Cdt::add(const Tri& t) {
    if (orient2d(points_[t[0]].data(), points_[t[1]].data(), points_[t[2]].data()) !=
        PredicateSign::positive) {
        throw std::logic_error("refusing to create an inverted or degenerate triangle");
    }
    int slot;
    if (!free_slots_.empty()) {
        slot = free_slots_.back();
        free_slots_.pop_back();
        tris_[static_cast<std::size_t>(slot)] = t;
        alive_[static_cast<std::size_t>(slot)] = 1;
    } else {
        slot = static_cast<int>(tris_.size());
        tris_.push_back(t);
        alive_.push_back(1);
    }
    for (int i = 0; i < 3; ++i) {
        if (!half_edges_.emplace(half_edge_key(t[i], t[(i + 1) % 3]), slot).second) {
            throw std::logic_error("an edge would be used twice in the same direction");
        }
    }
    if (recording_ != nullptr) recording_->added.push_back(t);
    return slot;
}

// Take a triangle out of the mesh and free its slot for reuse.
void Cdt::remove(int slot) {
    const Tri t = tris_[static_cast<std::size_t>(slot)];
    for (int i = 0; i < 3; ++i) half_edges_.erase(half_edge_key(t[i], t[(i + 1) % 3]));
    alive_[static_cast<std::size_t>(slot)] = 0;
    free_slots_.push_back(slot);
    if (recording_ != nullptr) {
        // A triangle made and removed inside the same change is left out of the record.
        // (A plain search is fine: one change touches only a few triangles.)
        auto& added = recording_->added;
        const auto found = std::find(added.begin(), added.end(), t);
        if (found != added.end()) {
            added.erase(found);
        } else {
            recording_->removed.push_back(t);
        }
    }
}

int Cdt::triangle_at(int a, int b) const {
    const auto found = half_edges_.find(half_edge_key(a, b));
    return found == half_edges_.end() ? -1 : found->second;
}

int Cdt::apex(int slot, int a) const {
    const Tri& t = tris_[static_cast<std::size_t>(slot)];
    const int i = t[0] == a ? 0 : t[1] == a ? 1 : 2;
    return t[(i + 2) % 3];
}

// Whether p is inside the triangle or on its edge: p is not on the right of any edge.
bool Cdt::triangle_contains(int slot, const Point& p) const {
    const Tri& t = tris_[static_cast<std::size_t>(slot)];
    for (int i = 0; i < 3; ++i) {
        if (orient2d(points_[t[i]].data(), points_[t[(i + 1) % 3]].data(), p.data()) ==
            PredicateSign::negative) {
            return false;
        }
    }
    return true;
}

// Find a triangle that contains p, starting from triangle `seed`.
int Cdt::locate(const Point& p, int seed) const {
    // Walk: if p is on the far side of one of this triangle's edges, step into
    // the neighbour across that edge, and repeat until p is inside.
    int current = seed;
    for (std::size_t step = 0; current >= 0 && step < tris_.size(); ++step) {
        const Tri& t = tris_[static_cast<std::size_t>(current)];
        int next = -1;
        bool inside = true;
        for (std::size_t k = 0; k < 3; ++k) {
            // Start with a different edge each time, so the walk cannot go round in circles.
            const int i = static_cast<int>((k + step) % 3);
            const int a = t[i], b = t[(i + 1) % 3];
            if (orient2d(points_[a].data(), points_[b].data(), p.data()) == PredicateSign::negative) {
                inside = false;
                next = triangle_at(b, a);
                if (next >= 0) break;
            }
        }
        if (inside) return current;
        current = next;
    }
    // The walk can get stuck at a hole or a dent in the boundary, because
    // there is no triangle to step into. Then check every triangle.
    for (std::size_t slot = 0; slot < tris_.size(); ++slot) {
        if (alive_[slot] && triangle_contains(static_cast<int>(slot), p)) {
            return static_cast<int>(slot);
        }
    }
    return -1;
}

// Bowyer-Watson, part 1: find the triangles a new point removes. These are
// the triangles whose circumcircle has the point strictly inside. They touch
// each other, so start at the triangle that holds the point and spread to
// neighbours, never across a segment.
Cdt::Cavity Cdt::cavity(int point, int located) const {
    const double* p = points_[static_cast<std::size_t>(point)].data();
    Cavity result;
    std::set<int> tested;
    std::vector<int> stack{located};
    while (!stack.empty()) {
        const int slot = stack.back();
        stack.pop_back();
        if (!tested.insert(slot).second) continue;
        const Tri& t = tris_[static_cast<std::size_t>(slot)];
        const auto sign = incircle(points_[t[0]].data(), points_[t[1]].data(), points_[t[2]].data(), p);
        // A circumcircle that passes exactly through the point does not count.
        // That way every edge around the hole has the point strictly on its
        // left, and no new triangle can come out flat.
        if (sign != PredicateSign::positive) continue;
        result.triangles.push_back(slot);
        for (int i = 0; i < 3; ++i) {
            const int a = t[i], b = t[(i + 1) % 3];
            if (is_segment(a, b)) continue;  // the cavity never grows across a segment
            const int neighbour = triangle_at(b, a);
            if (neighbour >= 0) stack.push_back(neighbour);
        }
    }
    if (result.triangles.empty()) throw std::logic_error("point location found no conflict cavity");

    // The edges around the hole: edges of removed triangles whose neighbour stays
    // (or is missing, or lies across a segment).
    const std::set<int> members(result.triangles.begin(), result.triangles.end());
    for (const int slot : result.triangles) {
        const Tri& t = tris_[static_cast<std::size_t>(slot)];
        for (int i = 0; i < 3; ++i) {
            const int a = t[i], b = t[(i + 1) % 3];
            const int neighbour = triangle_at(b, a);
            if (neighbour < 0 || members.count(neighbour) == 0 || is_segment(a, b)) {
                result.boundary.push_back({a, b});
            }
        }
    }
    return result;
}

// Bowyer-Watson, part 2: remove the triangles and join every edge around the
// hole to the new point. Returns one of the new triangles.
int Cdt::fill_cavity(const Cavity& region, int point) {
    for (const int slot : region.triangles) remove(slot);
    int last = -1;
    for (const auto& [a, b] : region.boundary) last = add({a, b, point});
    return last;
}

// Edge flip: the two triangles on either side of edge ab form a four-sided
// shape a-d-b-c. Replace edge ab with the other diagonal, cd. Returns cd.
Edge Cdt::flip(int a, int b) {
    // Triangles (a,b,c) and (b,a,d) become (a,d,c) and (d,b,c).
    const int left = triangle_at(a, b), right = triangle_at(b, a);
    const int c = apex(left, a), d = apex(right, b);
    remove(left);
    remove(right);
    add({a, d, c});
    add({d, b, c});
    return edge_key(c, d);
}

// Whether edges ab and cd cross at a point inside both (not at an end, not
// just touching): c and d are on opposite sides of line ab, and a and b on
// opposite sides of line cd.
bool Cdt::proper_crossing(int a, int b, int c, int d) const {
    if (a == c || a == d || b == c || b == d) return false;
    const auto o = [this](int i, int j, int k) {
        return orient2d(points_[i].data(), points_[j].data(), points_[k].data());
    };
    const auto opposite = [](PredicateSign x, PredicateSign y) {
        return (x == PredicateSign::positive && y == PredicateSign::negative) ||
               (x == PredicateSign::negative && y == PredicateSign::positive);
    };
    return opposite(o(a, b, c), o(a, b, d)) && opposite(o(c, d, a), o(c, d, b));
}

// Step 3 for one segment ab (Sloan's method): if ab is not an edge of the
// mesh, other edges cross it. Flip those crossing edges away until ab appears.
void Cdt::recover_segment(const Edge& segment) {
    const auto [a, b] = segment;
    if (triangle_at(a, b) >= 0 || triangle_at(b, a) >= 0) return;
    std::deque<Edge> crossing;
    for (const Tri& t : triangles()) {
        for (int i = 0; i < 3; ++i) {
            const int u = t[i], v = t[(i + 1) % 3];
            if (u < v && proper_crossing(a, b, u, v)) crossing.push_back({u, v});
        }
    }
    if (crossing.empty()) throw std::invalid_argument("a segment passes through an input point");

    // Step by step:
    //   1. Take the next crossing edge uv. Its two triangles form a four-sided
    //      shape u-c-v-d.
    //   2. If that shape has a dent (it is not convex), flipping would make a
    //      flipped-over triangle. Put uv at the back of the queue and go on.
    //   3. Otherwise flip uv to cd. If cd still crosses ab, queue it too.
    //   4. Stop when nothing crosses ab: then ab is an edge.
    // Sloan showed that some crossing edge can always be flipped, so this ends.
    // "stalled" counts edges skipped in a row, as a guard against a bug.
    std::size_t stalled = 0;
    while (!crossing.empty()) {
        const auto [u, v] = crossing.front();
        crossing.pop_front();
        const int c = apex(triangle_at(u, v), u), d = apex(triangle_at(v, u), v);
        if (!proper_crossing(c, d, u, v)) {
            crossing.push_back({u, v});
            if (++stalled > crossing.size()) throw std::runtime_error("segment recovery stalled");
            continue;
        }
        stalled = 0;
        const Edge created = flip(u, v);
        if (proper_crossing(a, b, created.first, created.second)) crossing.push_back(created);
    }
}

// Step 4: remove every triangle outside the domain (outside the outer loop,
// or inside a hole).
void Cdt::remove_exterior(int real_point_count) {
    // Step by step:
    //   1. Triangles that touch a corner of the big starting triangle are
    //      surely outside. Mark them 0 (outside).
    //   2. Spread to neighbours. Crossing an ordinary edge keeps the mark;
    //      crossing a segment switches it (0 to 1 or 1 to 0). So inside the
    //      polygon the mark is 1, inside a hole it is 0 again, and so on.
    //   3. If a triangle is reached with two different marks, the segments do
    //      not form closed loops: stop with an error.
    //   4. Remove every triangle not marked 1.
    std::vector<int> parity(tris_.size(), -1);
    std::deque<int> queue;
    for (std::size_t slot = 0; slot < tris_.size(); ++slot) {
        const Tri& t = tris_[slot];
        if (alive_[slot] && std::max({t[0], t[1], t[2]}) >= real_point_count) {
            parity[slot] = 0;
            queue.push_back(static_cast<int>(slot));
        }
    }
    while (!queue.empty()) {
        const int slot = queue.front();
        queue.pop_front();
        const Tri& t = tris_[static_cast<std::size_t>(slot)];
        for (int i = 0; i < 3; ++i) {
            const int a = t[i], b = t[(i + 1) % 3];
            const int neighbour = triangle_at(b, a);
            if (neighbour < 0) continue;
            const int next = parity[static_cast<std::size_t>(slot)] ^ (is_segment(a, b) ? 1 : 0);
            int& assigned = parity[static_cast<std::size_t>(neighbour)];
            if (assigned < 0) {
                assigned = next;
                queue.push_back(neighbour);
            } else if (assigned != next) {
                throw std::invalid_argument("segments do not form closed loops");
            }
        }
    }
    for (std::size_t slot = 0; slot < tris_.size(); ++slot) {
        if (alive_[slot] && parity[slot] != 1) remove(static_cast<int>(slot));
    }
}

// Lawson flips: make every edge on the stack pass the empty-circle rule.
// Edge ab with triangles a-b-c and b-a-d fails when d is inside the
// circumcircle of a-b-c. Then flip ab to cd, and check the four outer edges
// of the new pair again, since they now have new neighbours. Segments are
// never flipped: they must stay.
void Cdt::legalize(std::vector<Edge> stack) {
    while (!stack.empty()) {
        const auto [a, b] = stack.back();
        stack.pop_back();
        if (is_segment(a, b)) continue;
        const int left = triangle_at(a, b), right = triangle_at(b, a);
        if (left < 0 || right < 0) continue;
        const int c = apex(left, a), d = apex(right, b);
        if (incircle(points_[a].data(), points_[b].data(), points_[c].data(), points_[d].data()) !=
            PredicateSign::positive) {
            continue;
        }
        // When d is inside that circle, the shape a-d-b-c never has a dent, so the flip is always safe.
        flip(a, b);
        stack.insert(stack.end(), {{a, d}, {d, b}, {b, c}, {c, a}});
    }
}

}  // namespace mesher
