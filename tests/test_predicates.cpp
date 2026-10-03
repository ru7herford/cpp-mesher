// Checks the exact predicates, the yes/no tests that steps 2 to 5 rely on.
//
// The exact predicates against exact integer arithmetic (tests/exact.hpp) on
// inputs chosen to be as close to degenerate as binary64 allows: points one
// unit in the last place off a line or a circle, exact ties, and the
// near-collinear grid of Kettner et al. For comparison it also counts how
// often the plain double-precision formula gets the sign wrong on the same inputs.

#include "exact.hpp"
#include "predicates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

using mesher::PredicateSign;
using Point = std::array<double, 2>;

namespace {

int failures = 0;

void check(bool condition, const std::string& what) {
    if (!condition) {
        if (failures < 20) std::printf("FAIL %s\n", what.c_str());
        ++failures;
    }
}

int value(PredicateSign s) { return static_cast<int>(s); }
int sign(double x) { return (x > 0.0) - (x < 0.0); }

int naive_orient(const Point& a, const Point& b, const Point& c) {
    return sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
}

int naive_incircle(const Point& a, const Point& b, const Point& c, const Point& d) {
    const double adx = a[0] - d[0], ady = a[1] - d[1], bdx = b[0] - d[0], bdy = b[1] - d[1];
    const double cdx = c[0] - d[0], cdy = c[1] - d[1];
    const double alift = adx * adx + ady * ady, blift = bdx * bdx + bdy * bdy, clift = cdx * cdx + cdy * cdy;
    return sign(adx * (bdy * clift - cdy * blift) - ady * (bdx * clift - cdx * blift) + alift * (bdx * cdy - cdx * bdy));
}

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

/// x moved by k units in the last place.
double ulps(double x, int k) {
    for (; k > 0; --k) x = std::nextafter(x, INFINITY);
    for (; k < 0; ++k) x = std::nextafter(x, -INFINITY);
    return x;
}

struct Tally {
    long cases = 0, mismatches = 0, naive_wrong = 0, zeros = 0;
};

void orient_case(const Point& a, const Point& b, const Point& c, Tally& t) {
    const int truth = exact::orient(a.data(), b.data(), c.data());
    const int got = value(mesher::orient2d(a.data(), b.data(), c.data()));
    ++t.cases;
    t.zeros += truth == 0;
    t.naive_wrong += naive_orient(a, b, c) != truth;
    if (got != truth) {
        ++t.mismatches;
        check(false, "orient2d disagrees with exact arithmetic");
    }
    // Swapping two points flips the sign; rotating them keeps it.
    check(value(mesher::orient2d(b.data(), a.data(), c.data())) == -got, "orient2d is not antisymmetric");
    check(value(mesher::orient2d(b.data(), c.data(), a.data())) == got, "orient2d changes under rotation");
}

void incircle_case(const Point& a, const Point& b, const Point& c, const Point& d, Tally& t) {
    const int truth = exact::incircle(a.data(), b.data(), c.data(), d.data());
    const int got = value(mesher::incircle(a.data(), b.data(), c.data(), d.data()));
    ++t.cases;
    t.zeros += truth == 0;
    t.naive_wrong += naive_incircle(a, b, c, d) != truth;
    if (got != truth) {
        ++t.mismatches;
        check(false, "incircle disagrees with exact arithmetic");
    }
    check(value(mesher::incircle(b.data(), a.data(), c.data(), d.data())) == -got, "incircle is not antisymmetric");
    check(value(mesher::incircle(b.data(), c.data(), a.data(), d.data())) == got, "incircle changes under rotation");
}

void print(const char* name, const Tally& t) {
    std::printf("%-44s %6ld cases, %ld exact ties, %ld mismatches; plain doubles get %ld signs wrong\n", name, t.cases,
                t.zeros, t.mismatches, t.naive_wrong);
}

}  // namespace

int main() {
    Random random{20261002};

    // Kettner et al.'s grid: p = (0.5 + i u, 0.5 + j u), u = 2^-53, against the line through (12, 12) and (24, 24).
    Tally grid;
    const double u = std::ldexp(1.0, -53);
    const Point q{12, 12}, r{24, 24};
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 256; ++j) orient_case({0.5 + i * u, 0.5 + j * u}, q, r, grid);
    }
    print("orient2d, near-collinear grid (256 x 256)", grid);

    // Random near-collinear triples: a point on the line through two others, rounded, then nudged by a few ulps.
    Tally collinear;
    for (int n = 0; n < 20000; ++n) {
        const double scale = std::ldexp(1.0, random.integer(-20, 20));
        const Point a{random.uniform(-1, 1) * scale, random.uniform(-1, 1) * scale};
        const Point b{random.uniform(-1, 1) * scale, random.uniform(-1, 1) * scale};
        const double t = random.uniform(-2, 3);
        const Point c{ulps(a[0] + t * (b[0] - a[0]), random.integer(-2, 2)), ulps(a[1] + t * (b[1] - a[1]), random.integer(-2, 2))};
        orient_case(a, b, c, collinear);
    }
    print("orient2d, random near-collinear triples", collinear);

    // Random near-cocircular quadruples: four points on one circle, rounded, the fourth nudged by a few ulps.
    Tally cocircular;
    for (int n = 0; n < 20000; ++n) {
        const double scale = std::ldexp(1.0, random.integer(-20, 20));
        const double offset = random.integer(0, 1) ? random.uniform(-1e3, 1e3) * scale : 0.0;
        const double radius = random.uniform(0.1, 1.0) * scale;
        std::array<double, 4> angle{};
        for (double& theta : angle) theta = random.uniform(0.0, 2.0 * std::acos(-1.0));
        std::sort(angle.begin(), angle.begin() + 3);  // a, b, c counter-clockwise
        std::array<Point, 4> p{};
        for (int k = 0; k < 4; ++k) p[k] = {offset + radius * std::cos(angle[k]), offset + radius * std::sin(angle[k])};
        p[3] = {ulps(p[3][0], random.integer(-3, 3)), ulps(p[3][1], random.integer(-3, 3))};
        if (exact::orient(p[0].data(), p[1].data(), p[2].data()) <= 0) continue;  // three points too close together
        incircle_case(p[0], p[1], p[2], p[3], cocircular);
    }
    print("incircle, random near-cocircular quadruples", cocircular);

    // Exact ties: integer points on the circle x^2 + y^2 = 25 and on lines, scaled and shifted by powers of two
    // (which keeps them exact). Every answer must be exactly zero.
    Tally ties;
    const std::array<Point, 8> circle{{{5, 0}, {4, 3}, {3, 4}, {0, 5}, {-3, 4}, {-5, 0}, {-4, -3}, {3, -4}}};
    for (const int power : {-60, -20, 0, 20, 60}) {
        for (const double shift : {0.0, 1024.0, -3.0e6}) {
            const auto at = [&](const Point& p) { return Point{std::ldexp(p[0], power) + std::ldexp(shift, power), std::ldexp(p[1], power) + std::ldexp(shift, power)}; };
            for (int i = 0; i + 3 < 8; ++i) incircle_case(at(circle[i]), at(circle[i + 1]), at(circle[i + 2]), at(circle[i + 3]), ties);
            for (int k = 1; k < 6; ++k) orient_case(at({0, 0}), at({1.0 * k, 2.0 * k}), at({-3.0 * k, -6.0 * k}), ties);
        }
    }
    check(ties.zeros == ties.cases, "an exact tie did not come out as zero");
    print("exact ties (cocircular and collinear)", ties);

    // Segment intersection on touching, overlapping and near-miss cases.
    struct SegmentCase {
        Point a, b, c, d;
        bool meet;
    };
    const double near = std::nextafter(1.0, 2.0);
    const SegmentCase segment_cases[] = {
        {{0, 0}, {2, 2}, {0, 2}, {2, 0}, true},         // proper crossing
        {{0, 0}, {1, 1}, {1, 1}, {2, 0}, true},         // shared end point
        {{0, 0}, {2, 0}, {1, 0}, {1, 1}, true},         // T: an end point on the other segment
        {{0, 0}, {2, 0}, {1, 0}, {3, 0}, true},         // collinear overlap
        {{0, 0}, {1, 0}, {2, 0}, {3, 0}, false},        // collinear, apart
        {{0, 0}, {1, 0}, {0, 1}, {1, 1}, false},        // parallel
        {{0, 0}, {2, 0}, {1, near - 1.0}, {1, 1}, false},  // stops 2^-52 above the first segment
    };
    int segment_checks = 0;
    for (const SegmentCase& s : segment_cases) {
        const bool meet = mesher::segments_intersect_2d(s.a.data(), s.b.data(), s.c.data(), s.d.data());
        const bool swapped = mesher::segments_intersect_2d(s.c.data(), s.d.data(), s.b.data(), s.a.data());
        check(meet == s.meet && swapped == s.meet, "segments_intersect_2d case " + std::to_string(segment_checks));
        ++segment_checks;
    }
    // Random near-touching pairs against a reference built on exact orientations.
    const auto exact_meet = [](const Point& a, const Point& b, const Point& c, const Point& d) {
        const auto o = [](const Point& p, const Point& q, const Point& r) { return exact::orient(p.data(), q.data(), r.data()); };
        const auto within = [](const Point& p, const Point& q, const Point& r) {
            return std::min(p[0], q[0]) <= r[0] && r[0] <= std::max(p[0], q[0]) && std::min(p[1], q[1]) <= r[1] && r[1] <= std::max(p[1], q[1]);
        };
        const int abc = o(a, b, c), abd = o(a, b, d), cda = o(c, d, a), cdb = o(c, d, b);
        if (abc && abd && cda && cdb) return abc != abd && cda != cdb;
        return (!abc && within(a, b, c)) || (!abd && within(a, b, d)) || (!cda && within(c, d, a)) || (!cdb && within(c, d, b));
    };
    for (int n = 0; n < 2000; ++n) {
        const Point a{random.uniform(-1, 1), random.uniform(-1, 1)}, b{random.uniform(-1, 1), random.uniform(-1, 1)};
        const double t = random.uniform(-0.2, 1.2);  // c lands near the line through a and b, sometimes past its ends
        const Point c{ulps(a[0] + t * (b[0] - a[0]), random.integer(-1, 1)), ulps(a[1] + t * (b[1] - a[1]), random.integer(-1, 1))};
        const Point d{random.uniform(-1, 1), random.uniform(-1, 1)};
        const bool meet = mesher::segments_intersect_2d(a.data(), b.data(), c.data(), d.data());
        check(meet == exact_meet(a, b, c, d), "segments_intersect_2d disagrees with exact arithmetic");
        ++segment_checks;
    }

    // Inputs outside the range where the arithmetic is exact are refused, not guessed.
    int refused = 0;
    for (const Point bad : {Point{NAN, 0}, Point{INFINITY, 0}, Point{1e300, 0}, Point{1e-300, 0}}) {
        const Point a{0, 0}, b{1, 0};
        try {
            mesher::orient2d(a.data(), b.data(), bad.data());
        } catch (const std::runtime_error&) {
            ++refused;
        }
    }
    check(refused == 4, "an out-of-range coordinate was not refused");

    const long cases = grid.cases + collinear.cases + cocircular.cases + ties.cases;
    const long mismatches = grid.mismatches + collinear.mismatches + cocircular.mismatches + ties.mismatches;
    std::printf("segment cases %d, refused out-of-range inputs %d\n", segment_checks, refused);
    std::printf("SUMMARY {\"test\": \"predicates\", \"orient_grid\": %ld, \"orient_grid_naive_wrong\": %ld, "
                "\"orient_random\": %ld, \"orient_random_naive_wrong\": %ld, \"incircle_random\": %ld, "
                "\"incircle_random_naive_wrong\": %ld, \"exact_ties\": %ld, \"segment_cases\": %d, "
                "\"cases\": %ld, \"mismatches\": %ld, \"failures\": %d}\n",
                grid.cases, grid.naive_wrong, collinear.cases, collinear.naive_wrong, cocircular.cases,
                cocircular.naive_wrong, ties.cases, segment_checks, cases, mismatches, failures);
    std::printf(failures == 0 ? "all checks passed\n" : "%d checks failed\n", failures);
    return failures == 0 ? 0 : 1;
}
