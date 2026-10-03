// The exact predicates (Shewchuk's idea, with a simpler error bound).
// Used by steps 2 to 5: every "which side?" and "inside the circle?" question.
//
// Each test is the sign of a small determinant. How the sign is made exact:
//   1. Compute the determinant with ordinary doubles. Fast, but rounded.
//   2. Also compute a bound on how much rounding could have changed it.
//   3. If the value is further from zero than the bound, its sign is right.
//      Return it. This is almost always the case.
//   4. Otherwise compute the determinant again with no rounding at all, as an
//      "expansion": a list of doubles whose exact sum is the true value.
//      Two tricks make that possible:
//        two-sum (Knuth):   a + b = s + e exactly, s the rounded sum, e a double;
//        product (Dekker):  a * b = m + e exactly, m the rounded product.
//      Keep every e instead of throwing it away, and nothing is lost.
//   5. The largest non-zero part of the expansion has the sign of the total.
// This needs every + and * rounded on its own, as plain IEEE doubles. A fused
// multiply-add rounds a * b + c only once, which breaks these tricks, so the
// build turns it off (-ffp-contract=off in CMakeLists.txt).

#include "predicates.hpp"

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace mesher {
namespace {

using Expansion = std::vector<double>;

constexpr double kSplitter = 134217729.0;  // 2^27 + 1, used to cut a double into two halves
// Coordinates must lie in this range, so that no exact product overflows to
// infinity and no error term becomes too small to store exactly.
constexpr double kMaximumAdmitted = 4.017345110647476e59;   // 2^198.
constexpr double kMinimumAdmitted = 3.054936363499605e-151;  // 2^-500.
// The rounding bound is gamma_n * P, with gamma_n = n eps / (1 - n eps) and
// eps = 2^-52 (see sign_is_certain). n is a generous count of rounded operations.
// Shewchuk's careful analysis gives about 3 eps (orient2d) and 10 eps
// (incircle); a looser bound only sends a few more cases down the exact route.
constexpr double kOrientOperations = 32.0;
constexpr double kIncircleOperations = 256.0;

/// Refuse to run if doubles are not rounded to nearest: the exact tricks below need it.
void ensure_environment() {
    static_assert(std::numeric_limits<double>::is_iec559);
    static_assert(std::numeric_limits<double>::radix == 2);
    static_assert(std::numeric_limits<double>::digits == 53);
    if (std::fegetround() != FE_TONEAREST) {
        throw std::runtime_error("exact predicates require round-to-nearest binary64");
    }
}

/// Refuse coordinates outside the range where the exact arithmetic is known to work.
void require_admitted(const double* point) {
    for (int index = 0; index < 2; ++index) {
        const double value = point[index];
        const double magnitude = std::abs(value);
        if (!std::isfinite(value)) {
            throw std::runtime_error("exact predicates require finite coordinates");
        }
        if (magnitude != 0.0 && (magnitude < kMinimumAdmitted || magnitude > kMaximumAdmitted)) {
            throw std::runtime_error("coordinate outside the exact predicate range");
        }
    }
}

// Sign of one double.
PredicateSign sign_of(double value) {
    if (value < 0.0) return PredicateSign::negative;
    if (value > 0.0) return PredicateSign::positive;
    return PredicateSign::zero;
}

/// Sign of an expansion. Its parts are stored smallest first and do not
/// overlap, so the last non-zero part is larger than all the others together
/// and decides the sign.
PredicateSign sign_of(const Expansion& expansion) {
    for (auto it = expansion.rbegin(); it != expansion.rend(); ++it) {
        if (*it < 0.0) return PredicateSign::negative;
        if (*it > 0.0) return PredicateSign::positive;
    }
    return PredicateSign::zero;
}

/// Two-sum (Knuth): sum = rounded first + second, and error = what rounding
/// lost, so first + second = sum + error exactly.
void two_sum(double first, double second, double& sum, double& error) {
    sum = first + second;
    const double virtual_second = sum - first;
    const double virtual_first = sum - virtual_second;
    error = (second - virtual_second) + (first - virtual_first);
}

/// The same for a difference: first - second = difference + error exactly.
void two_diff(double first, double second, double& difference, double& error) {
    difference = first - second;
    if (!std::isfinite(difference)) {
        throw std::runtime_error("coordinate difference outside the exact predicate range");
    }
    const double virtual_second = first - difference;
    const double virtual_first = difference + virtual_second;
    error = (first - virtual_first) + (virtual_second - second);
}

/// Cut a double into a high and a low half, each with at most 26 significant
/// bits, so that multiplying two halves never needs rounding (Dekker).
void split(double value, double& high, double& low) {
    const double product = kSplitter * value;
    if (!std::isfinite(product)) {
        throw std::runtime_error("coordinate outside the exact predicate range");
    }
    const double difference = product - value;
    high = product - difference;
    low = value - high;
}

/// Exact product (Dekker): product = rounded first * second, and error = what
/// rounding lost. Multiply the halves pairwise (each exact) and subtract them
/// from the rounded product to find the error.
void two_product(double first, double second, double& product, double& error) {
    product = first * second;
    if (!std::isfinite(product)) {
        throw std::runtime_error("product outside the exact predicate range");
    }
    if (first != 0.0 && second != 0.0 &&
        (product == 0.0 || std::fpclassify(product) == FP_SUBNORMAL)) {
        throw std::runtime_error("product underflows the exact predicate range");
    }
    double first_high, first_low, second_high, second_low;
    split(first, first_high, first_low);
    split(second, second_high, second_low);
    const double first_error = product - first_high * second_high;
    const double second_error = first_error - first_low * second_high;
    const double third_error = second_error - first_high * second_low;
    error = first_low * second_low - third_error;
    if (error != 0.0 && std::fpclassify(error) == FP_SUBNORMAL) {
        throw std::runtime_error("product residual underflows the exact predicate range");
    }
}

/// Add one double to an expansion: two-sum it with each part in turn, keep
/// every error, carry the sum. Zero parts are dropped.
Expansion grow_expansion(const Expansion& expansion, double scalar) {
    if (scalar == 0.0) return expansion;
    Expansion result;
    result.reserve(expansion.size() + 1);
    double accumulator = scalar;
    for (double component : expansion) {
        double next, error;
        two_sum(accumulator, component, next, error);
        if (error != 0.0) result.push_back(error);
        accumulator = next;
    }
    if (accumulator != 0.0) result.push_back(accumulator);
    return result;
}

// Sum of two expansions: add the parts of the second one by one.
Expansion expansion_sum(const Expansion& first, const Expansion& second) {
    Expansion result = first;
    for (double component : second) result = grow_expansion(result, component);
    return result;
}

// first - second: negate every part of second, then add.
Expansion expansion_difference(const Expansion& first, const Expansion& second) {
    Expansion negated;
    negated.reserve(second.size());
    for (double component : second) negated.push_back(-component);
    return expansion_sum(first, negated);
}

// Expansion times one double: exact product of each part, keeping both results.
Expansion expansion_scale(const Expansion& expansion, double scalar) {
    Expansion result;
    for (double component : expansion) {
        double product, error;
        two_product(component, scalar, product, error);
        if (error != 0.0) result = grow_expansion(result, error);
        if (product != 0.0) result = grow_expansion(result, product);
    }
    return result;
}

// Expansion times expansion: scale by each part of the second, and add up.
Expansion expansion_product(const Expansion& first, const Expansion& second) {
    Expansion result;
    for (double component : second) {
        result = expansion_sum(result, expansion_scale(first, component));
    }
    return result;
}

/// first - second as an expansion of at most two parts, with nothing lost.
Expansion exact_difference(double first, double second) {
    double head, tail;
    two_diff(first, second, head, tail);
    Expansion result;
    if (tail != 0.0) result.push_back(tail);
    if (head != 0.0) result.push_back(head);
    return result;
}

// Exact 2 x 2 determinant: a00 a11 - a01 a10.
Expansion determinant2(const Expansion& a00, const Expansion& a01, const Expansion& a10, const Expansion& a11) {
    return expansion_difference(expansion_product(a00, a11), expansion_product(a01, a10));
}

// Exact 3 x 3 determinant, expanded along the first row.
Expansion determinant3(const std::array<std::array<Expansion, 3>, 3>& m) {
    const Expansion first = expansion_product(m[0][0], determinant2(m[1][1], m[1][2], m[2][1], m[2][2]));
    const Expansion second = expansion_product(m[0][1], determinant2(m[1][0], m[1][2], m[2][0], m[2][2]));
    const Expansion third = expansion_product(m[0][2], determinant2(m[1][0], m[1][1], m[2][0], m[2][1]));
    return expansion_sum(expansion_difference(first, second), third);
}

// The same 3 x 3 determinant in ordinary (rounded) doubles: the fast first try.
double determinant_double3(const std::array<std::array<double, 3>, 3>& m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/// P for the rounding bound: the same expansion as the determinant, but with
/// every product made positive and every minus turned into plus. No term in
/// the sum can be larger than P, so rounding errors are at most a few eps times P.
double size_bound3(const std::array<std::array<double, 3>, 3>& m) {
    double result = 0.0;
    for (std::size_t column = 0; column < 3; ++column) {
        const std::size_t c1 = column == 0 ? 1 : 0;
        const std::size_t c2 = column == 2 ? 1 : 2;
        const double minor = std::abs(m[1][c1] * m[2][c2]) + std::abs(m[1][c2] * m[2][c1]);
        result += std::abs(m[0][column]) * minor;
    }
    return result;
}

/// Whether the fast value's sign can be trusted: |value| > gamma_n * P.
bool sign_is_certain(double value, double size_bound, double operations) {
    if (!std::isfinite(value) || !std::isfinite(size_bound)) return false;
    // A loose bound on purpose: a sign that is not trusted just takes the exact route.
    const double epsilon = std::numeric_limits<double>::epsilon();
    const double gamma = operations * epsilon / (1.0 - operations * epsilon);
    return std::abs(value) > gamma * size_bound;
}

// Whether value lies between first and second (either order, ends included).
bool between(double value, double first, double second) {
    return value >= std::min(first, second) && value <= std::max(first, second);
}

/// Whether p, already known to be on the line through a and b, lies between them.
bool within(const double* a, const double* b, const double* p) {
    return between(p[0], a[0], b[0]) && between(p[1], a[1], b[1]);
}

}  // namespace

// Orientation test: sign of (b - a) x (c - a), the 2 x 2 determinant
//   | bx - ax   by - ay |
//   | cx - ax   cy - ay |
PredicateSign orient2d(const double* a, const double* b, const double* c) {
    ensure_environment();
    require_admitted(a);
    require_admitted(b);
    require_admitted(c);
    const double bax = b[0] - a[0];
    const double bay = b[1] - a[1];
    const double cax = c[0] - a[0];
    const double cay = c[1] - a[1];
    const double fast = bax * cay - bay * cax;
    const double size_bound = std::abs(bax * cay) + std::abs(bay * cax);
    if (sign_is_certain(fast, size_bound, kOrientOperations)) return sign_of(fast);
    // Not sure: redo it exactly. Each coordinate difference becomes a
    // two-part expansion, so even the subtractions lose nothing.
    return sign_of(expansion_difference(
        expansion_product(exact_difference(b[0], a[0]), exact_difference(c[1], a[1])),
        expansion_product(exact_difference(b[1], a[1]), exact_difference(c[0], a[0]))
    ));
}

// Point-in-circle test: positive when d is inside the circle through a, b, c
// (counter-clockwise). The sign of a 3 x 3 determinant, one row per point:
//   | a-d  |a-d|^2 |
//   | b-d  |b-d|^2 |
//   | c-d  |c-d|^2 |
PredicateSign incircle(const double* a, const double* b, const double* c, const double* d) {
    ensure_environment();
    require_admitted(a);
    require_admitted(b);
    require_admitted(c);
    require_admitted(d);
    // One row per point a, b, c, measured from d: [x - dx, y - dy, (x - dx)^2 + (y - dy)^2].
    // The determinant is positive exactly when d is inside the circle through
    // a, b, c (taken counter-clockwise).
    const std::array<const double*, 3> rows = {a, b, c};
    std::array<std::array<double, 3>, 3> fast_matrix{};
    for (std::size_t row = 0; row < 3; ++row) {
        const double x = rows[row][0] - d[0];
        const double y = rows[row][1] - d[1];
        fast_matrix[row] = {x, y, x * x + y * y};
    }
    const double fast = determinant_double3(fast_matrix);
    if (sign_is_certain(fast, size_bound3(fast_matrix), kIncircleOperations)) return sign_of(fast);

    // Not sure: build the same matrix again as expansions and take the exact determinant.
    std::array<std::array<Expansion, 3>, 3> exact_matrix;
    for (std::size_t row = 0; row < 3; ++row) {
        const Expansion x = exact_difference(rows[row][0], d[0]);
        const Expansion y = exact_difference(rows[row][1], d[1]);
        exact_matrix[row] = {x, y, expansion_sum(expansion_product(x, x), expansion_product(y, y))};
    }
    return sign_of(determinant3(exact_matrix));
}

// Whether pieces ab and cd meet. If no three of the points are on one line,
// they cross exactly when c, d are on opposite sides of ab and a, b on
// opposite sides of cd. If some are on one line, they meet when such a point
// lies between the other piece's ends.
bool segments_intersect_2d(const double* a, const double* b, const double* c, const double* d) {
    const PredicateSign abc = orient2d(a, b, c);
    const PredicateSign abd = orient2d(a, b, d);
    const PredicateSign cda = orient2d(c, d, a);
    const PredicateSign cdb = orient2d(c, d, b);
    const auto zero = PredicateSign::zero;
    if (abc != zero && abd != zero && cda != zero && cdb != zero) return abc != abd && cda != cdb;
    return (abc == zero && within(a, b, c)) || (abd == zero && within(a, b, d)) ||
           (cda == zero && within(c, d, a)) || (cdb == zero && within(c, d, b));
}

}  // namespace mesher
