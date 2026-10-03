// The exact predicates: the yes/no geometric tests that steps 2 to 5 rely on.
//
// Each returns only a sign (negative, zero, positive), and that sign is always
// right, even when rounding would fool ordinary floating-point arithmetic.
// Zero means a real tie (points exactly on one line, or exactly on one
// circle); the caller decides what to do with it (see Cdt::cavity and Cdt::insert).

#pragma once

#include <cstdint>

namespace mesher {

enum class PredicateSign : std::int8_t {
    negative = -1,
    zero = 0,
    positive = 1,
};

/// Orientation test: positive when a, b, c turn counter-clockwise (c is on the
/// left of the line from a to b), negative when clockwise, zero when on one line.
PredicateSign orient2d(const double* a, const double* b, const double* c);

/// Point-in-circle test: positive when d is strictly inside the circumcircle of
/// a, b, c (given counter-clockwise), zero when exactly on it.
PredicateSign incircle(const double* a, const double* b, const double* c, const double* d);

/// Whether the line pieces ab and cd share at least one point (touching counts).
bool segments_intersect_2d(const double* a, const double* b, const double* c, const double* d);

}  // namespace mesher
