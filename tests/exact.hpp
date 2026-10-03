// Used by the tests only: a slow, simple second way to get exact signs.
//
// Exact reference arithmetic for the tests, independent of src/predicates.cpp.
//
// Every finite double is an integer times a power of two, so scaling all the
// coordinates of one predicate by a common power of two turns them into
// integers. The determinants are then evaluated with arbitrary-precision
// integers, which gives their exact sign: exact rational arithmetic, done
// the simple and slow way.

#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace exact {

/// A signed integer of any size: magnitude in 32-bit limbs, least significant first.
class Integer {
public:
    Integer() = default;

    /// value * 2^shift, which must be an integer.
    static Integer from_double(double value, int shift) {
        Integer result;
        if (value == 0.0) return result;
        if (!std::isfinite(value)) throw std::invalid_argument("exact::Integer needs a finite value");
        int exponent = 0;
        const double fraction = std::frexp(std::abs(value), &exponent);
        std::uint64_t mantissa = static_cast<std::uint64_t>(std::ldexp(fraction, 53));
        int power = exponent - 53 + shift;
        while (power < 0) {
            if (mantissa & 1u) throw std::invalid_argument("exact::Integer shift leaves a fraction");
            mantissa >>= 1;
            ++power;
        }
        result.limbs_ = {static_cast<std::uint32_t>(mantissa), static_cast<std::uint32_t>(mantissa >> 32)};
        result.trim();
        result = result.shifted(power);
        result.negative_ = value < 0.0;
        return result;
    }

    int sign() const { return limbs_.empty() ? 0 : (negative_ ? -1 : 1); }

    friend Integer operator+(const Integer& a, const Integer& b) {
        if (a.negative_ == b.negative_) return with_sign(add(a.limbs_, b.limbs_), a.negative_);
        if (compare(a.limbs_, b.limbs_) >= 0) return with_sign(subtract(a.limbs_, b.limbs_), a.negative_);
        return with_sign(subtract(b.limbs_, a.limbs_), b.negative_);
    }
    friend Integer operator-(const Integer& a, const Integer& b) {
        Integer negated = b;
        negated.negative_ = !b.negative_;
        return a + negated;
    }
    friend Integer operator*(const Integer& a, const Integer& b) {
        if (a.limbs_.empty() || b.limbs_.empty()) return {};
        std::vector<std::uint32_t> product(a.limbs_.size() + b.limbs_.size(), 0);
        for (std::size_t i = 0; i < a.limbs_.size(); ++i) {
            std::uint64_t carry = 0;
            for (std::size_t j = 0; j < b.limbs_.size(); ++j) {
                const std::uint64_t t = static_cast<std::uint64_t>(a.limbs_[i]) * b.limbs_[j] + product[i + j] + carry;
                product[i + j] = static_cast<std::uint32_t>(t);
                carry = t >> 32;
            }
            product[i + b.limbs_.size()] = static_cast<std::uint32_t>(carry);
        }
        return with_sign(std::move(product), a.negative_ != b.negative_);
    }

private:
    using Limbs = std::vector<std::uint32_t>;

    static Integer with_sign(Limbs limbs, bool negative) {
        Integer result;
        result.limbs_ = std::move(limbs);
        result.trim();
        result.negative_ = negative && !result.limbs_.empty();
        return result;
    }
    void trim() {
        while (!limbs_.empty() && limbs_.back() == 0) limbs_.pop_back();
    }
    Integer shifted(int bits) const {
        if (limbs_.empty() || bits == 0) return *this;
        Limbs out(static_cast<std::size_t>(bits / 32), 0);
        const int rest = bits % 32;
        std::uint32_t carry = 0;
        for (const std::uint32_t limb : limbs_) {
            out.push_back(rest ? (limb << rest) | carry : limb);
            carry = rest ? limb >> (32 - rest) : 0;
        }
        out.push_back(carry);
        return with_sign(std::move(out), negative_);
    }
    static int compare(const Limbs& a, const Limbs& b) {
        if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
        for (std::size_t i = a.size(); i-- > 0;) {
            if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
        }
        return 0;
    }
    static Limbs add(const Limbs& a, const Limbs& b) {
        Limbs out;
        std::uint64_t carry = 0;
        for (std::size_t i = 0; i < std::max(a.size(), b.size()); ++i) {
            const std::uint64_t t = carry + (i < a.size() ? a[i] : 0) + (i < b.size() ? b[i] : 0);
            out.push_back(static_cast<std::uint32_t>(t));
            carry = t >> 32;
        }
        out.push_back(static_cast<std::uint32_t>(carry));
        return out;
    }
    static Limbs subtract(const Limbs& a, const Limbs& b) {  // |a| >= |b|
        Limbs out;
        std::int64_t borrow = 0;
        for (std::size_t i = 0; i < a.size(); ++i) {
            std::int64_t t = static_cast<std::int64_t>(a[i]) - borrow - (i < b.size() ? b[i] : 0);
            borrow = t < 0 ? 1 : 0;
            out.push_back(static_cast<std::uint32_t>(t + (borrow << 32)));
        }
        return out;
    }

    bool negative_ = false;
    Limbs limbs_;
};

/// The power of two that turns every given coordinate into an integer.
inline int common_shift(std::initializer_list<const double*> points) {
    int lowest = 0;
    bool any = false;
    for (const double* p : points) {
        for (int k = 0; k < 2; ++k) {
            if (p[k] == 0.0) continue;
            int exponent = 0;
            const double fraction = std::frexp(std::abs(p[k]), &exponent);
            const auto mantissa = static_cast<std::uint64_t>(std::ldexp(fraction, 53));
            const int low = exponent - 53 + std::countr_zero(mantissa);
            lowest = any ? std::min(lowest, low) : low;
            any = true;
        }
    }
    return -lowest;
}

/// Exact sign of det(b-a, c-a): positive when a, b, c turn counter-clockwise.
inline int orient(const double* a, const double* b, const double* c) {
    const int s = common_shift({a, b, c});
    const auto I = [s](double v) { return Integer::from_double(v, s); };
    const Integer bax = I(b[0]) - I(a[0]), bay = I(b[1]) - I(a[1]);
    const Integer cax = I(c[0]) - I(a[0]), cay = I(c[1]) - I(a[1]);
    return (bax * cay - bay * cax).sign();
}

/// Exact sign of the in-circle determinant: positive when d is strictly inside
/// the circle through counter-clockwise a, b, c.
inline int incircle(const double* a, const double* b, const double* c, const double* d) {
    const int s = common_shift({a, b, c, d});
    const auto I = [s](double v) { return Integer::from_double(v, s); };
    const Integer adx = I(a[0]) - I(d[0]), ady = I(a[1]) - I(d[1]);
    const Integer bdx = I(b[0]) - I(d[0]), bdy = I(b[1]) - I(d[1]);
    const Integer cdx = I(c[0]) - I(d[0]), cdy = I(c[1]) - I(d[1]);
    const Integer alift = adx * adx + ady * ady;
    const Integer blift = bdx * bdx + bdy * bdy;
    const Integer clift = cdx * cdx + cdy * cdy;
    return (adx * (bdy * clift - cdy * blift) - ady * (bdx * clift - cdx * blift) +
            alift * (bdx * cdy - cdx * bdy)).sign();
}

}  // namespace exact
