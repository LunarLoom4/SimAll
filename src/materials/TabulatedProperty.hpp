// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/TabulatedProperty.hpp
// Phase  : 15 — Tabulated property interpolation.
//
//   1-D : value(T)              — linear or monotonic-cubic (PCHIP-style)
//   2-D : value(T, p)           — bilinear over a rectilinear (T, p) grid
//
// The table x-axis must be strictly increasing.  Lookup uses std::lower_bound
// (O(log N)).  Out-of-range queries are clamped to the table endpoints unless
// allowExtrapolate=true.
// =============================================================================
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace simall::materials {

enum class TableInterp { Linear, MonotonicCubic };

class TabulatedProperty1D {
public:
    TabulatedProperty1D() = default;
    TabulatedProperty1D(std::vector<double> x, std::vector<double> y,
                        TableInterp mode = TableInterp::Linear,
                        bool allowExtrapolate = false)
        : x_(std::move(x)), y_(std::move(y)),
          mode_(mode), extrap_(allowExtrapolate)
    {
        if (x_.size() < 2 || x_.size() != y_.size())
            throw std::invalid_argument("TabulatedProperty1D: bad sizes");
        for (std::size_t i = 1; i < x_.size(); ++i)
            if (!(x_[i] > x_[i-1]))
                throw std::invalid_argument("TabulatedProperty1D: x not strictly increasing");
        if (mode_ == TableInterp::MonotonicCubic) buildPchip_();
    }

    double evaluate(double xq) const noexcept {
        if (x_.empty()) return 0.0;
        if (!extrap_) {
            if (xq <= x_.front()) return y_.front();
            if (xq >= x_.back())  return y_.back();
        }
        auto it = std::lower_bound(x_.cbegin(), x_.cend(), xq);
        std::size_t hi = (it == x_.cend()) ? x_.size() - 1
                       : (it == x_.cbegin()) ? 1
                       : static_cast<std::size_t>(it - x_.cbegin());
        std::size_t lo = hi - 1;
        const double h = x_[hi] - x_[lo];
        const double t = (xq - x_[lo]) / h;

        if (mode_ == TableInterp::Linear)
            return y_[lo] + t * (y_[hi] - y_[lo]);

        // PCHIP cubic Hermite using stored slopes m_[lo], m_[hi].
        const double t2 = t * t;
        const double t3 = t2 * t;
        const double h00 =  2.0*t3 - 3.0*t2 + 1.0;
        const double h10 =        t3 - 2.0*t2 + t;
        const double h01 = -2.0*t3 + 3.0*t2;
        const double h11 =        t3 -      t2;
        return h00*y_[lo] + h10*h*m_[lo] + h01*y_[hi] + h11*h*m_[hi];
    }

    std::size_t size() const noexcept { return x_.size(); }

private:
    /// Fritsch–Carlson monotonic cubic slope construction.
    void buildPchip_() {
        const std::size_t n = x_.size();
        std::vector<double> d(n - 1), h(n - 1);
        for (std::size_t i = 0; i + 1 < n; ++i) {
            h[i] = x_[i+1] - x_[i];
            d[i] = (y_[i+1] - y_[i]) / h[i];
        }
        m_.assign(n, 0.0);
        m_.front() = d.front();
        m_.back()  = d.back();
        for (std::size_t i = 1; i + 1 < n; ++i) {
            if (d[i-1] * d[i] <= 0.0) { m_[i] = 0.0; continue; }
            const double w1 = 2.0 * h[i] + h[i-1];
            const double w2 = h[i] + 2.0 * h[i-1];
            m_[i] = (w1 + w2) / (w1/d[i-1] + w2/d[i]);
        }
    }

    std::vector<double> x_, y_, m_;
    TableInterp         mode_   = TableInterp::Linear;
    bool                extrap_ = false;
};

/// Bilinear interpolation over a rectilinear (T, p) grid.
class TabulatedProperty2D {
public:
    TabulatedProperty2D() = default;
    /// values is row-major:  values[i*np + j] corresponds to (T[i], p[j]).
    TabulatedProperty2D(std::vector<double> T, std::vector<double> P,
                        std::vector<double> values,
                        bool allowExtrapolate = false)
        : T_(std::move(T)), P_(std::move(P)), V_(std::move(values)),
          extrap_(allowExtrapolate)
    {
        if (T_.size() < 2 || P_.size() < 2)
            throw std::invalid_argument("TabulatedProperty2D: need >=2 axis points");
        if (V_.size() != T_.size() * P_.size())
            throw std::invalid_argument("TabulatedProperty2D: values size mismatch");
        for (std::size_t i = 1; i < T_.size(); ++i)
            if (!(T_[i] > T_[i-1]))
                throw std::invalid_argument("TabulatedProperty2D: T axis not increasing");
        for (std::size_t j = 1; j < P_.size(); ++j)
            if (!(P_[j] > P_[j-1]))
                throw std::invalid_argument("TabulatedProperty2D: P axis not increasing");
    }

    double evaluate(double T, double p) const noexcept {
        if (T_.empty() || P_.empty()) return 0.0;
        const auto [it0, it1, tt] = bracket_(T_, T);
        const auto [jp0, jp1, tp] = bracket_(P_, p);
        const std::size_t np = P_.size();
        const double v00 = V_[it0 * np + jp0];
        const double v01 = V_[it0 * np + jp1];
        const double v10 = V_[it1 * np + jp0];
        const double v11 = V_[it1 * np + jp1];
        return (1-tt)*((1-tp)*v00 + tp*v01) + tt*((1-tp)*v10 + tp*v11);
    }

private:
    struct Bracket { std::size_t lo; std::size_t hi; double t; };
    Bracket bracket_(const std::vector<double>& axis, double q) const noexcept {
        if (!extrap_) {
            if (q <= axis.front()) return {0, 1, 0.0};
            if (q >= axis.back())  return {axis.size()-2, axis.size()-1, 1.0};
        }
        auto it = std::lower_bound(axis.cbegin(), axis.cend(), q);
        std::size_t hi = (it == axis.cend()) ? axis.size()-1
                       : (it == axis.cbegin()) ? 1
                       : static_cast<std::size_t>(it - axis.cbegin());
        std::size_t lo = hi - 1;
        const double t = (q - axis[lo]) / (axis[hi] - axis[lo]);
        return {lo, hi, t};
    }

    std::vector<double> T_, P_, V_;
    bool                extrap_ = false;
};

}  // namespace simall::materials
