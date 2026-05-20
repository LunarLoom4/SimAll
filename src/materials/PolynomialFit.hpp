// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/PolynomialFit.hpp
// Phase  : 15 — Temperature-dependent property representation.
//
//   PolynomialFit<N>::evaluate(T) = a0 + a1 T + a2 T^2 + ... + aN T^N
//
// Polynomial form chosen to match ANSYS-Fluent / STAR-CCM+ convention for
// piecewise property curves over [Tmin, Tmax].  Outside the validity window
// the evaluator clamps to the endpoint values (configurable).
// =============================================================================
#pragma once

#include <array>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace simall::materials
{

class PolynomialFit
{
public:
    PolynomialFit() = default;
    PolynomialFit(std::vector<double> coeffs, double Tmin, double Tmax, bool clampOutside = true)
        : c_(std::move(coeffs)), tmin_(Tmin), tmax_(Tmax), clamp_(clampOutside)
    {
        if (c_.empty())
            throw std::invalid_argument("PolynomialFit: empty coefficients");
        if (!(tmax_ > tmin_))
            throw std::invalid_argument("PolynomialFit: Tmax must exceed Tmin");
    }

    /// Horner-form evaluation.
    double evaluate(double T) const noexcept
    {
        const double t = clamp_ ? (T < tmin_ ? tmin_ : (T > tmax_ ? tmax_ : T)) : T;
        double y = c_.back();
        for (std::size_t k = c_.size() - 1; k-- > 0;)
            y = y * t + c_[k];
        return y;
    }

    /// d/dT (Horner-style).
    double derivative(double T) const noexcept
    {
        if (c_.size() < 2)
            return 0.0;
        const double t = clamp_ ? (T < tmin_ ? tmin_ : (T > tmax_ ? tmax_ : T)) : T;
        const std::size_t n = c_.size() - 1;
        double y = c_[n] * static_cast<double>(n);
        for (std::size_t k = n - 1; k > 0; --k)
            y = y * t + c_[k] * static_cast<double>(k);
        return y;
    }

    double tmin() const noexcept { return tmin_; }
    double tmax() const noexcept { return tmax_; }
    std::size_t order() const noexcept { return c_.empty() ? 0 : c_.size() - 1; }
    const std::vector<double>& coefficients() const noexcept { return c_; }

private:
    std::vector<double> c_;
    double tmin_ = 0.0;
    double tmax_ = 0.0;
    bool clamp_ = true;
};

} // namespace simall::materials
