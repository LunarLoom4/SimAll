// =============================================================================
// SimAll Beta - tests/regression/RegressionFramework.hpp
// Week 19 - Regression / verification / performance / golden-hash utilities.
//
// Header-only.  Used by every test_*.cpp file under tests/regression/,
// tests/verification/, tests/performance/, tests/golden/.  Keeping it
// header-only avoids adding an extra link target and lets each regression
// executable stay small.
//
// Provides:
//   * relative_error / max_relative_error helpers
//   * FNV-1a 64-bit byte hash + double-array hash for golden checks
//   * Simple uniform-grid generator and 1-D analytic flow profiles
//   * RAII high-resolution timer for the performance harness
// =============================================================================
#pragma once

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace simall::regression
{

// ------------------------------------------------------------------ errors
[[nodiscard]] inline double relative_error(double a, double b)
{
    const double s = std::max(std::abs(a), std::abs(b));
    return s < 1e-30 ? std::abs(a - b) : std::abs(a - b) / s;
}

[[nodiscard]] inline double max_relative_error(const std::vector<double>& a,
                                               const std::vector<double>& b)
{
    double m = 0.0;
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i)
        m = std::max(m, relative_error(a[i], b[i]));
    return m;
}

[[nodiscard]] inline double l2_error(const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    const std::size_t n = std::min(a.size(), b.size());
    for (std::size_t i = 0; i < n; ++i) {
        const double d = a[i] - b[i];
        s += d * d;
    }
    return std::sqrt(s / double(std::max<std::size_t>(n, 1)));
}

// ------------------------------------------------------------------ hashing
[[nodiscard]] inline std::uint64_t fnv1a64(const void* data, std::size_t n)
{
    const auto* p = static_cast<const std::uint8_t*>(data);
    std::uint64_t h = 0xcbf29ce484222325ULL;
    for (std::size_t i = 0; i < n; ++i) {
        h ^= p[i];
        h *= 0x100000001b3ULL;
    }
    return h;
}

[[nodiscard]] inline std::uint64_t hash_doubles(const std::vector<double>& v,
                                                int quantize_digits = 6)
{
    // Quantise to avoid bit-level noise across compilers.
    std::vector<std::int64_t> q(v.size());
    const double scale = std::pow(10.0, double(quantize_digits));
    for (std::size_t i = 0; i < v.size(); ++i) {
        q[i] = std::int64_t(std::llround(v[i] * scale));
    }
    return fnv1a64(q.data(), q.size() * sizeof(std::int64_t));
}

// ------------------------------------------------------------------ grids
struct UniformGrid1D
{
    std::size_t n;
    double x0, x1;
    [[nodiscard]] double dx() const { return (x1 - x0) / double(n - 1); }
    [[nodiscard]] double x(std::size_t i) const { return x0 + double(i) * dx(); }
};

// Hagen-Poiseuille parabolic profile in a circular pipe of radius R.
// u(r) = u_max · (1 - r²/R²),  u_max = 2 u_mean.
[[nodiscard]] inline double pipe_profile(double r, double R, double uMean)
{
    return 2.0 * uMean * (1.0 - (r * r) / (R * R));
}

// Plane-Poiseuille channel profile (half-height h, mean velocity).
[[nodiscard]] inline double channel_profile(double y, double h, double uMean)
{
    return 1.5 * uMean * (1.0 - (y * y) / (h * h));
}

// Log-law for the inner region of a turbulent boundary layer.
// kappa ~ 0.41, B ~ 5.0.
[[nodiscard]] inline double log_law(double yPlus, double kappa = 0.41, double B = 5.0)
{
    return (1.0 / kappa) * std::log(std::max(yPlus, 1e-30)) + B;
}

// ------------------------------------------------------------------ timer
class HiResTimer
{
public:
    HiResTimer() { reset(); }
    void reset() { t0_ = std::chrono::high_resolution_clock::now(); }
    [[nodiscard]] double seconds() const
    {
        return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0_)
            .count();
    }

private:
    std::chrono::high_resolution_clock::time_point t0_;
};

struct PerfReport
{
    std::string name;
    double seconds = 0.0;
    std::size_t opsPerformed = 0;
    [[nodiscard]] double mops() const
    {
        return double(opsPerformed) / std::max(seconds, 1e-12) * 1e-6;
    }
};

} // namespace simall::regression
