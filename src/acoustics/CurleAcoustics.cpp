// =============================================================================
// SimAll Beta - Acoustics Subsystem
// File   : src/acoustics/CurleAcoustics.cpp
// =============================================================================
#include "acoustics/CurleAcoustics.hpp"

#include <cmath>

namespace simall::acoustics {

namespace {
constexpr double kPi    = 3.14159265358979323846;
constexpr double k4PiInv = 1.0 / (4.0 * kPi);
}  // namespace

void CurleAcoustics::push_step(double t,
                                 const std::vector<CurleSurfaceSample>& surf,
                                 const CurleObserver& obs) {
    double sum = 0.0;
    for (const auto& s : surf) {
        const double dx = obs.position[0] - s.position[0];
        const double dy = obs.position[1] - s.position[1];
        const double dz = obs.position[2] - s.position[2];
        const double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < 1e-30) continue;
        const double r  = std::sqrt(r2);
        const double ndotR = s.normal[0]*dx + s.normal[1]*dy + s.normal[2]*dz;
        // Kernel: (n · (x-y)) / r² · p · dS  (the 1/r² accounts for the
        // standard far-field projection in Curle's integral).
        sum += (ndotR / (r * r2)) * s.pressure * s.area;
    }
    times_.push_back(t);
    kernel_.push_back(sum);
}

std::vector<double> CurleAcoustics::observer_pressure_history() const {
    const std::size_t n = times_.size();
    std::vector<double> p(n, 0.0);
    if (n < 2) return p;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t ip = (i + 1 < n) ? i + 1 : i;
        const std::size_t im = (i > 0)     ? i - 1 : i;
        const double dt = times_[ip] - times_[im];
        if (dt > 1e-30) {
            p[i] = k4PiInv / c0_ * (kernel_[ip] - kernel_[im]) / dt;
        }
    }
    return p;
}

}  // namespace simall::acoustics
