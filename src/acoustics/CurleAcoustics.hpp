// =============================================================================
// SimAll Beta - Acoustics Subsystem
// File   : src/acoustics/CurleAcoustics.hpp
// Week   : 18
//
// Curle's acoustic analogy (Curle, 1955).  Computes the far-field pressure
// signature radiated by *low-Mach* surface dipoles — the standard model
// for solid-body broadband noise (cylinder/jet/airfoil trailing-edge
// noise) when retarded-time differences across the body are negligible.
//
//      p'(x, t) =  1/(4π c)  ·  ∂/∂t  ∫_S (n_j p · (x_j - y_j) / r²) dS
//
//   * x  observer location
//   * y  surface point
//   * r  = |x - y|
//   * p  surface pressure fluctuation
//   * n  outward surface normal
//
// This module accumulates the integrand over surface-pressure samples
// emitted by the solver (typically every wall face every step) and
// returns the observer time series after temporal differencing.
//
// FW-H (Ffowcs-Williams & Hawkings, 1969) — the monopole + dipole +
// Lighthill-stress generalisation — already lives in
// `acoustics/FwhSurface.{hpp,cpp}`.  Curle is the dipole-only limit of
// FW-H, useful for incompressible-solver post-processing.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace simall::acoustics {

struct CurleSurfaceSample {
    std::array<double, 3> position;
    std::array<double, 3> normal;
    double                pressure = 0.0;
    double                area     = 0.0;
};

struct CurleObserver {
    std::array<double, 3> position;
};

class CurleAcoustics {
public:
    /// `soundSpeed` is the ambient c_0 used in the leading 1/(4π c) factor.
    explicit CurleAcoustics(double soundSpeed = 343.0) : c0_(soundSpeed) {}

    /// Append one time-step snapshot of the surface pressure field for the
    /// chosen observer.  After all steps are supplied, call
    /// `observer_pressure_history()` to get the differenced p'(t).
    void push_step(double                             time,
                    const std::vector<CurleSurfaceSample>& surface,
                    const CurleObserver&               observer);

    /// Returns p'(t) = (1/(4πc)) · ∂/∂t  of the integrated kernel.
    [[nodiscard]] std::vector<double> observer_pressure_history() const;

    [[nodiscard]] const std::vector<double>& times() const noexcept { return times_; }

    void clear() { times_.clear(); kernel_.clear(); }

private:
    double              c0_;
    std::vector<double> times_;
    std::vector<double> kernel_;     // raw kernel value before time derivative
};

}  // namespace simall::acoustics
