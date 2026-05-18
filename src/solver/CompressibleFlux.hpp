// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/CompressibleFlux.hpp
// Phase  : 6.8 — Density-based compressible flux for the inviscid
// Euler/Navier-Stokes equations:
//
//   ∂U/∂t + ∇·F(U) = 0
//   U = (ρ, ρu, ρv, ρw, ρE)ᵀ
//
// Implements the **HLLC** approximate Riemann solver (Toro 1994/2009) for
// arbitrary-orientation polyhedral faces. HLLC restores the contact and
// shear waves that HLL averages away, which is essential for viscous
// boundary layers in compressible flow.
//
// The flux routine is mesh-agnostic — given left/right primitive states
// (ρ, u, v, w, p) and the face normal, it returns a 5-component conservative
// flux. The caller (CompressibleAlgorithm) handles second-order MUSCL
// reconstruction with a Venkatakrishnan limiter and assembles the residual.
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"
#include <array>

namespace simall::solver {

struct PrimState {
    double rho = 1.0;
    double u   = 0.0;
    double v   = 0.0;
    double w   = 0.0;
    double p   = 1.0e5;
};

/// Conservative state U = (ρ, ρu, ρv, ρw, ρE).
using ConsFlux = std::array<double, 5>;

class CompressibleFlux {
public:
    explicit CompressibleFlux(double gamma = 1.4) : gamma_(gamma) {}

    /// HLLC numerical flux through a face with outward unit normal n and
    /// area magnitude A. The returned flux is already multiplied by A.
    ConsFlux hllc(const PrimState& L, const PrimState& R,
                  const util::Vec3d& n, double A) const;

    /// Rotate primitive velocity (u,v,w) into face-aligned frame so the
    /// shock-normal axis is x. Used by 1-D Riemann inside hllc().
    static void rotate_to_face(const util::Vec3d& n, double& u, double& v, double& w,
                               double& tx_x, double& tx_y, double& tx_z,
                               double& ty_x, double& ty_y, double& ty_z);

    double gamma() const noexcept { return gamma_; }

private:
    double gamma_ = 1.4;
};

}  // namespace simall::solver
