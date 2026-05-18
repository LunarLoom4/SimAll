// =============================================================================
// SimAll Beta - Adjoint Subsystem
// File   : src/adjoint/ContinuousAdjoint.hpp
// Week   : 18
//
// Continuous adjoint for shape sensitivity of incompressible CFD problems.
// Approximates the surface sensitivity field  G(x) = (τψ · n - p_ψ) · û
// from primal flow (u, p) and adjoint flow (uψ, pψ) on a closed boundary
// patch.  Used to drive gradient-based shape optimisation (FFD / morphing).
//
// This module does NOT solve the adjoint Navier-Stokes equations — the
// solver subsystem does that and feeds the converged adjoint fields here.
// What this module provides:
//
//   * `assemble_shape_sensitivity()` — integrates the surface kernel
//     over a list of design-boundary faces, returns per-node grad-J;
//   * `assemble_volume_sensitivity()` — analogous volume integral for
//     topology optimisation (drives DensityMethod elsewhere);
//   * `accumulate_node_sensitivity()` — face → node area-weighted scatter.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::adjoint {

using NodeIdx = std::uint32_t;

struct SurfaceFace {
    std::array<NodeIdx, 4>    nodes;              // tri uses [0..2], quad uses [0..3]
    std::uint8_t              nNodes = 3;
    std::array<double, 3>     normal{0,0,1};      // outward unit normal
    double                    area   = 0.0;
};

struct AdjointFlowSample {
    std::array<double, 3> uPsi{};                  // adjoint velocity
    double                pPsi = 0.0;              // adjoint pressure
    double                mu   = 0.0;              // dynamic viscosity at face
};

/// Per-face shape sensitivity (scalar: dJ/dx_n along the local normal).
[[nodiscard]] std::vector<double> assemble_shape_sensitivity(
        const std::vector<SurfaceFace>&        faces,
        const std::vector<AdjointFlowSample>&  perFace);

/// Volume kernel for topology optimisation density α(x):  s(x) = -ψ_p · ∇p
/// approximated by face-pressure jumps over a control volume.
[[nodiscard]] std::vector<double> assemble_volume_sensitivity(
        const std::vector<double>& pressureCell,
        const std::vector<double>& adjointPressureCell);

/// Scatter face sensitivities onto vertex/node sensitivities, weighting by
/// face area so the resulting per-node gradient is suitable for an
/// FFD / morpher step.
[[nodiscard]] std::vector<double> accumulate_node_sensitivity(
        std::size_t                            nNodes,
        const std::vector<SurfaceFace>&        faces,
        const std::vector<double>&             faceSens);

}  // namespace simall::adjoint
