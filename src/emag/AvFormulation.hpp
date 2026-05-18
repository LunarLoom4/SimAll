// =============================================================================
// SimAll Beta - Electromagnetics Subsystem
// File   : src/emag/AvFormulation.hpp
// Week   : 18
//
// Magnetic vector potential / electric scalar potential (A-V) formulation
// for low-frequency electromagnetics in conducting media (eddy currents,
// MHD pre-coupling).  Solves
//
//      ∇ × (1/μ ∇ × A) + σ ∂A/∂t + σ ∇V = J_s
//      ∇ · (σ ∂A/∂t + σ ∇V) = 0
//
// in primitive (A, V) variables.  This file ships the *assembler* — it
// produces the per-element 4×4 (scalar V) and 3×3 (vector A) elemental
// stiffness contributions for a linear-tetrahedral discretisation.  The
// solver subsystem assembles them into a global CSR system.
//
// The element matrices are derived from
//      K_A^e = ∫_Ω (1/μ) ∇N_i · ∇N_j dV
//      M_A^e = ∫_Ω σ N_i N_j dV
//      G^e   = ∫_Ω σ N_i ∇N_j dV
// using closed-form integration on a linear tetrahedron.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace simall::emag {

struct Tet4 {
    std::array<std::uint32_t, 4>           nodes;
    std::array<std::array<double, 3>, 4>   coords;
    double                                   mu  = 1.25663706e-6;   // permeability
    double                                   sigma = 1.0;              // conductivity
};

struct ElementMatrices {
    std::array<std::array<double, 4>, 4>   K_V;        // ∫ σ ∇N·∇N dV
    std::array<std::array<double, 4>, 4>   M_A;        // ∫ σ N N dV (for ∂A/∂t)
    std::array<std::array<double, 4>, 4>   K_A;        // ∫ (1/μ) ∇N·∇N dV
    double                                  volume = 0.0;
};

[[nodiscard]] ElementMatrices assemble_tet4(const Tet4& e);

/// Convenience: build the V-only Laplace matrix (DC current conduction)
/// for an entire tetrahedral mesh.  Returns CSR (rowPtr, colIdx, values)
/// and the per-row diagonal index for quick Jacobi preconditioning.
struct CsrFromTets {
    std::size_t              n = 0;
    std::vector<std::size_t> rowPtr;
    std::vector<std::size_t> colIdx;
    std::vector<double>      values;
};
[[nodiscard]] CsrFromTets build_av_scalar_system(const std::vector<Tet4>& mesh,
                                                  std::size_t               nNodes);

}  // namespace simall::emag
