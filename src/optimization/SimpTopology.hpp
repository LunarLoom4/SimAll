// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/SimpTopology.hpp
// Phase  : 13 — SIMP (Solid Isotropic Material with Penalization)
// density-based topology optimisation with Helmholtz PDE filter and
// Optimality-Criteria (OC) update.
//
// Design variable ρ ∈ [0,1] per cell, penalised effective property
// K_eff(ρ) = K_min + (K_max − K_min) ρ^p, p ∈ [1, 5].
//
// Filtering: Helmholtz PDE   −r² ∇²ρ̃ + ρ̃ = ρ
// Volume constraint: Σ ρ̃_e V_e ≤ V_target
// Update: classical OC with bisection on Lagrange multiplier.
//
// Sensitivities (∂J/∂ρ_e) are supplied by the caller (typically from the
// adjoint solver). This module provides the filter + bound-handling +
// projection-update without prescribing the physics.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstddef>
#include <memory>

namespace simall::solver { class ILinearSolver; struct LinearSolverConfig; }

namespace simall::optimization {

struct SimpConfig {
    double penalty       = 3.0;     ///< SIMP p
    double K_min         = 1e-9;    ///< Ersatz-material lower bound
    double K_max         = 1.0;     ///< Solid-phase property
    double filterRadius  = 0.0;     ///< Helmholtz r (0 ⇒ disable filtering)
    double volumeFraction= 0.4;     ///< Σρ̃·V / Vtot ≤ this
    double move          = 0.2;     ///< OC move-limit
    double damping       = 0.5;     ///< OC damping exponent η
    double rhoMin        = 1e-3;    ///< Lower density bound
};

class SimpTopology {
public:
    /// Initialise design field. Creates registry fields:
    ///   "rho"      — raw design variable
    ///   "rhoTilde" — filtered density
    ///   "K_simp"   — penalised property field
    bool initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    const SimpConfig& cfg);

    /// One design iteration:
    ///   1. apply Helmholtz filter            : rho → rhoTilde
    ///   2. update penalised property         : K_simp = K_min+(K_max-K_min) rhoTilde^p
    ///   3. consume external sensitivity dC   : "dC_drho"
    ///   4. apply chain-rule to filtered      : dC_drhoTilde
    ///   5. OC update with bisection          : new rho
    /// Returns the L1 design change ‖Δρ‖₁ / N.
    double step(const meshing::Mesh& mesh, solver::FieldRegistry& fields);

    const SimpConfig& config() const noexcept { return cfg_; }

private:
    void apply_helmholtz_filter(const meshing::Mesh& mesh,
                                solver::FieldRegistry& fields);

    SimpConfig cfg_{};
    const meshing::Mesh* mesh_ = nullptr;
    std::unique_ptr<solver::ILinearSolver> filterSolver_;
};

}  // namespace simall::optimization
