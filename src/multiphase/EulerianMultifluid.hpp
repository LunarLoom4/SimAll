// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/EulerianMultifluid.hpp
// Phase  : 12.9 — N-phase Eulerian Multifluid model.
//
// Each phase k ∈ {0,…,N-1} carries its own:
//   - volume fraction  α_k
//   - velocity field   U_k
//   - temperature      T_k
//   - mass-conservation source ṁ_kl (from external phase-change models)
//
// Per-phase continuity:
//   ∂(α_k ρ_k)/∂t + ∇·(α_k ρ_k U_k) = Σ_l ṁ_lk
//
// Per-phase momentum:
//   ∂(α_k ρ_k U_k)/∂t + ∇·(α_k ρ_k U_k ⊗ U_k)
//     = -α_k ∇p + ∇·(α_k τ_k) + α_k ρ_k g + M_k
//   where  M_k = Σ_l K_kl (U_l - U_k)        (interphase drag transfer)
//
// Per-phase energy:
//   ∂(α_k ρ_k h_k)/∂t + ∇·(α_k ρ_k U_k h_k)
//     = α_k Dp/Dt + ∇·(α_k λ_k ∇T_k) + Σ_l h_kl A_kl (T_l - T_k) + L ṁ_kl
//
// Drag closure: Schiller-Naumann for spheres, Ishii-Zuber for swarms.
// Each phase's α-equation solved by ScalarTransport; momentum and energy
// equations solved with the standard solver::SimpleAlgorithm/PISO loop on
// per-phase fields registered as "U_<name>", "T_<name>", etc.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"
#include "solver/LinearSolvers.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <string>
#include <vector>

namespace simall::multiphase {

struct EulerianPhase {
    std::string name;
    double rho;                     // kg/m³
    double mu;                      // Pa·s
    double cp;                      // J/(kg K)
    double k_thermal;               // W/(m K)
    double diameter;                // m   (representative particle/bubble)
    bool   dispersed = true;        // false → continuous carrier phase
};

struct DragPair {
    std::size_t i, j;               // phase indices
    double      K;                  // last computed K_ij = ρ_m C_d (3/4) |U_r| / d_p
};

class EulerianMultifluid {
public:
    EulerianMultifluid(meshing::Mesh& mesh,
                       solver::FieldRegistry& fields,
                       solver::ILinearSolver& linear);

    void configure(std::vector<EulerianPhase> phases,
                   const std::vector<solver::BoundarySpec>& bcs,
                   double gx = 0.0, double gy = 0.0, double gz = -9.81);

    /// Per-time-step coupled iteration.  Solves α_k transport equations,
    /// updates interphase drag coefficients K_ij, applies them as explicit
    /// momentum source terms on each per-phase velocity, and re-evaluates
    /// per-phase thermal-energy source/sink couplings.  The full velocity
    /// linear solve is delegated to the parent SimpleAlgorithm/PISO driver
    /// per-phase (this module wires the closure couplings).
    void step(double dt);

    std::size_t num_phases() const noexcept { return phases_.size(); }
    const std::vector<EulerianPhase>& phases() const noexcept { return phases_; }
    const std::vector<DragPair>&      drag_pairs() const noexcept { return drag_; }

private:
    void update_mixture_density();
    void update_drag_coefficients();
    void apply_interphase_momentum();
    void apply_interphase_energy();

    meshing::Mesh&             mesh_;
    solver::FieldRegistry&     F_;
    solver::ILinearSolver&     lin_;
    std::vector<EulerianPhase> phases_;
    std::vector<solver::BoundarySpec> bcs_;
    double g_[3]{0,0,-9.81};

    // Per-phase α-transport equations.
    std::vector<std::unique_ptr<solver::ScalarTransport>> alphaEq_;
    // Interphase drag coefficients K_ij (one per ordered pair i<j).
    std::vector<DragPair>      drag_;
};

}  // namespace simall::multiphase
