// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/MixtureModel.hpp
// Phase  : 12.8 — Algebraic Slip Mixture model (Manninen, Taivassalo &
// Kallio 1996).  Single momentum equation in the mixture velocity, with
// secondary-phase velocities reconstructed via an algebraic-slip relation.
//
//   ρ_m = Σ α_k ρ_k                    (mixture density)
//   μ_m = Σ α_k μ_k                    (mixture viscosity)
//   U_m = Σ α_k ρ_k U_k / ρ_m         (mass-weighted mixture velocity)
//   U_kr = U_k - U_m                  (relative velocity of phase k)
//
//   Algebraic slip (drag balance against pressure-gradient + gravity):
//       U_kr  =  τ_p · (ρ_k - ρ_m) / ρ_k  · ( g - (U_m·∇)U_m - ∂U_m/∂t )
//       τ_p   =  ρ_k d²_k / (18 μ_m f_d)     (Stokes relaxation time)
//       f_d   =  1 + 0.15 Re_p^0.687         (Schiller-Naumann drag)
//
// Reads:  fields["U_mix"], fields["alpha_k"] (k=1..N-1; α_0 = 1-Σ)
// Writes: fields["rho_mix"], fields["mu_mix"], fields["U_kr_k"] per phase,
//         and slip-induced drift flux contributions for the α-equations.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <string>
#include <vector>

namespace simall::multiphase
{

struct MixturePhase
{
    std::string name;
    double rho;      // kg/m³
    double mu;       // Pa·s
    double diameter; // m   (representative particle/bubble dia)
};

class MixtureModel
{
public:
    void initialize(const meshing::Mesh& mesh,
                    std::vector<MixturePhase> phases,
                    double gx = 0.0,
                    double gy = 0.0,
                    double gz = -9.81);

    /// One outer iteration of mixture-property update + algebraic-slip
    /// velocity reconstruction.  Returns max |U_kr| observed (diagnostic).
    double apply(solver::FieldRegistry& fields);

    std::size_t num_phases() const noexcept { return phases_.size(); }
    const std::vector<MixturePhase>& phases() const noexcept { return phases_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    std::vector<MixturePhase> phases_;
    double g_[3]{0, 0, -9.81};
};

} // namespace simall::multiphase
