// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/DriftFlux.hpp
// Phase  : 12.8 — Drift-Flux model (Zuber & Findlay 1965) for two-phase
// dispersed flows (gas bubbles in liquid, oil droplets in water, etc.).
//
//   j   =  α_g U_g  +  α_l U_l                              (total volumetric flux)
//   j_g =  α_g U_g
//   U_gj=  U_g - j                                          (drift velocity)
//   ⟨j_g⟩=  C_0 ⟨α_g⟩ ⟨j⟩  +  ⟨α_g⟩ ⟨U_gj⟩                  (Z-F constitutive law)
//
// Z-F constants for vertical bubbly flow with Wallis correlation:
//   C_0   = 1.2                                              (distribution param)
//   U_gj  = 1.41 ( σ g (ρ_l-ρ_g) / ρ_l² )^{1/4}              (Harmathy)
//
// Reads:  fields["U_mix"], fields["alpha_g"]
// Writes: fields["U_drift_g"] (= U_gj),
//         fields["U_g"]       (= U_mix + (1-α_g)/ρ_m·...) reconstructed phase velocity.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

namespace simall::multiphase
{

struct DriftFluxProps
{
    double rhoLiquid = 998.2;
    double rhoGas = 1.225;
    double surfaceTension = 0.0728; // N/m (water-air)
    double gravity = 9.81;
    double C0 = 1.2; // distribution parameter (Zuber-Findlay)
    // Drift velocity model:
    enum class DriftLaw
    {
        Harmathy,
        IshiiChawla,
        UserConstant
    };
    DriftLaw driftLaw = DriftLaw::Harmathy;
    double Ugj_const = 0.0; // used if driftLaw == UserConstant
};

class DriftFlux
{
public:
    void initialize(const meshing::Mesh& mesh, DriftFluxProps props);

    /// One outer iteration: computes drift-velocity field and reconstructs
    /// phase velocities. Returns max |U_gj|.
    double apply(solver::FieldRegistry& fields);

    const DriftFluxProps& props() const noexcept { return p_; }

private:
    const meshing::Mesh* mesh_ = nullptr;
    DriftFluxProps p_{};
};

} // namespace simall::multiphase
