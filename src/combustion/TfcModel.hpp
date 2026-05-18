// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/TfcModel.hpp
// Phase  : 11.16 — Turbulent Flame-speed Closure (TFC) of Zimont (1979).
//
// Transports a progress variable c ∈ [0, 1] (0 = fresh, 1 = burnt):
//
//   ∂(ρ̄ c̃)/∂t + ∇·(ρ̄ Ũ c̃) = ∇·(ρ̄ D_t ∇c̃) + ρ_u · U_t · |∇c̃|
//
// where U_t is the Zimont turbulent burning velocity:
//
//   U_t = A · u'^{3/4} · S_L^{1/2} · α_u^{-1/4} · L_t^{1/4}
//
// with:
//   A     = 0.52  (model constant)
//   u'    = √(2 k/3)
//   S_L   = laminar flame speed
//   α_u   = unburnt thermal diffusivity = λ/(ρ_u c_p)
//   L_t   = integral length scale ≈ k^{3/2}/ε
//
// The reaction-rate source is ω̇_c = ρ_u U_t |∇c|.
//
// Reference: Zimont, "Theory of turbulent combustion of a homogeneous
//            fuel mixture at high Reynolds numbers", Combust. Explos.
//            Shock Waves 15, 305-311 (1979).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/ScalarTransport.hpp"
#include "solver/Solver.hpp"

#include <memory>
#include <vector>

namespace simall::combustion {

struct TfcProps {
    double A         = 0.52;
    double S_L       = 0.40;        // laminar flame speed [m/s]
    double rho_u     = 1.18;        // unburnt density   [kg/m³]
    double rho_b     = 0.18;        // burnt density
    double alpha_u   = 2.0e-5;      // unburnt thermal diff. [m²/s]
    double H_combust = 50.0e6;      // J/kg fuel LHV
    double D_t       = 1.0e-4;      // turbulent diffusivity for c-equation
};

class TfcModel {
public:
    TfcModel(meshing::Mesh& mesh,
             solver::FieldRegistry& fields,
             solver::ILinearSolver& linear);

    void configure(TfcProps props,
                   const std::vector<solver::BoundarySpec>& bcs);

    /// One outer iteration: compute U_t, build ω̇_c source, transport c.
    /// Reads "k", "epsilon" (Favre-mean turbulence) and writes "S_progress",
    /// "S_combustion", "rho_mix", "T_combust".
    void step();

    const TfcProps& props() const noexcept { return p_; }

private:
    meshing::Mesh&                         mesh_;
    solver::FieldRegistry&                 F_;
    solver::ILinearSolver&                 lin_;
    TfcProps                               p_{};
    std::vector<solver::BoundarySpec>      bcs_;
    std::unique_ptr<solver::ScalarTransport> cEq_;
};

}  // namespace simall::combustion
