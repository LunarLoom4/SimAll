// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/GEquation.hpp
// Phase  : 11.15 — Level-set G-equation for premixed turbulent combustion.
//
// The instantaneous flame surface is the iso-level G(x,t) = G_0; the field
// G is advected and propagated normal to itself at the laminar burning
// velocity s_L:
//
//   ∂G/∂t  +  (Ũ · ∇G)  =  s_T |∇G|           (Favre-mean form)
//
// where s_T is the turbulent burning velocity from Damköhler 1940:
//
//   s_T / s_L  =  1  +  C_DA · (u'_rms / s_L)
//
// Reference:
//   Peters, "Turbulent Combustion", Cambridge UP (2000), §2.4.
//
// G is reinitialised by Sussman-style iterations to remain a signed-distance
// function (|∇G| → 1).  Reinitialisation runs every reinitInterval steps.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "solver/Solver.hpp"

namespace simall::combustion
{

struct GEquationParams
{
    double s_laminar = 0.40;     // m/s  laminar flame speed
    double C_damkohler = 2.0;    // Damköhler constant in s_T law
    int reinitInterval = 5;      // outer steps between reinitialisations
    int reinitSubSteps = 10;     // pseudo-time sub-steps per reinit
    double reinitDtFactor = 0.4; // CFL fraction on local cell size
    double G_iso = 0.0;
};

class GEquation
{
public:
    void initialize(const meshing::Mesh& mesh,
                    solver::FieldRegistry& fields,
                    GEquationParams params = {});

    /// Advance G by dt using the G-equation; reads "U" (vector) and
    /// "u_rms" (scalar, sub-grid turbulent fluctuation).  Returns the
    /// L2 norm of the change in G.
    double step(double dt);

    /// Returns indicator field (0 → unburnt, 1 → burnt) by Heaviside H(G_iso-G).
    void update_progress_indicator();

    const GEquationParams& params() const noexcept { return p_; }

private:
    void reinitialise();
    double cell_size(meshing::CellId c) const;

    const meshing::Mesh* mesh_ = nullptr;
    solver::FieldRegistry* fields_ = nullptr;
    GEquationParams p_{};
    int stepCounter_ = 0;
};

} // namespace simall::combustion
