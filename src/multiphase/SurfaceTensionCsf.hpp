// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/SurfaceTensionCsf.hpp
// Phase  : 12.5 — Continuum Surface Force (CSF) model of Brackbill, Kothe
// & Zemach (JCP 1992) for VOF/level-set multiphase momentum coupling.
//
//   F_sigma = σ κ ∇α              (per unit volume)
//   κ       = -∇ · (∇α / |∇α|)    (interface curvature)
//
// Optional wall-adhesion contact-angle BC tilts the gradient at boundary
// faces toward the wall normal according to Young's relation.
//
// Outputs:
//   - "S_SurfaceT" (vector) body force fed into the momentum equation.
//   - "kappa"      (scalar) curvature, diagnostic.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <string>
#include <unordered_map>

namespace simall::multiphase
{

struct CsfProps
{
    double sigma = 0.072; // N/m   (water/air @ 25°C)
    std::string alphaField = "alpha";
    /// Per-wall-zone equilibrium contact angle [rad]. Missing zones default
    /// to 90° (no adhesion adjustment).
    std::unordered_map<meshing::ZoneId, double> contactAngles;
};

class SurfaceTensionCsf
{
public:
    void initialize(const meshing::Mesh& mesh, solver::FieldRegistry& fields, CsfProps props);
    void compute(solver::FieldRegistry& fields);

private:
    const meshing::Mesh* mesh_ = nullptr;
    CsfProps p_{};
};

} // namespace simall::multiphase
