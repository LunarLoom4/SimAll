// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/ForcesReport.hpp
// Phase  : 17.2 — Aerodynamic / hydrodynamic forces and moments on a set
// of boundary zones (typically the wetted surface of a body). Splits
// pressure and viscous contributions and projects onto user-defined lift,
// drag, and side axes, producing engineering coefficients
//
//     C_L = F · l_hat / (½ ρ V_ref² A_ref)
//     C_D = F · d_hat / (½ ρ V_ref² A_ref)
//     C_M = M · m_hat / (½ ρ V_ref² A_ref L_ref)
//
// Reads cell-centred ρ, U, p, μ_eff from FieldRegistry. Wall shear is
// approximated by μ_eff (∂u_t/∂y)|_wall using the owner cell as the
// off-wall sample point. Suitable for both incompressible (μ_eff = μ + μ_t)
// and compressible workflows.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <string>
#include <vector>

namespace simall::io {

struct ForceReference {
    double      rhoRef     = 1.225;       // [kg/m³]
    double      velRef     = 1.0;         // [m/s]
    double      areaRef    = 1.0;         // [m²]
    double      lengthRef  = 1.0;         // [m] (for moment coefficient)
    util::Vec3d momentCenter{0, 0, 0};
    util::Vec3d liftAxis  {0, 1, 0};
    util::Vec3d dragAxis  {1, 0, 0};
    util::Vec3d sideAxis  {0, 0, 1};
};

struct ForceResult {
    util::Vec3d Fpressure  {0, 0, 0};
    util::Vec3d Fviscous   {0, 0, 0};
    util::Vec3d Mpressure  {0, 0, 0};
    util::Vec3d Mviscous   {0, 0, 0};
    double      CL = 0, CD = 0, CS = 0;
    double      CMl = 0, CMd = 0, CMs = 0;
    double      area = 0;
};

class ForcesReport {
public:
    /// Compute total forces & moments over the union of `zones`.
    /// `muEffField` is the FieldRegistry scalar name for the effective
    /// viscosity (default "mu_eff" — fall back to "mu" if not present).
    static ForceResult compute(const meshing::Mesh& mesh,
                               const solver::FieldRegistry& fields,
                               const std::vector<meshing::ZoneId>& zones,
                               const ForceReference& ref,
                               const std::string& muEffField = "mu_eff");
};

}  // namespace simall::io
