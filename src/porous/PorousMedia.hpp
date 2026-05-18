// =============================================================================
// SimAll Beta - Porous Media Subsystem
// File   : src/porous/PorousMedia.hpp
// Phase  : 12.5 — Volume-averaged Darcy-Forchheimer model for flow through
//          packed beds, filters, catalytic reactors, radiators, foams …
//
// Momentum sink in cells flagged as porous:
//
//     S_porous = -(μ D + ½ ρ |U| F) · U     [N/m³]
//
//   D : viscous resistance tensor  [1/m²]
//   F : inertial resistance tensor [1/m]
//
// Both can be ANISOTROPIC (diagonal in a local frame); a per-zone rotation
// brings the local D, F into the global Cartesian frame.
//
// The source is written into FieldRegistry vector "S_porous" — the solver
// momentum assembly already supports vector source channels.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"
#include "utilities/MathTypes.hpp"

#include <cstdint>
#include <vector>

namespace simall::porous {

/// Diagonal D / F in a local Cartesian frame (e1,e2,e3 → orthonormal).
struct PorousZone {
    meshing::ZoneId zone   = 0;             // mesh cell zone (must already exist)
    util::Vec3d     dDiag  {0,0,0};         // viscous resistance per-axis [1/m²]
    util::Vec3d     fDiag  {0,0,0};         // inertial resistance per-axis [1/m]
    util::Vec3d     axis1  {1,0,0};
    util::Vec3d     axis2  {0,1,0};
    util::Vec3d     axis3  {0,0,1};
};

class PorousMedia {
public:
    void add_zone(PorousZone z) { zones_.push_back(z); }
    void clear()                { zones_.clear(); }

    /// Build per-cell tensors from the supplied zone definitions.
    /// `cellZone` is a per-cell zone-id array (size = nCells); typically
    /// produced by the mesh-builder or by ImmersedBoundary tagging.
    void initialize(const meshing::Mesh& mesh,
                    const std::vector<meshing::ZoneId>& cellZone);

    /// Update fields["S_porous"] from current fields["U"] and fluid props.
    /// rho is the working density (or mixture-averaged value); mu the
    /// dynamic viscosity.
    void apply(solver::FieldRegistry& fields, double rho, double mu) const;

private:
    const meshing::Mesh*    mesh_ = nullptr;
    std::vector<PorousZone> zones_;
    // Per-cell 3x3 matrices D and F stored as 9 doubles each (row-major).
    util::aligned_vector<double> D_;
    util::aligned_vector<double> F_;
    std::vector<std::uint8_t>    active_;  // 1 if cell belongs to any zone
};

}  // namespace simall::porous
