// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/PrismLayerExtruder.hpp
// Phase  : 5.4 — boundary-layer prism extrusion.
//
// Given a triangulated wall surface (boundary face zone) and a polyhedral
// volume mesh, extrudes N anisotropic prism layers normal to the wall before
// transitioning to the volumetric (tet/poly) mesh. Layer thicknesses follow
// a geometric progression:
//
//   t_i = t_0 · r^i  ,    r = growth ratio (typ. 1.2)
//   t_0 = y+ μ / (ρ u_τ)  (Phase 5.4 first-layer-from-y+ formula)
//
// Each input wall triangle becomes a stack of (N) wedges (3-node prism).
// The new prism cells are inserted INTO the existing Mesh via
// ConnectivityBuilder, preserving the polyhedral data layout.
//
// Fully 3-D, no symmetry assumptions, supports arbitrarily curved walls.
// =============================================================================
#pragma once

#include "MeshStorage.hpp"

#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace simall::meshing
{

struct PrismLayerOptions
{
    int nLayers = 5;
    double firstLayerHeight = 1.0e-4;
    double growthRatio = 1.2;
};

class PrismLayerExtruder
{
public:
    /// Generate N anisotropic prism layers offset INWARD from a wall surface.
    /// `wallTris` are triangle vertex indices into `nodes`. The result is
    /// stored as a new Mesh whose first N · |wallTris| cells are the prism
    /// stack; downstream meshers concatenate this with the volumetric mesh.
    void extrude(const std::vector<util::Vec3d>& wallNodes,
                 const std::vector<std::array<std::uint32_t, 3>>& wallTris,
                 const std::vector<util::Vec3d>& outwardNormalsPerNode,
                 PrismLayerOptions opt,
                 Mesh& out);

    /// First-layer height from y+ for a known shear velocity.
    static double first_layer_from_yplus(double yPlus, double nu, double uTau)
    {
        return yPlus * nu / std::max(uTau, 1e-30);
    }
};

} // namespace simall::meshing
