// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/OversetInterpolation.hpp
// Phase  : 5.9 — Chimera/Overset donor-search and field interpolation.
//
// Given two independent meshes (background + foreground), this module:
//
//   1. Builds an AABB-tree over the donor mesh.
//   2. For every "receptor" cell in the target mesh, locates the donor
//      cell containing its centroid (point-in-cell test using face plane
//      half-spaces; assumes convex cells, valid for hex/tet/prism/pyr).
//   3. Computes inverse-distance weights from the receptor centroid to the
//      donor centroid and its first-ring face neighbours (Sherman-Lozano-
//      Stuart "ID3"-style stencil for arbitrary polyhedra).
//   4. Caches receptor→{donor, weights} pairs.
//   5. Applies scalar or vector interpolation on demand.
//
// Holes (cells overlapping a body region in the donor) are flagged by
// inheriting a user-supplied mask field ("cutTag") so background cells
// inside the body are excluded.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace simall::meshing {

struct OversetStencil {
    std::int64_t                 receptor = -1;
    std::array<std::int64_t, 8>  donors{-1,-1,-1,-1,-1,-1,-1,-1};
    std::array<double, 8>        weights{0,0,0,0,0,0,0,0};
    int                          count = 0;
};

class OversetInterpolation {
public:
    /// Pre-compute receptor→donor stencils. `holeMaskField` (if non-empty)
    /// is consulted on the donor mesh to skip cells marked as inside a body.
    bool build(const Mesh& donor,
               const Mesh& receptor,
               const std::string& holeMaskField,
               const solver::FieldRegistry& donorFields);

    /// Interpolate a scalar field from donor to receptor.
    /// Receptors with no donor are left untouched.
    void interpolate_scalar(const solver::ScalarField& srcDonor,
                            solver::ScalarField&       dstReceptor) const;
    void interpolate_vector(const solver::VectorField& srcDonor,
                            solver::VectorField&       dstReceptor) const;

    const std::vector<OversetStencil>& stencils() const noexcept { return stencils_; }

private:
    struct BvhNode {
        util::BoundingBox aabb;
        int left  = -1;
        int right = -1;
        int first = -1;
        int count = 0;
    };
    int  build_bvh(int first, int count);
    void cell_aabb(CellId c, util::BoundingBox& out) const;
    bool point_in_cell(const util::Vec3d& p, CellId c) const;
    CellId locate(const util::Vec3d& p) const;

    const Mesh*  donor_    = nullptr;
    const Mesh*  receptor_ = nullptr;
    std::vector<BvhNode>      bvh_;
    std::vector<int>          bvhIdx_;
    std::vector<util::BoundingBox> cellBox_;
    std::vector<OversetStencil>    stencils_;
    std::vector<std::uint8_t>      holeMask_;
};

}  // namespace simall::meshing
