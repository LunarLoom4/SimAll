// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/CutCellHexMesher.hpp
// Phase  : 5.8 — Cartesian background-hex mesher with surface classification
// for cut-cell and immersed-boundary workflows.
//
// Given:
//   - axis-aligned bounding box [min,max]
//   - integer resolution (nx, ny, nz)
//   - a triangle soup (std::vector<Triangle>) describing the surface
//
// Produces:
//   - meshing::Mesh populated with hex cells inside-or-cut by the surface
//   - per-cell signed-distance scalar "phi" (negative inside fluid)
//   - per-cell tag "cutTag" ∈ {0=outside, 1=fluid, 2=cut}
//
// The signed-distance evaluation uses a BVH over the triangle soup for
// O(log T) queries per cell. Inside/outside classification uses parity of
// ray intersections from each cell centre along the +x axis.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <vector>

namespace simall::meshing {

struct CutTriangle {
    std::array<util::Vec3d, 3> v;
};

struct CutCellParams {
    util::Vec3d bboxMin{0, 0, 0};
    util::Vec3d bboxMax{1, 1, 1};
    int nx = 16, ny = 16, nz = 16;
    /// Discard cells whose signed distance > skinThickness above the surface
    /// (i.e. far outside). Set to large positive to keep the full background.
    double skinThickness = 0.0;
};

class CutCellHexMesher {
public:
    /// Build the cut-cell mesh. Returns false on degenerate input.
    bool build(const CutCellParams& params,
               const std::vector<CutTriangle>& tris,
               Mesh& outMesh,
               solver::FieldRegistry& outFields);

private:
    struct BvhNode {
        util::BoundingBox aabb;
        int               left  = -1;
        int               right = -1;
        int               triFirst = -1;
        int               triCount = 0;
    };
    static int build_bvh(std::vector<BvhNode>& nodes,
                         std::vector<int>& triIdx,
                         const std::vector<CutTriangle>& tris,
                         int first, int count);
    static double point_tri_dist_sq(const util::Vec3d& p, const CutTriangle& t);
    static bool ray_tri_intersect(const util::Vec3d& o, const util::Vec3d& d,
                                  const CutTriangle& t, double& tOut);
};

}  // namespace simall::meshing
