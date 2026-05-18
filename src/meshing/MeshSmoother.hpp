// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshSmoother.hpp
// Phase  : 4.7 — Geometric mesh smoothers for repair after morphing /
// adaptation. Implemented variants:
//
//   - Laplace        : x_i ← (1/N) Σ x_j        (neighbours of node i)
//   - LengthWeighted : x_i ← Σ (w_j x_j) / Σ w_j with w_j = 1/||x_i - x_j||
//   - Optimization   : minimize Σ (1/μ_face) where μ_face is a per-face
//                      shape-quality metric (face-area / edge-length²); a
//                      few gradient-descent steps per pass.
//
// Boundary nodes are pinned by default; an optional "boundaryFreeZones"
// list allows sliding boundary motion (projection back to its mesh plane
// is intentionally NOT done here — that is delegated to the CAD layer).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include <cstdint>
#include <vector>

namespace simall::meshing {

enum class SmoothMode { Laplace, LengthWeighted, Optimization };

struct SmoothParams {
    SmoothMode  mode      = SmoothMode::Laplace;
    int         iterations = 5;
    double      relax     = 0.5;          // 0..1 under-relaxation
    /// If set, only nodes whose every face stays in any of these zones
    /// (or interior) are moved — used to lock pre-defined wall patches.
    std::vector<ZoneId> pinnedBoundaryZones;
};

class MeshSmoother {
public:
    /// Smooth `mesh` nodes in place. Recomputes cell volumes & face areas
    /// at the end; returns true if all cells remain non-degenerate.
    static bool smooth(Mesh& mesh, const SmoothParams& p);

private:
    static void build_node_neighbours(const Mesh& m,
                                      std::vector<std::int32_t>& offsets,
                                      std::vector<NodeId>&       indices);
};

}  // namespace simall::meshing
