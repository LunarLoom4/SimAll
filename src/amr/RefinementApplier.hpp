// =============================================================================
// SimAll Beta - AMR Subsystem
// File   : src/amr/RefinementApplier.hpp
// Phase  : 16.2 — Consume RefineFlag array produced by AdaptiveRefinement
//          and emit a new (refined) hex mesh from an existing axis-aligned
//          hex / Cartesian background.
//
// Currently supports:
//   - Cartesian / octree-leaf cells (each cell becomes 8 children on refine,
//     8 sibling cells collapse to 1 parent on coarsen if all marked).
//   - Conformal junction faces are handled by emitting a polyhedral T-face
//     (variable face node count) when a refined cell meets a coarse neighbour
//     — the FaceStorage CSR layout supports this natively (see MeshStorage.hpp).
//
// More general unstructured refinement (tet / poly cell-splitting via
// templates) is reserved for Pass 8+.
// =============================================================================
#pragma once

#include "amr/AdaptiveRefinement.hpp"
#include "meshing/MeshStorage.hpp"

#include <vector>

namespace simall::amr {

class RefinementApplier {
public:
    /// `level` carries the current refinement level per cell of `in` (use
    /// zeros if you have never refined before). On return `outLevel` carries
    /// the per-cell level of the new mesh.
    void apply(const meshing::Mesh& in,
               const std::vector<RefineFlag>& flags,
               const std::vector<std::int32_t>& level,
               meshing::Mesh& out,
               std::vector<std::int32_t>& outLevel);
};

}  // namespace simall::amr
