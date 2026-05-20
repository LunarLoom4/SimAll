// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/SnappyHexMesher.hpp
// Phase  : 23 Pass 10
//
// SnappyHexMesher orchestrates the three classic "snap-from-STL" stages
// using SimAll's existing building blocks:
//
//    1. CASTELLATION
//         OctreeMesher generates a watertight stepped hex background mesh
//         from the bounding box of the input STL surface.  Cells inside the
//         closed surface (and BOUNDARY leaves intersected by triangles) are
//         retained.
//
//    2. SNAPPING
//         Boundary nodes of the castellated hex mesh are projected onto the
//         nearest point on the STL.  Projection is iterated with a relax
//         factor so a node is moved a fraction of the remaining distance
//         per iteration, producing a smooth convergence onto the surface
//         without inverting interior cells.
//
//    3. LAYER ADDITION   (opt.nLayers > 0)
//         After snapping, each boundary quad face is triangulated and
//         passed to the existing PrismLayerExtruder, which emits an
//         independent prism-layer Mesh.  The returned `prismLayerOut`
//         mesh is intentionally NOT stitched into `backgroundOut` here;
//         the caller couples them via the `mergeMeshes` utility (P23
//         work-item #2) once that lands.  This keeps the snap pipeline
//         pure (no in-place topology mutation of the hex mesh) and lets
//         users opt out of layer addition entirely.
//
// All geometry is general 3-D.  No symmetry / 2-D specialisation.
// =============================================================================
#pragma once

#include "MeshStorage.hpp"
#include "StlImporter.hpp"

#include <cstddef>

namespace simall::meshing {

struct SnappyHexOptions {
    // --- Castellation ------------------------------------------------------
    int    maxDepth          = 6;     // forwarded to OctreeMeshOptions
    int    minDepth          = 3;

    // --- Snapping ----------------------------------------------------------
    bool   enableSnapping    = true;
    int    nSnapIters        = 5;
    /// Maximum distance a node may snap from its castellated position,
    /// expressed as a multiple of the finest cell edge length.  Nodes
    /// further than this from the STL are left alone (prevents wild
    /// projections through thin features).
    double snapMaxDistFrac   = 2.0;

    // --- Layer addition ----------------------------------------------------
    int    nLayers           = 0;     // 0 -> skip the layer-addition phase
    double firstLayerHeight  = 1.0e-4;
    double layerGrowthRatio  = 1.2;
};

struct SnappyHexStats {
    std::size_t backgroundCells   = 0;
    std::size_t backgroundNodes   = 0;
    std::size_t boundaryNodes     = 0;   // candidate snap nodes
    std::size_t snappedNodes      = 0;   // actually moved (within budget)
    int         snapItersRun      = 0;
    std::size_t prismLayerCells   = 0;
    std::size_t prismLayerNodes   = 0;
};

class SnappyHexMesher {
public:
    /// Run the full pipeline.  `backgroundOut` receives the castellated +
    /// snapped hex mesh.  `prismLayerOut` receives the prism-layer mesh
    /// (empty if `opt.nLayers == 0`).  Returns per-stage statistics.
    SnappyHexStats mesh(const StlSurface&  surface,
                        SnappyHexOptions   opt,
                        Mesh&              backgroundOut,
                        Mesh&              prismLayerOut);

    /// Castellation-only convenience entry point (no snapping, no layers).
    SnappyHexStats mesh_background_only(const StlSurface& surface,
                                         SnappyHexOptions  opt,
                                         Mesh&             backgroundOut);
};

}  // namespace simall::meshing
