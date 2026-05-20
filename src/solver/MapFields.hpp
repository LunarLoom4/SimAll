// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/MapFields.hpp
// Phase  : 23 Pass 12
//
// Cross-mesh field mapper. Interpolates every scalar and vector field stored
// in a source FieldRegistry from a source Mesh onto a target Mesh, writing
// the result into a destination FieldRegistry sized to the target mesh's
// cell count. Two strategies are supported:
//
//   * MapMethod::Nearest          - copy the value from the source cell whose
//                                   centroid is closest to the target cell
//                                   centroid (zero-order, conservative-ish,
//                                   defensive default).
//
//   * MapMethod::InverseDistance  - inverse-distance weighting using the
//                                   k nearest source cells (1/d^p weights).
//                                   First-order accurate on smooth fields.
//
// Both meshes MUST have valid cell centroids (i.e. compute_geometry must
// have been called by the caller). The source registry's fields must be
// sized to source.cells().size(); under-sized fields are silently skipped.
//
// Plan : Phase 23, work-item #3 (mapFields cross-mesh field mapper).
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"

#include <cstddef>

namespace simall::solver {

enum class MapMethod {
    Nearest,
    InverseDistance
};

struct MapFieldsOptions {
    MapMethod   method     = MapMethod::Nearest;
    std::size_t kNeighbors = 4;     ///< IDW only.
    double      idwPower   = 2.0;   ///< 1/d^p exponent.
};

struct MapFieldsStats {
    std::size_t scalarFieldsMapped = 0;
    std::size_t vectorFieldsMapped = 0;
    std::size_t sourceCells        = 0;
    std::size_t targetCells        = 0;
    /// Number of target queries that fell back to a linear scan because
    /// the spatial-hash ring search did not find enough candidates.
    std::size_t fallbackQueries    = 0;
};

/// Map every field in `srcFields` from `srcMesh` onto `tgtMesh`,
/// writing results into `dstFields` (sized to the target cell count).
/// Existing entries in `dstFields` for the same field names are overwritten.
MapFieldsStats map_fields(const meshing::Mesh& srcMesh,
                          const FieldRegistry& srcFields,
                          const meshing::Mesh& tgtMesh,
                          FieldRegistry&       dstFields,
                          const MapFieldsOptions& opts = {});

}  // namespace simall::solver
