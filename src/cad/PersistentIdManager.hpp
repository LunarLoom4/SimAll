// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/PersistentIdManager.hpp
// Phase  : 4.5 (Persistent IDs across geometry modifications)
//
// Geometry edits (heal, defeature, boolean) destroy OCC TopoDS_Shape* and
// replace them with new ones. SimAll's persistent IDs MUST survive these
// transformations so that named selections, BCs, and result probes stay
// glued to the correct topology.
//
// Strategy: before an edit, snapshot a fingerprint per entity (id, centroid,
// area/length, bounding box, type). After the edit, walk the new topology
// and match by spatial proximity + type + size; transfer the old id when a
// unique match within tolerance is found. Unmatched entities get fresh ids.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"
#include "cad/TopologyGraph.hpp"
#include "utilities/MathTypes.hpp"

#include <unordered_map>
#include <vector>

namespace simall::cad
{

struct PersistentFingerprint
{
    util::PersistentId id;
    TopologyType type;
    util::Vec3d centroid;
    util::BoundingBox bbox;
    double measure = 0.0; // length for edges, area for faces, volume for solids
};

struct RemapResult
{
    std::unordered_map<util::PersistentId, util::PersistentId> oldToNew;
    std::vector<util::PersistentId> unmatchedOld;
    std::vector<util::PersistentId> freshlyAllocated;
};

class PersistentIdManager
{
public:
    /// Capture a fingerprint of every face / edge / vertex / solid in shape.
    /// Call BEFORE the edit. Replaces any previously stored snapshot.
    void snapshot(const ShapeHandle& shape);

    /// After the edit, attempt to transfer ids from the snapshot onto
    /// `newShape`. Returns the mapping and a list of orphaned old ids.
    RemapResult remapTo(ShapeHandle& newShape,
                        double positionTol = 1e-4,
                        double measureRelTol = 0.05);

    std::size_t snapshotSize() const noexcept { return prints_.size(); }

private:
    std::vector<PersistentFingerprint> prints_;
};

} // namespace simall::cad
