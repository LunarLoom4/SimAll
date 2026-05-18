// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/TopologyGraph.hpp
// Phase  : 4 (OPENCASCADE CAD ENGINE) — owns EXACT geometry.
//
// OpenCASCADE owns: TopoDS_Shape, BRep, Geom_Surface, etc.
// VTK NEVER owns exact geometry (Principle 2 of Master Blueprint).
//
// TopologyGraph is the persistent-id graph extracted from a TopoDS_Shape.
// Persistent IDs survive re-tessellation, healing, defeaturing — they are
// the bridge from picking (vtkCellId → TopologyFaceId → solver boundary).
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"

#include <atomic>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::cad {

enum class TopologyType { Vertex, Edge, Wire, Face, Shell, Solid, Compound };

struct TopologyNode {
    util::PersistentId         id;
    TopologyType               type;
    std::vector<util::PersistentId> children;
    std::vector<util::PersistentId> parents;
    util::BoundingBox          aabb;
    void*                      occHandle = nullptr;  // opaque TopoDS_Shape ptr
    std::string                userName;             // for named selections
};

class TopologyGraph {
public:
    util::PersistentId add(TopologyType t, void* occHandle);
    void               link(util::PersistentId parent, util::PersistentId child);

    TopologyNode*       find(util::PersistentId id);
    const TopologyNode* find(util::PersistentId id) const;

    std::vector<util::PersistentId> by_type(TopologyType t) const;

    std::size_t size() const noexcept { return nodes_.size(); }
    const auto& nodes() const noexcept { return nodes_; }
    /// Mutable accessor — for internal use by CAD operations (e.g.
    /// PersistentIdManager) that re-key nodes after a destructive edit.
    auto&       mutable_nodes()    noexcept { return nodes_; }

private:
    std::unordered_map<util::PersistentId, TopologyNode> nodes_;
    std::atomic<util::PersistentId>                      next_id_{1};
};

}  // namespace simall::cad
