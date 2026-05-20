// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/SceneGraph.hpp
// Phase  : 3 / Week 14 (scene-management trio: SceneGraph)
//
// Hierarchical node tree with affine transforms, visibility flags, and lazy
// bounding-box propagation.  Nodes can carry zero or more ActorIds — the
// SceneGraph is *organisational*, the heavy actor lifetime is owned by the
// ActorRegistry.  Independent of VTK so it can be flattened to any
// back-end (VTK, OSPRay, OpenGL ES, software rasteriser).
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <array>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::visualization
{

/// Affine 4x4 in row-major layout. Identity by default.
struct Transform
{
    std::array<double, 16> m{1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

    static Transform identity() noexcept { return Transform{}; }
    static Transform translation(const util::Vec3d& t) noexcept
    {
        Transform tr;
        tr.m[3] = t.x;
        tr.m[7] = t.y;
        tr.m[11] = t.z;
        return tr;
    }
    static Transform scale(double s) noexcept
    {
        Transform tr;
        tr.m[0] = tr.m[5] = tr.m[10] = s;
        return tr;
    }
    static Transform multiply(const Transform& a, const Transform& b) noexcept;

    util::Vec3d apply(const util::Vec3d& p) const noexcept;
};

class SceneGraph
{
public:
    SceneGraph();

    /// Create a child node.  Pass `kInvalidNodeId` (the default) to attach
    /// directly under the root.  Throws std::out_of_range if `parent` does
    /// not exist.
    NodeId create_node(const std::string& name, NodeId parent = kInvalidNodeId);

    /// Detach a node and all its descendants.  Silently no-ops on root or
    /// on missing ids — chainable from controllers.
    void destroy_node(NodeId id);

    /// Reparent.  Refuses cycles (returns false); succeeds otherwise.
    bool reparent(NodeId child, NodeId newParent);

    void set_transform(NodeId id, const Transform& xform);
    void set_visible(NodeId id, bool visible);
    void set_local_bbox(NodeId id, const util::BoundingBox& box);

    /// Attach/remove actors on a node.  Both are O(1) average.
    void attach_actor(NodeId id, ActorId actor);
    void detach_actor(NodeId id, ActorId actor);

    NodeId root() const noexcept { return root_; }
    bool contains(NodeId id) const noexcept;
    std::string name(NodeId id) const;
    bool visible(NodeId id) const;                 // effective (parents AND)
    Transform world_transform(NodeId id) const;    // composed root→leaf
    util::BoundingBox world_bbox(NodeId id) const; // recursive union
    std::vector<NodeId> children(NodeId id) const;
    std::vector<ActorId> actors(NodeId id) const;

    /// DFS visit (pre-order) of every node, useful for renderers.  visit(id)
    /// returning false skips that node's subtree.
    void traverse(const std::function<bool(NodeId)>& visit) const;

    std::size_t node_count() const noexcept;

private:
    struct Node
    {
        std::string name;
        NodeId parent = kInvalidNodeId;
        std::vector<NodeId> children;
        std::vector<ActorId> actors;
        Transform local;
        util::BoundingBox localBox; // user-supplied
        bool visible = true;
    };

    bool descends_from(NodeId candidate, NodeId ancestor) const;
    util::BoundingBox accumulate_bbox(NodeId id, const Transform& parentXf) const;

    mutable std::mutex mu_;
    std::unordered_map<NodeId, Node> nodes_;
    NodeId root_ = 1;
    NodeId nextId_ = 2;
};

} // namespace simall::visualization
