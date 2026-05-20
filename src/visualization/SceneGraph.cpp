// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/SceneGraph.cpp
// =============================================================================
#include "visualization/SceneGraph.hpp"

#include <algorithm>
#include <stdexcept>

namespace simall::visualization
{

// ---------------------------------------------------------------------------
//  Transform helpers
// ---------------------------------------------------------------------------
Transform Transform::multiply(const Transform& a, const Transform& b) noexcept
{
    Transform r;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            double s = 0.0;
            for (int k = 0; k < 4; ++k) {
                s += a.m[i * 4 + k] * b.m[k * 4 + j];
            }
            r.m[i * 4 + j] = s;
        }
    }
    return r;
}

util::Vec3d Transform::apply(const util::Vec3d& p) const noexcept
{
    const double x = m[0] * p.x + m[1] * p.y + m[2] * p.z + m[3];
    const double y = m[4] * p.x + m[5] * p.y + m[6] * p.z + m[7];
    const double z = m[8] * p.x + m[9] * p.y + m[10] * p.z + m[11];
    const double w = m[12] * p.x + m[13] * p.y + m[14] * p.z + m[15];
    if (w != 0.0 && w != 1.0) {
        return {x / w, y / w, z / w};
    }
    return {x, y, z};
}

// ---------------------------------------------------------------------------
//  SceneGraph
// ---------------------------------------------------------------------------
SceneGraph::SceneGraph()
{
    Node root;
    root.name = "<root>";
    nodes_.emplace(root_, std::move(root));
}

NodeId SceneGraph::create_node(const std::string& name, NodeId parent)
{
    std::lock_guard lock(mu_);
    NodeId effectiveParent = (parent == kInvalidNodeId) ? root_ : parent;
    if (!nodes_.count(effectiveParent)) {
        throw std::out_of_range("SceneGraph::create_node: parent not found");
    }
    const NodeId id = nextId_++;
    Node n;
    n.name = name;
    n.parent = effectiveParent;
    nodes_.emplace(id, std::move(n));
    nodes_[effectiveParent].children.push_back(id);
    return id;
}

void SceneGraph::destroy_node(NodeId id)
{
    std::lock_guard lock(mu_);
    if (id == root_ || !nodes_.count(id))
        return;
    // Collect descendants iteratively to avoid recursion stack blow-ups on
    // pathological trees.
    std::vector<NodeId> stack{id};
    std::vector<NodeId> kill;
    while (!stack.empty()) {
        NodeId cur = stack.back();
        stack.pop_back();
        kill.push_back(cur);
        auto it = nodes_.find(cur);
        if (it == nodes_.end())
            continue;
        for (NodeId c : it->second.children)
            stack.push_back(c);
    }
    const NodeId parent = nodes_[id].parent;
    auto& siblings = nodes_[parent].children;
    siblings.erase(std::remove(siblings.begin(), siblings.end(), id), siblings.end());
    for (NodeId k : kill)
        nodes_.erase(k);
}

bool SceneGraph::reparent(NodeId child, NodeId newParent)
{
    std::lock_guard lock(mu_);
    if (child == root_ || !nodes_.count(child) || !nodes_.count(newParent)) {
        return false;
    }
    if (descends_from(newParent, child))
        return false; // cycle
    const NodeId oldParent = nodes_[child].parent;
    if (oldParent == newParent)
        return true;
    auto& oldSib = nodes_[oldParent].children;
    oldSib.erase(std::remove(oldSib.begin(), oldSib.end(), child), oldSib.end());
    nodes_[newParent].children.push_back(child);
    nodes_[child].parent = newParent;
    return true;
}

void SceneGraph::set_transform(NodeId id, const Transform& xf)
{
    std::lock_guard lock(mu_);
    if (auto it = nodes_.find(id); it != nodes_.end())
        it->second.local = xf;
}

void SceneGraph::set_visible(NodeId id, bool v)
{
    std::lock_guard lock(mu_);
    if (auto it = nodes_.find(id); it != nodes_.end())
        it->second.visible = v;
}

void SceneGraph::set_local_bbox(NodeId id, const util::BoundingBox& b)
{
    std::lock_guard lock(mu_);
    if (auto it = nodes_.find(id); it != nodes_.end())
        it->second.localBox = b;
}

void SceneGraph::attach_actor(NodeId id, ActorId actor)
{
    std::lock_guard lock(mu_);
    auto it = nodes_.find(id);
    if (it == nodes_.end() || actor == kInvalidActorId)
        return;
    auto& v = it->second.actors;
    if (std::find(v.begin(), v.end(), actor) == v.end())
        v.push_back(actor);
}

void SceneGraph::detach_actor(NodeId id, ActorId actor)
{
    std::lock_guard lock(mu_);
    auto it = nodes_.find(id);
    if (it == nodes_.end())
        return;
    auto& v = it->second.actors;
    v.erase(std::remove(v.begin(), v.end(), actor), v.end());
}

bool SceneGraph::contains(NodeId id) const noexcept
{
    std::lock_guard lock(mu_);
    return nodes_.count(id) != 0;
}

std::string SceneGraph::name(NodeId id) const
{
    std::lock_guard lock(mu_);
    auto it = nodes_.find(id);
    return it == nodes_.end() ? std::string{} : it->second.name;
}

bool SceneGraph::visible(NodeId id) const
{
    std::lock_guard lock(mu_);
    NodeId cur = id;
    while (cur != kInvalidNodeId) {
        auto it = nodes_.find(cur);
        if (it == nodes_.end())
            return false;
        if (!it->second.visible)
            return false;
        if (cur == root_)
            break;
        cur = it->second.parent;
    }
    return true;
}

Transform SceneGraph::world_transform(NodeId id) const
{
    std::lock_guard lock(mu_);
    Transform acc = Transform::identity();
    // Build chain root→…→id, then compose left-to-right.
    std::vector<NodeId> chain;
    NodeId cur = id;
    while (cur != kInvalidNodeId && nodes_.count(cur)) {
        chain.push_back(cur);
        if (cur == root_)
            break;
        cur = nodes_.at(cur).parent;
    }
    for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
        acc = Transform::multiply(acc, nodes_.at(*it).local);
    }
    return acc;
}

util::BoundingBox SceneGraph::world_bbox(NodeId id) const
{
    std::lock_guard lock(mu_);
    if (!nodes_.count(id))
        return {};
    Transform xf = Transform::identity();
    // Recompute from root because we hold the lock and do not want to
    // double-lock world_transform.
    std::vector<NodeId> chain;
    NodeId cur = nodes_.at(id).parent;
    while (cur != kInvalidNodeId && nodes_.count(cur)) {
        chain.push_back(cur);
        if (cur == root_)
            break;
        cur = nodes_.at(cur).parent;
    }
    for (auto it = chain.rbegin(); it != chain.rend(); ++it) {
        xf = Transform::multiply(xf, nodes_.at(*it).local);
    }
    return accumulate_bbox(id, xf);
}

std::vector<NodeId> SceneGraph::children(NodeId id) const
{
    std::lock_guard lock(mu_);
    auto it = nodes_.find(id);
    return it == nodes_.end() ? std::vector<NodeId>{} : it->second.children;
}

std::vector<ActorId> SceneGraph::actors(NodeId id) const
{
    std::lock_guard lock(mu_);
    auto it = nodes_.find(id);
    return it == nodes_.end() ? std::vector<ActorId>{} : it->second.actors;
}

void SceneGraph::traverse(const std::function<bool(NodeId)>& visit) const
{
    std::lock_guard lock(mu_);
    if (!visit)
        return;
    std::vector<NodeId> stack{root_};
    while (!stack.empty()) {
        NodeId cur = stack.back();
        stack.pop_back();
        if (!visit(cur))
            continue;
        auto it = nodes_.find(cur);
        if (it == nodes_.end())
            continue;
        // Push children in reverse so first child is visited first.
        const auto& c = it->second.children;
        for (auto cit = c.rbegin(); cit != c.rend(); ++cit)
            stack.push_back(*cit);
    }
}

std::size_t SceneGraph::node_count() const noexcept
{
    std::lock_guard lock(mu_);
    return nodes_.size();
}

// ---------------------------------------------------------------------------
//  Internals
// ---------------------------------------------------------------------------
bool SceneGraph::descends_from(NodeId candidate, NodeId ancestor) const
{
    NodeId cur = candidate;
    while (cur != kInvalidNodeId) {
        if (cur == ancestor)
            return true;
        if (cur == root_)
            return false;
        auto it = nodes_.find(cur);
        if (it == nodes_.end())
            return false;
        cur = it->second.parent;
    }
    return false;
}

util::BoundingBox SceneGraph::accumulate_bbox(NodeId id, const Transform& parentXf) const
{
    util::BoundingBox out;
    auto it = nodes_.find(id);
    if (it == nodes_.end())
        return out;
    const Transform xf = Transform::multiply(parentXf, it->second.local);
    if (it->second.localBox.valid()) {
        // Expand the 8 corners through the affine.
        const auto& b = it->second.localBox;
        const util::Vec3d corners[8] = {{b.min.x, b.min.y, b.min.z},
                                        {b.max.x, b.min.y, b.min.z},
                                        {b.min.x, b.max.y, b.min.z},
                                        {b.max.x, b.max.y, b.min.z},
                                        {b.min.x, b.min.y, b.max.z},
                                        {b.max.x, b.min.y, b.max.z},
                                        {b.min.x, b.max.y, b.max.z},
                                        {b.max.x, b.max.y, b.max.z}};
        for (auto& c : corners)
            out.expand(xf.apply(c));
    }
    for (NodeId c : it->second.children)
        out.expand(accumulate_bbox(c, xf));
    return out;
}

} // namespace simall::visualization
