// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/ProjectNode.hpp
// Phase  : 21 (PROJECT GRAPH)
//
// Node of the project dependency graph. When a node is invalidated, all
// downstream nodes (mesh→solver→results) are marked stale via a topological
// propagation. The GUI workflow tree observes these state changes through
// the EventBus.
// =============================================================================
#pragma once

#include <atomic>
#include <memory>
#include <string>
#include <vector>

namespace simall::io {

enum class NodeKind {
    Project, Geometry, NamedSelection, Material, Mesh, Physics,
    BoundaryCondition, Initialization, Solution, Result
};

enum class NodeState { Empty, Valid, Stale, Computing, Failed };

class ProjectNode {
public:
    ProjectNode(NodeKind kind, std::string name) : kind_(kind), name_(std::move(name)) {}

    NodeKind     kind()  const noexcept { return kind_; }
    const auto&  name()  const noexcept { return name_; }
    NodeState    state() const noexcept { return state_.load(); }

    void add_child(std::unique_ptr<ProjectNode> c) {
        c->parent_ = this; children_.push_back(std::move(c));
    }
    const auto& children() const noexcept { return children_; }
    ProjectNode* parent() const noexcept { return parent_; }

    /// Mark this node and all downstream consumers stale.
    void invalidate() {
        state_.store(NodeState::Stale);
        for (auto& d : downstream_) d->invalidate();
    }

    void add_downstream(ProjectNode* dep) { downstream_.push_back(dep); }
    void set_state(NodeState s) { state_.store(s); }

private:
    NodeKind                                kind_;
    std::string                             name_;
    std::atomic<NodeState>                  state_{NodeState::Empty};
    ProjectNode*                            parent_ = nullptr;
    std::vector<std::unique_ptr<ProjectNode>> children_;
    std::vector<ProjectNode*>               downstream_;
};

}  // namespace simall::io
