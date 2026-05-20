// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/ProjectGraph.hpp
// Week   : 17
//
// Project dependency invalidation graph.  Extends `ProjectNode` (a tree of
// workflow steps) with:
//
//   * Hash-based change detection — every node carries a 64-bit content
//     hash; mutating commit() recomputes it and propagates Stale to all
//     transitive downstream consumers.
//
//   * Topological ordering — `topological_order()` returns nodes in safe
//     compute order; used by the solver pipeline driver to know which
//     stages can run in parallel.
//
//   * Observer pattern — callers subscribe to per-node state-change events
//     (used by the GUI workflow tree to refresh its icons).
//
//   * Deterministic IDs — each node gets a stable `uuid` so the graph
//     survives `.simall` save/load round-trips.
// =============================================================================
#pragma once

#include "io/ProjectNode.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::io
{

[[nodiscard]] std::uint64_t hash_bytes(const void* data,
                                       std::size_t n,
                                       std::uint64_t seed = 1469598103934665603ULL) noexcept;

class ProjectGraph
{
public:
    ProjectGraph();

    // -- node operations -----------------------------------------------------
    ProjectNode* add_node(NodeKind kind, std::string name, std::string uuid = {});
    ProjectNode* find_by_uuid(std::string_view uuid) const noexcept;
    [[nodiscard]] const std::string& uuid_of(const ProjectNode* n) const;

    void link(ProjectNode* from, ProjectNode* to);
    void commit(ProjectNode* n, const void* contentBytes, std::size_t n_bytes);

    // -- queries -------------------------------------------------------------
    [[nodiscard]] std::vector<ProjectNode*> topological_order() const;
    [[nodiscard]] std::vector<ProjectNode*> stale_nodes() const;
    [[nodiscard]] std::uint64_t hash_of(const ProjectNode* n) const;

    // -- observer ------------------------------------------------------------
    using StateObserver = std::function<void(ProjectNode*, NodeState old, NodeState now)>;
    int subscribe(StateObserver obs);
    void unsubscribe(int id);

private:
    struct NodeMeta
    {
        std::string uuid;
        std::uint64_t contentHash = 0;
    };
    std::unordered_map<const ProjectNode*, NodeMeta> meta_;
    std::vector<std::unique_ptr<ProjectNode>> owned_;
    std::unordered_map<int, StateObserver> observers_;
    int nextObsId_ = 1;
};

} // namespace simall::io
