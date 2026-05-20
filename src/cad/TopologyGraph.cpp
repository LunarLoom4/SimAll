#include "cad/TopologyGraph.hpp"

namespace simall::cad
{

util::PersistentId TopologyGraph::add(TopologyType t, void* occ)
{
    auto id = next_id_.fetch_add(1, std::memory_order_relaxed);
    nodes_[id] = TopologyNode{id, t, {}, {}, {}, occ, {}};
    return id;
}

void TopologyGraph::link(util::PersistentId parent, util::PersistentId child)
{
    auto* p = find(parent);
    auto* c = find(child);
    if (!p || !c)
        return;
    p->children.push_back(child);
    c->parents.push_back(parent);
}

TopologyNode* TopologyGraph::find(util::PersistentId id)
{
    auto it = nodes_.find(id);
    return it == nodes_.end() ? nullptr : &it->second;
}
const TopologyNode* TopologyGraph::find(util::PersistentId id) const
{
    auto it = nodes_.find(id);
    return it == nodes_.end() ? nullptr : &it->second;
}

std::vector<util::PersistentId> TopologyGraph::by_type(TopologyType t) const
{
    std::vector<util::PersistentId> out;
    for (const auto& [id, n] : nodes_)
        if (n.type == t)
            out.push_back(id);
    return out;
}

} // namespace simall::cad
