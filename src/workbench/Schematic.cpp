// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Schematic.cpp
// Phase  : 22 Pass 22.1
// =============================================================================
#include "workbench/Schematic.hpp"

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>

namespace simall::workbench
{

// ---------------------------------------------------------------------------
// Cells
// ---------------------------------------------------------------------------
CellId Schematic::add_cell(CellKind kind, std::string label)
{
    const CellId id = next_cell_id_++;
    cells_.emplace_back(id, kind, std::move(label));
    return id;
}

bool Schematic::remove_cell(CellId id, std::vector<CellLink>* removed_links_out)
{
    // Resolve cell -- bail if unknown.
    const auto cit =
        std::find_if(cells_.begin(), cells_.end(), [id](const Cell& c) { return c.id() == id; });
    if (cit == cells_.end())
        return false;

    // Snapshot + erase every incident link (incoming or outgoing).
    if (removed_links_out)
        removed_links_out->clear();
    for (auto it = links_.begin(); it != links_.end();) {
        if (it->from_cell == id || it->to_cell == id) {
            if (removed_links_out)
                removed_links_out->push_back(*it);
            it = links_.erase(it);
        } else {
            ++it;
        }
    }
    cells_.erase(cit);
    return true;
}

bool Schematic::restore_cell(Cell c)
{
    // Refuse id collisions -- undo would otherwise silently shadow an
    // unrelated cell minted in the meantime.
    if (cell(c.id()) != nullptr)
        return false;
    const CellId restored_id = c.id();
    cells_.push_back(std::move(c));
    if (restored_id >= next_cell_id_)
        next_cell_id_ = restored_id + 1;
    return true;
}

Cell* Schematic::cell(CellId id) noexcept
{
    for (auto& c : cells_) {
        if (c.id() == id)
            return &c;
    }
    return nullptr;
}

const Cell* Schematic::cell(CellId id) const noexcept
{
    for (const auto& c : cells_) {
        if (c.id() == id)
            return &c;
    }
    return nullptr;
}


// ---------------------------------------------------------------------------
// Cycle detection -- DFS forward from `to` and check whether we can reach
// `from`.  If yes, adding from -> to would close a cycle.
// ---------------------------------------------------------------------------
bool Schematic::would_create_cycle(CellId from, CellId to) const
{
    if (from == to)
        return true; // self-loop

    std::unordered_set<CellId> visited;
    std::vector<CellId> stack{to};

    while (!stack.empty()) {
        const CellId cur = stack.back();
        stack.pop_back();
        if (!visited.insert(cur).second)
            continue;

        for (const auto& l : links_) {
            if (l.from_cell == cur) {
                if (l.to_cell == from)
                    return true;
                stack.push_back(l.to_cell);
            }
        }
    }
    return false;
}


// ---------------------------------------------------------------------------
// Links
// ---------------------------------------------------------------------------
bool Schematic::add_link(const CellLink& link)
{
    const Cell* from = cell(link.from_cell);
    const Cell* to = cell(link.to_cell);
    if (!from || !to)
        return false;

    const CellPort* fp = from->find_port(link.from_port);
    const CellPort* tp = to->find_port(link.to_port);
    if (!fp || !tp)
        return false;

    if (fp->direction != PortDirection::Output)
        return false;
    if (tp->direction != PortDirection::Input)
        return false;

    if (!fp->data_type.empty() && !tp->data_type.empty() && fp->data_type != tp->data_type) {
        return false;
    }

    // Duplicate edge?
    if (std::find(links_.begin(), links_.end(), link) != links_.end()) {
        return false;
    }

    // Single fan-in on input ports.
    for (const auto& l : links_) {
        if (l.to_cell == link.to_cell && l.to_port == link.to_port) {
            return false;
        }
    }

    if (would_create_cycle(link.from_cell, link.to_cell)) {
        return false;
    }

    links_.push_back(link);
    return true;
}

bool Schematic::remove_link(const CellLink& link)
{
    const auto it = std::find(links_.begin(), links_.end(), link);
    if (it == links_.end())
        return false;
    links_.erase(it);
    return true;
}


// ---------------------------------------------------------------------------
// Graph queries
// ---------------------------------------------------------------------------
std::vector<CellId> Schematic::upstream(CellId id) const
{
    std::vector<CellId> out;
    for (const auto& l : links_) {
        if (l.to_cell == id)
            out.push_back(l.from_cell);
    }
    return out;
}

std::vector<CellId> Schematic::downstream(CellId id) const
{
    std::vector<CellId> out;
    for (const auto& l : links_) {
        if (l.from_cell == id)
            out.push_back(l.to_cell);
    }
    return out;
}

bool Schematic::inputs_satisfied(CellId id) const
{
    const Cell* c = cell(id);
    if (!c)
        return false;

    for (const auto& p : c->ports()) {
        if (p.direction != PortDirection::Input)
            continue;
        if (!p.required)
            continue;

        // Look for at least one link feeding this input.
        bool wired = false;
        for (const auto& l : links_) {
            if (l.to_cell == id && l.to_port == p.id) {
                wired = true;
                break;
            }
        }
        if (!wired)
            return false;
    }
    return true;
}


// ---------------------------------------------------------------------------
// Kahn's topological sort.  Ties broken by CellId ascending for a stable,
// deterministic output -- the tests depend on it.
// ---------------------------------------------------------------------------
std::vector<CellId> Schematic::topological_order() const
{
    std::unordered_map<CellId, std::size_t> indegree;
    indegree.reserve(cells_.size());
    for (const auto& c : cells_)
        indegree[c.id()] = 0;
    for (const auto& l : links_)
        ++indegree[l.to_cell];

    auto cmp = [](CellId a, CellId b) { return a > b; }; // min-heap
    std::priority_queue<CellId, std::vector<CellId>, decltype(cmp)> ready(cmp);
    for (const auto& [id, n] : indegree) {
        if (n == 0)
            ready.push(id);
    }

    std::vector<CellId> order;
    order.reserve(cells_.size());
    while (!ready.empty()) {
        const CellId cur = ready.top();
        ready.pop();
        order.push_back(cur);
        for (const auto& l : links_) {
            if (l.from_cell == cur) {
                if (--indegree[l.to_cell] == 0)
                    ready.push(l.to_cell);
            }
        }
    }

    if (order.size() != cells_.size())
        return {}; // cycle (unreachable)
    return order;
}

} // namespace simall::workbench
