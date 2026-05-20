// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/WorkflowEngine.cpp
// Phase  : 22 Pass 22.2
// =============================================================================
#include "workbench/WorkflowEngine.hpp"

#include "workbench/CellAdapterRegistry.hpp"

#include <string>
#include <unordered_set>

namespace simall::workbench {

// ---------------------------------------------------------------------------
// Edge insertion / removal
// ---------------------------------------------------------------------------
bool WorkflowEngine::connect(const CellLink& link) {
    if (!schematic_->add_link(link)) return false;
    // Downstream-invalidation rule: a brand-new upstream edge always
    // invalidates the sink (and everyone downstream of it).  We can NOT
    // simply call recompute_all here -- doing so would silently demote
    // Failed cells if their upstream became OK, and the user expects
    // Failed to be sticky.  mark_modified is the correct primitive.
    state_->mark_modified(link.to_cell);
    return true;
}

bool WorkflowEngine::disconnect(const CellLink& link) {
    if (!schematic_->remove_link(link)) return false;
    // Removing an edge can leave the sink with an unwired required
    // input.  mark_modified will re-derive correctly (Unfulfilled when
    // a required input is missing).
    state_->mark_modified(link.to_cell);
    return true;
}


// ---------------------------------------------------------------------------
// Readiness predicates
// ---------------------------------------------------------------------------
bool WorkflowEngine::is_ready(CellId id) const {
    if (!schematic_->inputs_satisfied(id)) return false;
    for (const CellId u : schematic_->upstream(id)) {
        const Cell* parent = schematic_->cell(u);
        if (!parent) return false;
        if (parent->state() != CellState::UpToDate) return false;
    }
    return true;
}

bool WorkflowEngine::has_failed_ancestor(CellId id) const {
    std::unordered_set<CellId> seen;
    std::vector<CellId>        stack{id};
    while (!stack.empty()) {
        const CellId cur = stack.back();
        stack.pop_back();
        if (!seen.insert(cur).second) continue;
        for (const CellId u : schematic_->upstream(cur)) {
            const Cell* p = schematic_->cell(u);
            if (!p) continue;
            if (p->state() == CellState::Failed) return true;
            stack.push_back(u);
        }
    }
    return false;
}


// ---------------------------------------------------------------------------
// Refresh policy
// ---------------------------------------------------------------------------
std::vector<CellId>
WorkflowEngine::refresh_plan(RefreshPolicy policy) const {
    const auto order = schematic_->topological_order();
    std::vector<CellId> plan;
    plan.reserve(order.size());

    for (const CellId id : order) {
        const Cell* c = schematic_->cell(id);
        if (!c) continue;
        if (c->state() != CellState::RefreshRequired) continue;

        if (has_flag(policy, RefreshPolicy::StopOnFailed)
            && has_failed_ancestor(id)) {
            continue;
        }

        const bool include_stale = has_flag(policy, RefreshPolicy::IncludeStale);
        if (!include_stale && !is_ready(id)) {
            continue;
        }
        plan.push_back(id);
    }
    return plan;
}


std::optional<CellId>
WorkflowEngine::next_refreshable(RefreshPolicy policy) const {
    const auto plan = refresh_plan(policy);
    if (plan.empty()) return std::nullopt;
    return plan.front();
}


bool WorkflowEngine::refresh_one(CellId id,
                                 std::function<bool(CellId)> runner) {
    const bool ok = runner ? runner(id) : false;
    if (ok) state_->mark_solved(id);
    else    state_->mark_failed(id);
    return ok;
}


// ---------------------------------------------------------------------------
// Registry-backed overload -- the canonical "Run cell" entry point.
// ---------------------------------------------------------------------------
bool WorkflowEngine::refresh_one(CellId id,
                                 const CellAdapterRegistry& registry,
                                 ExecutionContext&          ctx) {
    Cell* cell = schematic_->cell(id);
    if (!cell) {
        ctx.error("workflow: unknown cell id");
        return false;
    }

    const std::string& aid = cell->adapter_id();
    if (aid.empty()) {
        ctx.error("workflow: cell has no adapter_id; cannot refresh");
        state_->mark_failed(id);
        return false;
    }

    auto adapter = registry.make(aid);
    if (!adapter) {
        ctx.error(std::string("workflow: no adapter registered for id '")
                  + aid + "'");
        state_->mark_failed(id);
        return false;
    }

    // Sanity-check the kind match.  We only WARN -- a "Custom" adapter
    // is intentionally polymorphic and may legitimately serve any kind.
    if (adapter->kind() != cell->kind()
        && adapter->kind() != CellKind::Custom) {
        ctx.warn(std::string("workflow: adapter kind mismatch for '")
                 + aid + "'");
    }

    const bool ok = adapter->execute(*cell, ctx);
    if (!ok) {
        const auto err = adapter->last_error();
        if (!err.empty()) ctx.error(err);
        state_->mark_failed(id);
        return false;
    }
    state_->mark_solved(id);
    return true;
}


// ---------------------------------------------------------------------------
// Deterministic full-graph order
// ---------------------------------------------------------------------------
std::vector<CellId> WorkflowEngine::evaluation_order() const {
    return schematic_->topological_order();
}

}  // namespace simall::workbench
