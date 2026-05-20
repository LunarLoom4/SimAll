// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/StateMachine.cpp
// Phase  : 22 Pass 22.1
//
// Implementation note on the propagation rules:
//
//   derive_state(c) walks c's *direct* upstream cells and combines their
//   states with c's wiring to produce the natural settled state for c.
//   The numeric ordering of CellState (UpToDate=0 < RefreshRequired=1 <
//   Unfulfilled=2 < Failed=3) lets a single std::max() express the rule
//   "the worst upstream state dominates".
//
//   * Unfulfilled inputs       -> Unfulfilled (hard).
//   * Failed upstream parent   -> Unfulfilled (we cannot derive output
//                                from broken input -- but the cell itself
//                                is not Failed, the parent is).
//   * RefreshRequired upstream -> RefreshRequired (downstream is stale).
//   * Otherwise                -> preserve the cell's own state, which is
//                                set by the user via mark_solved /
//                                mark_failed / mark_modified.
//
//   mark_solved / mark_failed / mark_modified set the cell's own state,
//   then sweep the transitive downstream subgraph in topological order so
//   each child sees parents that have already settled this pass.
// =============================================================================
#include "workbench/StateMachine.hpp"

#include <queue>
#include <unordered_set>

namespace simall::workbench
{

CellState StateMachine::derive_state(CellId id) const
{
    if (!schematic_->inputs_satisfied(id)) {
        return CellState::Unfulfilled;
    }

    CellState worst = CellState::UpToDate;
    for (const CellId u : schematic_->upstream(id)) {
        const Cell* parent = schematic_->cell(u);
        if (!parent)
            continue;

        switch (parent->state()) {
        case CellState::Failed:
            // Broken producer -> child can't compute -> Unfulfilled.
            worst = std::max(worst, CellState::Unfulfilled);
            break;
        case CellState::Unfulfilled:
            worst = std::max(worst, CellState::Unfulfilled);
            break;
        case CellState::RefreshRequired:
            worst = std::max(worst, CellState::RefreshRequired);
            break;
        case CellState::UpToDate:
            break;
        }
    }
    return worst;
}


void StateMachine::recompute_all()
{
    // Pin user-asserted states (UpToDate / Failed) and re-derive the rest.
    const auto order = schematic_->topological_order();
    for (const CellId id : order) {
        Cell* c = schematic_->cell(id);
        if (!c)
            continue;

        // Failed is a sticky terminal user signal -- only mark_modified /
        // mark_solved can clear it.  Everything else (Unfulfilled /
        // RefreshRequired / UpToDate) gets re-derived from the graph.
        if (c->state() == CellState::Failed)
            continue;

        const CellState derived = derive_state(id);

        // If the cell was previously UpToDate and the derived state is
        // still ok, keep UpToDate -- our parents haven't shifted.
        if (c->state() == CellState::UpToDate && derived == CellState::UpToDate) {
            continue;
        }
        c->set_state(derived == CellState::UpToDate
                         ? CellState::RefreshRequired // ready but not yet solved
                         : derived);
    }
}


namespace
{

// BFS over the downstream subgraph from `seed` (exclusive); invokes `fn`
// on every transitive child in some topological-friendly order.
template <typename Fn> void visit_downstream(const Schematic& s, CellId seed, Fn fn)
{
    std::unordered_set<CellId> seen;
    std::queue<CellId> q;
    for (const CellId d : s.downstream(seed))
        q.push(d);

    while (!q.empty()) {
        const CellId cur = q.front();
        q.pop();
        if (!seen.insert(cur).second)
            continue;
        fn(cur);
        for (const CellId d : s.downstream(cur))
            q.push(d);
    }
}

} // namespace


void StateMachine::mark_modified(CellId id)
{
    Cell* c = schematic_->cell(id);
    if (!c)
        return;

    const CellState own =
        schematic_->inputs_satisfied(id) ? CellState::RefreshRequired : CellState::Unfulfilled;
    c->set_state(own);

    visit_downstream(*schematic_, id, [&](CellId child) {
        Cell* k = schematic_->cell(child);
        if (!k)
            return;
        if (k->state() == CellState::Failed)
            return; // sticky
        k->set_state(derive_state(child));
        if (k->state() == CellState::UpToDate) {
            // Parent moved; child must re-run.
            k->set_state(CellState::RefreshRequired);
        }
    });
}


void StateMachine::mark_solved(CellId id)
{
    Cell* c = schematic_->cell(id);
    if (!c)
        return;

    c->set_state(CellState::UpToDate);

    visit_downstream(*schematic_, id, [&](CellId child) {
        Cell* k = schematic_->cell(child);
        if (!k)
            return;
        if (k->state() == CellState::Failed)
            return; // sticky
        const CellState derived = derive_state(child);
        // If the child was UpToDate, the freshly solved parent
        // invalidates it -> RefreshRequired.  Otherwise honour the
        // derived state.
        if (k->state() == CellState::UpToDate && derived == CellState::UpToDate) {
            k->set_state(CellState::RefreshRequired);
        } else {
            k->set_state(derived);
        }
    });
}


void StateMachine::mark_failed(CellId id)
{
    Cell* c = schematic_->cell(id);
    if (!c)
        return;

    c->set_state(CellState::Failed);

    visit_downstream(*schematic_, id, [&](CellId child) {
        Cell* k = schematic_->cell(child);
        if (!k)
            return;
        if (k->state() == CellState::Failed)
            return; // sticky
        // A failed ancestor poisons everyone downstream -> Unfulfilled.
        k->set_state(CellState::Unfulfilled);
    });
}

} // namespace simall::workbench
