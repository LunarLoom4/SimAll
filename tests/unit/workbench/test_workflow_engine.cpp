// =============================================================================
// SimAll Beta — Workbench Unit Tests
// File   : tests/unit/workbench/test_workflow_engine.cpp
// Phase  : 22 Pass 22.2
//
// Validates the WorkflowEngine layer on top of Schematic + StateMachine:
//   * connect() / disconnect() invalidate downstream on every edit.
//   * refresh_plan() respects ReadyOnly / IncludeStale / StopOnFailed
//     policies and returns a deterministic topological order.
//   * next_refreshable() picks the first plan entry.
//   * refresh_one() dispatches to mark_solved / mark_failed correctly.
//   * evaluation_order() is stable across runs.
// =============================================================================
#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/RefreshPolicy.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"
#include "workbench/WorkflowEngine.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <optional>
#include <vector>

using namespace simall::workbench;

namespace
{

// A 5-cell pipeline: G -> M -> Su -> So -> R, all with "x" data type.
struct Pipe
{
    Schematic s;
    StateMachine sm{s};
    WorkflowEngine eng{s, sm};
    CellId g{}, m{}, su{}, so{}, r{};
    CellId build()
    {
        g = s.add_cell(CellKind::Geometry, "G");
        m = s.add_cell(CellKind::Mesh, "M");
        su = s.add_cell(CellKind::Setup, "Su");
        so = s.add_cell(CellKind::Solution, "So");
        r = s.add_cell(CellKind::Results, "R");

        const PortId g_o = s.cell(g)->add_output("x", "x");
        const PortId m_i = s.cell(m)->add_input("x", "x");
        const PortId m_o = s.cell(m)->add_output("x", "x");
        const PortId su_i = s.cell(su)->add_input("x", "x");
        const PortId su_o = s.cell(su)->add_output("x", "x");
        const PortId so_i = s.cell(so)->add_input("x", "x");
        const PortId so_o = s.cell(so)->add_output("x", "x");
        const PortId r_i = s.cell(r)->add_input("x", "x");

        REQUIRE(eng.connect({g, g_o, m, m_i}));
        REQUIRE(eng.connect({m, m_o, su, su_i}));
        REQUIRE(eng.connect({su, su_o, so, so_i}));
        REQUIRE(eng.connect({so, so_o, r, r_i}));
        return r;
    }
};

} // namespace


// ---------------------------------------------------------------------------
// RefreshPolicy enum
// ---------------------------------------------------------------------------
TEST_CASE("RefreshPolicy stringifies and supports flag composition", "[workbench][refresh-policy]")
{
    REQUIRE(to_string(RefreshPolicy::ReadyOnly) == "ReadyOnly");
    REQUIRE(to_string(RefreshPolicy::IncludeStale) == "IncludeStale");
    REQUIRE(to_string(RefreshPolicy::StopOnFailed) == "StopOnFailed");
    REQUIRE(to_string(RefreshPolicy::UpdateProject) == "UpdateProject");

    REQUIRE(has_flag(RefreshPolicy::UpdateProject, RefreshPolicy::IncludeStale));
    REQUIRE(has_flag(RefreshPolicy::UpdateProject, RefreshPolicy::StopOnFailed));
    REQUIRE_FALSE(has_flag(RefreshPolicy::ReadyOnly, RefreshPolicy::IncludeStale));
    REQUIRE_FALSE(has_flag(RefreshPolicy::IncludeStale, RefreshPolicy::StopOnFailed));
}


// ---------------------------------------------------------------------------
// connect() / disconnect() side effects
// ---------------------------------------------------------------------------
TEST_CASE("connect() invalidates the target cell beyond add_link", "[workbench][engine][connect]")
{
    Pipe p;
    p.build();
    // Solve everything so the downstream is in the most-OK state.
    for (const CellId id : {p.g, p.m, p.su, p.so, p.r})
        p.sm.mark_solved(id);
    REQUIRE(p.s.cell(p.r)->state() == CellState::UpToDate);

    // Add a fresh sink to mesh and connect it.  m is currently UpToDate;
    // a new wiring on m must demote m to RefreshRequired *and* propagate.
    const CellId t = p.s.add_cell(CellKind::Custom, "T");
    const PortId to = p.s.cell(p.m)->add_output("y", "y");
    const PortId ti = p.s.cell(t)->add_input("y", "y");
    REQUIRE(p.eng.connect({p.m, to, t, ti}));

    REQUIRE(p.s.cell(t)->state() == CellState::RefreshRequired);
    // m itself stays UpToDate -- adding an *outgoing* edge doesn't dirty
    // the source.  But the new sink (t) is invalidated as required.
    REQUIRE(p.s.cell(p.m)->state() == CellState::UpToDate);
}


TEST_CASE("disconnect() reverts the sink to Unfulfilled when its required input vanishes",
          "[workbench][engine][disconnect]")
{
    Pipe p;
    p.build();
    for (const CellId id : {p.g, p.m, p.su, p.so, p.r})
        p.sm.mark_solved(id);

    // Drop g -> m.
    const PortId g_o = p.s.cell(p.g)->ports()[0].id;
    const PortId m_i = p.s.cell(p.m)->ports()[0].id;
    REQUIRE(p.eng.disconnect({p.g, g_o, p.m, m_i}));

    REQUIRE(p.s.cell(p.m)->state() == CellState::Unfulfilled);
    // su loses its only ready upstream chain -> also Unfulfilled.
    REQUIRE(p.s.cell(p.su)->state() == CellState::Unfulfilled);
    REQUIRE(p.s.cell(p.r)->state() == CellState::Unfulfilled);
}


// ---------------------------------------------------------------------------
// refresh_plan() with ReadyOnly
// ---------------------------------------------------------------------------
TEST_CASE("refresh_plan(ReadyOnly) returns only cells with all-UpToDate parents",
          "[workbench][engine][refresh-plan]")
{
    Pipe p;
    p.build();
    p.sm.recompute_all(); // every cell -> RefreshRequired

    // Only g has no parents, so only g is "ready".
    auto plan = p.eng.refresh_plan(RefreshPolicy::ReadyOnly);
    REQUIRE(plan == std::vector<CellId>{p.g});

    REQUIRE(p.eng.next_refreshable() == std::optional<CellId>{p.g});

    p.sm.mark_solved(p.g); // g done
    plan = p.eng.refresh_plan(RefreshPolicy::ReadyOnly);
    REQUIRE(plan == std::vector<CellId>{p.m});

    p.sm.mark_solved(p.m);
    p.sm.mark_solved(p.su);
    p.sm.mark_solved(p.so);
    p.sm.mark_solved(p.r);
    REQUIRE(p.eng.refresh_plan(RefreshPolicy::ReadyOnly).empty());
    REQUIRE(p.eng.next_refreshable() == std::nullopt);
}


// ---------------------------------------------------------------------------
// refresh_plan() with IncludeStale -- full topo chain.
// ---------------------------------------------------------------------------
TEST_CASE("refresh_plan(IncludeStale) returns every RefreshRequired cell in topo order",
          "[workbench][engine][refresh-plan]")
{
    Pipe p;
    p.build();
    p.sm.recompute_all(); // all RefreshRequired

    const auto plan = p.eng.refresh_plan(RefreshPolicy::IncludeStale);
    REQUIRE(plan == std::vector<CellId>{p.g, p.m, p.su, p.so, p.r});
}


// ---------------------------------------------------------------------------
// refresh_plan() with StopOnFailed -- prune cells under a Failed ancestor.
// ---------------------------------------------------------------------------
TEST_CASE("refresh_plan(StopOnFailed) prunes cells whose ancestor failed",
          "[workbench][engine][refresh-plan][failed]")
{
    Pipe p;
    p.build();
    for (const CellId id : {p.g, p.m, p.su, p.so, p.r})
        p.sm.mark_solved(id);
    p.sm.mark_failed(p.su); // -> su=Failed, so & r = Unfulfilled

    // Now mark g modified so it shows up as RefreshRequired.
    p.sm.mark_modified(p.g);
    // g, m are RefreshRequired (not under su).
    // so, r are downstream of failed su -> Unfulfilled, not in plan anyway.

    const auto plan_with = p.eng.refresh_plan(RefreshPolicy::UpdateProject);
    const auto plan_without = p.eng.refresh_plan(RefreshPolicy::IncludeStale);

    // Without StopOnFailed: g and m show up.  With StopOnFailed: m would
    // also show up unless m had a failed ancestor (it doesn't -- only su
    // is Failed, and su is *downstream* of m).  So both plans are equal.
    REQUIRE(plan_with == std::vector<CellId>{p.g, p.m});
    REQUIRE(plan_without == std::vector<CellId>{p.g, p.m});

    // Now mark g failed.  Under StopOnFailed, m must drop out (its
    // ancestor is failed).  Under plain IncludeStale, m survives.
    p.sm.mark_failed(p.g);
    p.sm.mark_modified(p.m); // m is now RefreshRequired (own state asserted)
    REQUIRE(p.s.cell(p.m)->state() == CellState::RefreshRequired);

    const auto plan_stop = p.eng.refresh_plan(RefreshPolicy::UpdateProject);
    const auto plan_open = p.eng.refresh_plan(RefreshPolicy::IncludeStale);
    REQUIRE(plan_stop.empty()); // m pruned
    REQUIRE(plan_open == std::vector<CellId>{p.m});
}


// ---------------------------------------------------------------------------
// refresh_one() driver
// ---------------------------------------------------------------------------
TEST_CASE("refresh_one() marks Solved on true and Failed on false",
          "[workbench][engine][refresh-one]")
{
    Pipe p;
    p.build();
    p.sm.recompute_all();
    REQUIRE(p.s.cell(p.g)->state() == CellState::RefreshRequired);

    REQUIRE(p.eng.refresh_one(p.g, [](CellId) { return true; }));
    REQUIRE(p.s.cell(p.g)->state() == CellState::UpToDate);

    REQUIRE_FALSE(p.eng.refresh_one(p.m, [](CellId) { return false; }));
    REQUIRE(p.s.cell(p.m)->state() == CellState::Failed);

    // Downstream of failed m is Unfulfilled.
    REQUIRE(p.s.cell(p.su)->state() == CellState::Unfulfilled);
}


// ---------------------------------------------------------------------------
// evaluation_order() is identical across calls (deterministic).
// ---------------------------------------------------------------------------
TEST_CASE("evaluation_order() is deterministic", "[workbench][engine][order]")
{
    Pipe p;
    p.build();
    const auto a = p.eng.evaluation_order();
    const auto b = p.eng.evaluation_order();
    REQUIRE(a == b);
    REQUIRE(a == std::vector<CellId>{p.g, p.m, p.su, p.so, p.r});
}
