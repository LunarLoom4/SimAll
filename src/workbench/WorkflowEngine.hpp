// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/WorkflowEngine.hpp
// Phase  : 22 Pass 22.2
//
// WorkflowEngine -- the *active* layer that sits on top of the (otherwise
// passive) Schematic + StateMachine pair.  It exposes three categories of
// behaviour:
//
//   (1) EDGE INSERTION HOOKS
//       connect()    -- wraps Schematic::add_link() and immediately
//                       triggers downstream invalidation on the target
//                       cell so the freshly-wired sink can never linger
//                       at UpToDate from a previous, unrelated config.
//       disconnect() -- wraps remove_link() and re-derives the target's
//                       state (it may now lack a required input).
//
//   (2) REFRESH POLICY
//       refresh_plan(policy)  -- returns a deterministic, topologically
//                                ordered list of CellIds that need a
//                                refresh under the chosen policy.
//       next_refreshable()    -- the first ReadyOnly cell, or nullopt.
//       refresh_one(id, fn)   -- runs the user-supplied lambda for one
//                                cell, then calls mark_solved /
//                                mark_failed based on its bool return.
//
//   (3) DETERMINISTIC EVALUATION ORDER
//       evaluation_order()    -- the full graph's topological order
//                                with stable min-heap tie-break, so
//                                "refresh project" loops produce the
//                                same output across runs.
//
// The engine owns NEITHER the schematic NOR the state machine -- they
// must outlive the engine.  All methods are single-threaded; the caller
// serializes if needed.
// =============================================================================
#pragma once

#include "workbench/CellAdapter.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/RefreshPolicy.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"

#include <functional>
#include <optional>
#include <vector>

namespace simall::workbench
{

class CellAdapterRegistry;

class WorkflowEngine
{
public:
    WorkflowEngine(Schematic& s, StateMachine& sm) noexcept : schematic_(&s), state_(&sm) {}

    // -- edge insertion -------------------------------------------------
    // Returns whatever Schematic::add_link() returned.  On success, the
    // *target* cell (link.to_cell) is marked modified so any previously
    // UpToDate state cannot survive a new upstream wiring.
    bool connect(const CellLink& link);

    // Mirror of the above for link removal.
    bool disconnect(const CellLink& link);

    // -- refresh policy -------------------------------------------------
    [[nodiscard]] std::vector<CellId> refresh_plan(
        RefreshPolicy policy = RefreshPolicy::ReadyOnly) const;

    [[nodiscard]] std::optional<CellId> next_refreshable(
        RefreshPolicy policy = RefreshPolicy::ReadyOnly) const;

    // Convenience driver.  `runner` is a caller-provided callable
    // (CellId) -> bool, where `true` means "I solved it" and `false`
    // means "I failed".  This is the smallest possible loop body for a
    // "Run cell" button.  Returns whatever the runner returned.
    bool refresh_one(CellId id, std::function<bool(CellId)> runner);

    // Pass 22.5 registry-aware overload.  Looks up the cell's
    // `adapter_id` in the supplied registry, instantiates a fresh
    // adapter, and runs it inside `ctx`.  Outcome dispatch:
    //   * adapter returns true                -> mark_solved, return true
    //   * adapter returns false               -> mark_failed, return false
    //   * adapter_id missing OR not registered -> mark_failed + ctx.error
    //   * cell unknown                        -> returns false without state change
    // The adapter pointer escapes only through the engine; callers do
    // NOT own it.  `last_error()` from the adapter is forwarded into
    // `ctx.error()` so the GUI gets a single canonical place to read.
    bool refresh_one(CellId id, const CellAdapterRegistry& registry, ExecutionContext& ctx);

    // -- deterministic order -------------------------------------------
    [[nodiscard]] std::vector<CellId> evaluation_order() const;

    // -- accessors ------------------------------------------------------
    // Read-only views of the underlying objects.  Useful for command
    // implementations that need to format human-readable labels without
    // re-plumbing a Schematic pointer through their constructors.
    [[nodiscard]] const Schematic& schematic() const noexcept { return *schematic_; }
    [[nodiscard]] const StateMachine& state() const noexcept { return *state_; }

private:
    // A cell is "ready" under ReadyOnly iff its required inputs are wired
    // AND every upstream cell is UpToDate (so the runner can execute
    // immediately without waiting on a parent refresh).
    [[nodiscard]] bool is_ready(CellId) const;

    // A cell has a Failed transitive ancestor.  Used for StopOnFailed.
    [[nodiscard]] bool has_failed_ancestor(CellId) const;

    Schematic* schematic_;
    StateMachine* state_;
};

} // namespace simall::workbench
