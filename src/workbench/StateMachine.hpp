// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/StateMachine.hpp
// Phase  : 22 Pass 22.1
//
// StateMachine -- deterministic propagation of CellState across a Schematic
// in response to user signals (mark_modified / mark_solved / mark_failed)
// and to changes in connectivity (recompute_all).
//
// Single owner: the Schematic must outlive this object.  No threading; the
// caller serializes access if needed.
// =============================================================================
#pragma once

#include "workbench/Schematic.hpp"
#include "workbench/Workbench.hpp"

namespace simall::workbench {

class StateMachine {
public:
    explicit StateMachine(Schematic& s) noexcept : schematic_(&s) {}

    // Recompute every cell's state from scratch using the current graph
    // and the user signal already encoded in each cell (UpToDate /
    // Failed are preserved; Unfulfilled / RefreshRequired are derived).
    // Visits cells in topological order so children see settled parents.
    void recompute_all();

    // The cell's data was edited externally.  It becomes RefreshRequired
    // (or Unfulfilled if its inputs are not satisfied) and the same
    // information is cascaded to every transitive downstream cell.
    void mark_modified(CellId);

    // The cell ran successfully.  It becomes UpToDate.  Every transitive
    // downstream cell that was UpToDate or Failed now drops to
    // RefreshRequired -- their inputs just changed under their feet.
    // Downstream cells that were already Unfulfilled stay that way until
    // their required inputs are wired.
    void mark_solved(CellId);

    // The cell errored.  It becomes Failed.  Every transitive downstream
    // cell falls to Unfulfilled (their producer is broken).
    void mark_failed(CellId);

private:
    // The state a cell *would* settle into given (a) whether its required
    // inputs are wired and (b) its current upstream cells' states.  Does
    // NOT touch the cell itself; pure function.
    [[nodiscard]] CellState derive_state(CellId) const;

    Schematic* schematic_;
};

}  // namespace simall::workbench
