// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Commands.hpp
// Phase  : 22 Pass 22.4
//
// Reversible (core::ICommand) wrappers around every Schematic / Workflow
// mutation a user can perform from the schematic panel.  Each command:
//
//   * Snapshots the pre-change state needed to restore the schematic.
//   * Drives the change through Schematic / WorkflowEngine (NOT through
//     direct member writes), so state propagation rules stay consistent
//     with the rest of the workbench.
//   * Optionally records a `Do` / `Undo` line into a ChangeJournal so a
//     UI panel can render session history independent of the redo bucket.
//   * Implements core::ICommand::description() with a stable, human-
//     readable string ("Add Mesh cell 'M3'", "Connect G.out -> M.in").
//
// Public surface is intentionally limited to free factory functions
// returning `std::unique_ptr<core::ICommand>`; the concrete classes live
// in the .cpp so call-site code never depends on their layouts.
// =============================================================================
#pragma once

#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/Workbench.hpp"

#include <memory>
#include <string>

namespace simall::core      { class ICommand; }
namespace simall::workbench { class Schematic; class StateMachine;
                              class WorkflowEngine; class ChangeJournal; }

namespace simall::workbench {

// Add a new (empty-port) cell.  After execute(), the minted CellId is
// retrievable via the post-condition pointer if the caller passed one.
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_add_cell_command(Schematic&    s,
                      CellKind      kind,
                      std::string   label,
                      CellId*       out_minted_id = nullptr,
                      ChangeJournal* journal      = nullptr);

// Remove a cell *and every incident link*.  Undo restores the cell, its
// ports, and every previously-incident link verbatim.  Returns false at
// execute() time if the id is unknown.
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_remove_cell_command(Schematic&     s,
                         CellId         id,
                         ChangeJournal* journal = nullptr);

// Connect from_cell.from_port -> to_cell.to_port through the
// WorkflowEngine so propagation rules fire on both execute / undo.
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_connect_command(WorkflowEngine& eng,
                     CellLink        link,
                     ChangeJournal*  journal = nullptr);

// Inverse of make_connect_command.
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_disconnect_command(WorkflowEngine& eng,
                        CellLink        link,
                        ChangeJournal*  journal = nullptr);

// Rename a cell.
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_set_cell_label_command(Schematic&     s,
                            CellId         id,
                            std::string    new_label,
                            ChangeJournal* journal = nullptr);

// Explicit user override of a cell's CellState (e.g. mark UpToDate after
// importing pre-computed results).  Does NOT re-propagate downstream --
// that is left to the WorkflowEngine on the next mutation.
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_set_cell_state_command(Schematic&     s,
                            CellId         id,
                            CellState      new_state,
                            ChangeJournal* journal = nullptr);

// Bind (or rebind) a cell to a back-end adapter via its dotted id, e.g.
// "cad.import.step".  An empty string unbinds.  Reversible.  Captured
// as a value snapshot in the command so undo/redo are pure-data and
// safe to replay against a registry that mutates between sessions.
// (Pass 22.5.)
[[nodiscard]] std::unique_ptr<simall::core::ICommand>
make_set_cell_adapter_command(Schematic&     s,
                              CellId         id,
                              std::string    new_adapter_id,
                              ChangeJournal* journal = nullptr);

}  // namespace simall::workbench
