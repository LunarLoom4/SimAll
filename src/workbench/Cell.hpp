// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Cell.hpp
// Phase  : 22 Pass 22.1
//
// Cell -- a single node in the project schematic.  Owns its ports, its
// kind, its current state, and a human-readable label.  Cells are minted
// (and ID-stamped) by Schematic::add_cell().
// =============================================================================
#pragma once

#include "workbench/CellPort.hpp"
#include "workbench/Workbench.hpp"

#include <string>
#include <vector>

namespace simall::workbench
{

class Cell
{
public:
    Cell(CellId id, CellKind kind, std::string label);

    // ----- identity ----------------------------------------------------
    [[nodiscard]] CellId id() const noexcept { return id_; }
    [[nodiscard]] CellKind kind() const noexcept { return kind_; }
    [[nodiscard]] const std::string& label() const noexcept { return label_; }
    void set_label(std::string l) { label_ = std::move(l); }

    // ----- state -------------------------------------------------------
    [[nodiscard]] CellState state() const noexcept { return state_; }
    void set_state(CellState s) noexcept { state_ = s; }

    // ----- adapter binding (Pass 22.5) ---------------------------------
    // Each cell may carry an opaque adapter identifier string, e.g.
    // "cad.import.step" or "solver.run.simple".  WorkflowEngine looks
    // this id up in a CellAdapterRegistry to obtain the concrete
    // back-end driver that turns the cell's `RefreshRequired` state
    // into a real CAD/Mesh/Solver/Results invocation.  An empty string
    // means "no adapter" -- the cell is treated as a passive data node
    // and the engine reports a failure when asked to run it.
    [[nodiscard]] const std::string& adapter_id() const noexcept { return adapter_id_; }
    void set_adapter_id(std::string id) { adapter_id_ = std::move(id); }

    // ----- ports -------------------------------------------------------
    // Returns the minted port id.  Names need not be unique across input
    // and output, but must be unique within their direction (debug check).
    PortId add_input(std::string name, std::string data_type, bool required = true);
    PortId add_output(std::string name, std::string data_type);

    [[nodiscard]] const std::vector<CellPort>& ports() const noexcept { return ports_; }
    [[nodiscard]] const CellPort* find_port(PortId) const noexcept;

private:
    CellId id_;
    CellKind kind_;
    std::string label_;
    std::string adapter_id_;
    CellState state_{CellState::Unfulfilled};
    std::vector<CellPort> ports_;
    PortId next_port_id_{0};
};

} // namespace simall::workbench
