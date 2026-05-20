// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Schematic.hpp
// Phase  : 22 Pass 22.1
//
// Schematic -- the directed-acyclic project graph that owns every Cell and
// every CellLink.  Provides:
//
//   * add_cell()           -- mint a Cell with auto-allocated CellId.
//   * add_link()           -- validated link insertion (cycle-free, type-
//                             checked, fan-in-limited on inputs).
//   * topological_order()  -- Kahn's algorithm; deterministic on ties.
//   * upstream / downstream-- direct neighbour queries.
//   * inputs_satisfied()   -- whether every *required* input port on a
//                             cell is wired (used by the StateMachine).
//
// No I/O, no UI, no threading.  StateMachine sits on top of this class.
// =============================================================================
#pragma once

#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/Workbench.hpp"

#include <cstddef>
#include <vector>

namespace simall::workbench {

class Schematic {
public:
    Schematic() = default;

    // ----- cells -------------------------------------------------------
    CellId add_cell(CellKind kind, std::string label);

    // Removes a cell *and every incident link*.  Returns false if the id
    // is unknown.  When `removed_links_out` is non-null it is populated
    // (in registration order) with every link that was implicitly
    // dropped so an undo command can restore them verbatim via
    // `restore_link`.  Note: state propagation (downstream invalidation
    // of formerly-dependent cells) is left to the caller / WorkflowEngine
    // -- Schematic itself stays pure-graph.
    bool remove_cell(CellId id,
                     std::vector<CellLink>* removed_links_out = nullptr);

    // Re-inserts a cell that was previously held by `remove_cell` (or any
    // externally-constructed Cell carrying a stable id).  Fails when the
    // id is already in use.  Bumps next_cell_id_ above the restored id
    // so subsequent `add_cell` calls cannot collide.  This is the only
    // path that injects a non-minted id; commands/undo use it exclusively.
    bool restore_cell(Cell cell);

    // Convenience wrapper for symmetry with restore_cell.  Returns the
    // same bool as add_link.
    bool restore_link(const CellLink& link) { return add_link(link); }

    [[nodiscard]] Cell*                     cell(CellId)       noexcept;
    [[nodiscard]] const Cell*               cell(CellId) const noexcept;
    [[nodiscard]] const std::vector<Cell>&  cells()      const noexcept { return cells_; }
    [[nodiscard]] std::size_t               size()       const noexcept { return cells_.size(); }

    // ----- links -------------------------------------------------------
    // Rejects (returns false) when any of the following holds:
    //   * either cell or port id is unknown
    //   * from_port is not an Output / to_port is not an Input
    //   * data_type strings disagree (both non-empty, mismatched)
    //   * the same link already exists
    //   * to_port already has an upstream link (single-fan-in inputs)
    //   * insertion would create a cycle
    bool add_link(const CellLink& link);
    bool remove_link(const CellLink& link);

    [[nodiscard]] const std::vector<CellLink>& links() const noexcept { return links_; }

    // ----- graph queries ----------------------------------------------
    // Kahn's algorithm.  Returns an order in which every dependency
    // appears before its dependants.  Empty result iff the graph is
    // somehow cyclic (should be unreachable since add_link guards it).
    [[nodiscard]] std::vector<CellId> topological_order() const;

    [[nodiscard]] std::vector<CellId> upstream  (CellId) const;
    [[nodiscard]] std::vector<CellId> downstream(CellId) const;

    [[nodiscard]] bool inputs_satisfied(CellId) const;

private:
    [[nodiscard]] bool would_create_cycle(CellId from, CellId to) const;

    std::vector<Cell>     cells_;
    std::vector<CellLink> links_;
    CellId                next_cell_id_{1};  // 0 reserved as kInvalidCellId
};

}  // namespace simall::workbench
