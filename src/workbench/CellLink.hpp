// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/CellLink.hpp
// Phase  : 22 Pass 22.1
//
// CellLink -- directed edge from one cell's output port to another cell's
// input port.  Pure POD; ownership belongs to Schematic.
// =============================================================================
#pragma once

#include "workbench/Workbench.hpp"

namespace simall::workbench
{

struct CellLink
{
    CellId from_cell{kInvalidCellId};
    PortId from_port{kInvalidPortId};
    CellId to_cell{kInvalidCellId};
    PortId to_port{kInvalidPortId};

    [[nodiscard]] friend bool operator==(const CellLink&, const CellLink&) = default;
};

} // namespace simall::workbench
