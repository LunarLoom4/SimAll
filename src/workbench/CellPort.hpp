// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/CellPort.hpp
// Phase  : 22 Pass 22.1
//
// CellPort -- a typed input or output socket on a workbench Cell.  Ports
// participate in link validation (output -> input, matching data_type) and
// in the "inputs satisfied" check used by the StateMachine.
// =============================================================================
#pragma once

#include "workbench/Workbench.hpp"

#include <string>

namespace simall::workbench {

struct CellPort {
    PortId        id{kInvalidPortId};
    PortDirection direction{PortDirection::Input};
    std::string   name;
    std::string   data_type;   // free-form tag: "geometry", "mesh", "case", ...
    bool          required{true};   // input-only: blocks state if unwired
};

}  // namespace simall::workbench
