// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Workbench.cpp
// Phase  : 22 Pass 22.1
//
// Diagnostic stringification for the workbench enums.
// =============================================================================
#include "workbench/Workbench.hpp"

namespace simall::workbench {

std::string_view to_string(CellKind k) noexcept {
    switch (k) {
        case CellKind::Geometry: return "Geometry";
        case CellKind::Mesh:     return "Mesh";
        case CellKind::Setup:    return "Setup";
        case CellKind::Solution: return "Solution";
        case CellKind::Results:  return "Results";
        case CellKind::Custom:   return "Custom";
    }
    return "<unknown CellKind>";
}

std::string_view to_string(CellState s) noexcept {
    switch (s) {
        case CellState::UpToDate:        return "UpToDate";
        case CellState::RefreshRequired: return "RefreshRequired";
        case CellState::Unfulfilled:     return "Unfulfilled";
        case CellState::Failed:          return "Failed";
    }
    return "<unknown CellState>";
}

std::string_view to_string(PortDirection d) noexcept {
    switch (d) {
        case PortDirection::Input:  return "Input";
        case PortDirection::Output: return "Output";
    }
    return "<unknown PortDirection>";
}

}  // namespace simall::workbench
