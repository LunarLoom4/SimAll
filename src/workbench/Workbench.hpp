// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/Workbench.hpp
// Phase  : 22 Pass 22.1 (Project Schematic + Workflow State Machine)
//
// Core enums and identifier types shared by every workbench data class.
// Mirrors the ANSYS-Workbench project-schematic concept (Fluent-TG Ch.1,
// CFD-SG "Workflow"): a directed acyclic graph of "cells" representing
// project steps (Geometry -> Mesh -> Setup -> Solution -> Results) where
// each cell exposes typed ports and carries a propagation-aware state.
//
// This header is intentionally pure data + utility -- no UI, no Qt, no I/O.
// =============================================================================
#pragma once

#include <cstdint>
#include <limits>
#include <string_view>

namespace simall::workbench
{

// ---------------------------------------------------------------------------
// Identifier types.  IDs are local to a Schematic instance and are minted by
// the schematic when add_cell() / add_input() / add_output() is called.
// ---------------------------------------------------------------------------
using CellId = std::uint32_t;
using PortId = std::uint16_t;

inline constexpr CellId kInvalidCellId = 0;
inline constexpr PortId kInvalidPortId = std::numeric_limits<PortId>::max();


// ---------------------------------------------------------------------------
// CellKind -- the canonical Workbench-style project steps, plus a Custom
// escape hatch for plugin-defined cells.
// ---------------------------------------------------------------------------
enum class CellKind : std::uint8_t
{
    Geometry,
    Mesh,
    Setup,
    Solution,
    Results,
    Custom,
};


// ---------------------------------------------------------------------------
// CellState -- the four propagation states from the Workbench UI.  Ordering
// is deliberate: numeric value reflects "how bad", so std::max(a, b) gives
// the more conservative of two upstream states when collapsing required
// inputs into a single child state in StateMachine.
//
//   Unfulfilled     : a required input is missing OR an upstream cell is
//                     Unfulfilled or Failed.  Cell cannot run.
//   RefreshRequired : inputs are wired and reachable, but data is stale --
//                     either this cell has not been run since its inputs
//                     last changed, or it was just freshly marked modified.
//   UpToDate        : last run succeeded and inputs have not changed since.
//   Failed          : last run errored.  Treated as a hard-stop for
//                     downstream cells (they fall back to Unfulfilled).
// ---------------------------------------------------------------------------
enum class CellState : std::uint8_t
{
    UpToDate = 0,
    RefreshRequired = 1,
    Unfulfilled = 2,
    Failed = 3,
};


// ---------------------------------------------------------------------------
// PortDirection -- self-explanatory.  An "Input" port is wired by exactly
// one upstream output (no fan-in); an "Output" port may feed many inputs.
// ---------------------------------------------------------------------------
enum class PortDirection : std::uint8_t
{
    Input,
    Output,
};


// ---------------------------------------------------------------------------
// Diagnostic stringification.  Kept noexcept + string_view so it is safe in
// logger hot paths and gtest/Catch2 assertion messages.
// ---------------------------------------------------------------------------
[[nodiscard]] std::string_view to_string(CellKind) noexcept;
[[nodiscard]] std::string_view to_string(CellState) noexcept;
[[nodiscard]] std::string_view to_string(PortDirection) noexcept;

} // namespace simall::workbench
