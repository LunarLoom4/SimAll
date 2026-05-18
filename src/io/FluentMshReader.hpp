// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/FluentMshReader.hpp
// Week   : 17
//
// Reader for the ANSYS Fluent ASCII mesh format (.msh / .cas legacy).
// Implements the subset that covers ~all single-phase CFD models:
//
//   (0  ...comment...)               — discarded
//   (1  ...header...)                — discarded
//   (2  dimension)                   — 2 or 3
//   (10 (zone-id first last type ND) (x y [z] ...))      nodes
//   (12 (zone-id first last type elem-type))             cells
//   (13 (zone-id first last bc-type face-type) (...))   faces + connectivity
//   (39 (zone-id zone-type zone-name) ())               zone names
//
// Binary `.msh` is identified by `(2010 ...)` markers — we report it as an
// error rather than silently emitting wrong topology.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"
#include <string>

namespace simall::io {

struct FluentReadResult {
    bool          ok        = false;
    std::string   error;
    ImportedMesh  mesh;
    std::uint8_t  dimension = 0;
};

[[nodiscard]] FluentReadResult read_fluent_msh(const std::string& path);

// Test hook: parse from an in-memory string.
[[nodiscard]] FluentReadResult parse_fluent_msh_string(const std::string& text,
                                                       std::string sourceHint = "<string>");

}  // namespace simall::io
