// =============================================================================
// SimAll Beta — IO Subsystem
// File   : src/io/SchematicJson.hpp
// Phase  : 22 Pass 22.3
//
// JSON persistence for a workbench::Schematic (+ optional StateMachine).
// Stable, human-readable, hand-edit-friendly format intended to live inside
// the larger `.simall` project container or stand alone as a `.swb.json`.
//
// Format (schema "simall.workbench.schematic/v1"):
//   { "schema":"simall.workbench.schematic/v1",
//     "cells":[
//       {"id":N,"kind":"Geometry","label":"...","state":"UpToDate",
//        "ports":[{"id":N,"dir":"Input","name":"x","type":"x","required":true},...]}
//     ],
//     "links":[
//       {"from_cell":N,"from_port":N,"to_cell":N,"to_port":N}
//     ] }
//
// Round-trip contract: a Schematic saved with `schematic_to_json` and
// reloaded with `schematic_from_json` yields the same cell labels, kinds,
// ports (by direction + name + type + required), link topology, and
// (when a StateMachine is supplied on load) the same per-cell states.
//
// Note: cell-ids and port-ids are remapped on load — Schematic mints them
// fresh, so the JSON's numeric ids are only used to resolve link
// references during parsing.  The post-load Schematic is structurally
// equivalent, not bit-identical at the id level.
// =============================================================================
#pragma once

#include <string>

namespace simall::workbench { class Schematic; class StateMachine; }

namespace simall::io {

// In-memory string round trip.
[[nodiscard]] std::string schematic_to_json(const workbench::Schematic& s,
                                            const workbench::StateMachine* sm = nullptr);

// Repopulates `out` from JSON.  `out` is cleared / reset first.  When a
// non-null StateMachine pointer is supplied, the parsed per-cell states
// are also written back via Cell::set_state() so the engine sees the
// persisted snapshot.  Returns false on malformed JSON or on any
// validation failure raised by Schematic::add_link.
[[nodiscard]] bool schematic_from_json(const std::string& json,
                                       workbench::Schematic& out,
                                       workbench::StateMachine* sm = nullptr);

// File round trip.
[[nodiscard]] bool save_schematic_json(const std::string& path,
                                       const workbench::Schematic& s,
                                       const workbench::StateMachine* sm = nullptr);
[[nodiscard]] bool load_schematic_json(const std::string& path,
                                       workbench::Schematic& out,
                                       workbench::StateMachine* sm = nullptr);

}  // namespace simall::io
