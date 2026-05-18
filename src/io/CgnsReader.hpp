// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/CgnsReader.hpp
// Week   : 17
//
// CGNS (CFD General Notation System) reader.  Two backends:
//
//   * SIMALL_HAVE_CGNS=1 — links libcgns (HDF5-backed ADF) and reads any
//     CGNS file produced by Fluent / STAR-CCM+ / Pointwise / etc.
//
//   * Fallback (always present) — reads the self-contained "CGNS-native"
//     mini-format the SimAll writer emits.  Magic `"CGNS-NATIVE\n"`,
//     little-endian length-prefixed records (zone name, coords arrays,
//     element sections, boundary patches).  Used by tests / regression
//     so the round-trip path is testable without libcgns.
//
// The reader sniffs the magic header and dispatches.  Returned `ImportedMesh`
// has the same shape regardless of backend.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"
#include <string>

namespace simall::io {

struct CgnsReadResult {
    bool          ok = false;
    std::string   error;
    std::string   backend;             // "libcgns" or "cgns_native"
    ImportedMesh  mesh;
};

[[nodiscard]] bool           cgns_libcgns_available() noexcept;
[[nodiscard]] CgnsReadResult read_cgns(const std::string& path);

}  // namespace simall::io
