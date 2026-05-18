// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/CgnsWriter.hpp
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"
#include <string>

namespace simall::io {

struct CgnsWriteResult {
    bool        ok = false;
    std::string error;
    std::string backend;
};

[[nodiscard]] CgnsWriteResult write_cgns_native(const std::string& path,
                                                 const ImportedMesh& mesh);

// If libcgns is linked, this writes a real CGNS file; otherwise it
// forwards to write_cgns_native() and labels the backend "cgns_native".
[[nodiscard]] CgnsWriteResult write_cgns(const std::string& path,
                                          const ImportedMesh& mesh);

}  // namespace simall::io
