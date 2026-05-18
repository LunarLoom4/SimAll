// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/TecplotWriter.hpp
// Week   : 17
//
// Tecplot ASCII (.dat) writer for surface / volume meshes plus an
// arbitrary number of cell- or node-centred fields.  Emits Tecplot 360
// `FEPOLYHEDRON` / `FETETRAHEDRON` / `FEBRICK` zones in BLOCK packing.
//
// Cell-centred fields are declared via `VARLOCATION=([varIndex]=CELLCENTERED)`.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"
#include <string>

namespace simall::io {

struct TecplotWriteResult {
    bool        ok = false;
    std::string error;
};

[[nodiscard]] TecplotWriteResult write_tecplot_ascii(const std::string& path,
                                                      const ImportedMesh& mesh,
                                                      const FieldFrame*   fields = nullptr,
                                                      const std::string&  title  = "SimAll Beta");

}  // namespace simall::io
