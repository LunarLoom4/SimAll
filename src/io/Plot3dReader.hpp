// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/Plot3dReader.hpp
// Week   : 17
//
// Reader for NASA's Plot3D structured-grid format.  Two flavours:
//
//   * ASCII (multi-block, formatted):
//        nblocks
//        ni1 nj1 nk1   ni2 nj2 nk2 ...
//        <coords block 1: x[1..N] y[1..N] z[1..N]>
//        <coords block 2: ...>
//
//   * Unformatted Fortran binary: identical layout but with the standard
//     4-byte record-length prefix/suffix around every record.
//
// We sniff the format by inspecting the first 8 bytes (printable ASCII →
// formatted, otherwise binary).  Both real and double precision are
// supported via an optional flag on the reader options.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"
#include <string>

namespace simall::io {

struct Plot3dOptions {
    bool   doublePrecision = false;
    bool   forceBinary     = false;   // skip auto-detection
    bool   forceFormatted  = false;
};

struct Plot3dReadResult {
    bool            ok = false;
    std::string     error;
    StructuredGrid  grid;
};

[[nodiscard]] Plot3dReadResult read_plot3d(const std::string& path,
                                           Plot3dOptions opts = {});
[[nodiscard]] Plot3dReadResult parse_plot3d_formatted(const std::string& text,
                                                      std::string sourceHint = "<string>");

}  // namespace simall::io
