// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/BrepIo.hpp
// Phase  : 4.1 (CAD I/O — native OCC BREP)
//
// OpenCASCADE's native lossless format. Round-trips exact BRep + topology
// without translation loss, so SimAll uses it for cache snapshots between
// pipeline stages (e.g. mesh-time snapshot of healed geometry).
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <string>

namespace simall::cad
{

class BrepIo
{
public:
    /// Read a .brep file. Throws on IO error.
    ShapeHandle read(const std::string& path);

    /// Write a ShapeHandle to disk in OCC BREP format. Returns false on IO
    /// error. Compatible with OCC 7.x.
    bool write(const ShapeHandle& shape, const std::string& path) const;
};

} // namespace simall::cad
