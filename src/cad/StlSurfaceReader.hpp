// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/StlSurfaceReader.hpp
// Phase  : 4.1 (CAD readers — discrete STL)
//
// STL is a triangulated surface format with no exact geometry. The reader
// loads triangles into a TopoDS_Shape compound of triangular faces, which is
// useful for visualization and meshing but cannot be used by booleans or
// healing operations that demand BRep continuity.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <string>

namespace simall::cad {

struct StlReadOptions {
    bool   mergeCoincidentVertices = true;
    double mergeTolerance          = 1.0e-6;
    bool   buildTopology           = false;  // construct edge/face TopoDS connectivity
};

class StlSurfaceReader {
public:
    /// Read an STL (ASCII or binary). Throws on IO failure.
    ShapeHandle read(const std::string& path, const StlReadOptions& opts = {});

    /// Read STL directly into a TriangleMesh (skips TopoDS construction).
    TriangleMesh readTriangles(const std::string& path,
                               const StlReadOptions& opts = {});
};

}  // namespace simall::cad
