// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/ObjReader.hpp
// Phase  : 4.1 — Wavefront OBJ surface importer.
//
// Wavefront .obj is a polygonal text format widely used as a CAD/mesh
// exchange container.  Like STL, OBJ carries no exact-geometry information,
// so it is decoded into a triangulated TopoDS_Compound (faces are converted
// into triangle fans for n-gons with n > 3) suitable for visualization and
// downstream meshing.  Vertex normals (vn) and texture coordinates (vt)
// are tolerated but not preserved by the CAD graph.
//
// The reader honours the spec's Section 4.2 directive that OBJ be a
// supported import format.  It deliberately keeps the parser self-contained
// (no Assimp dependency) so the CAD module remains buildable without a
// heavy 3rd-party mesh loader.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <string>

namespace simall::cad {

struct ObjReadOptions {
    bool   mergeCoincidentVertices = true;
    double mergeTolerance          = 1.0e-6;
    bool   triangulateNGons        = true;
    bool   recomputeNormals        = true;
};

class ObjReader {
public:
    /// Read an OBJ file and return a ShapeHandle whose TopoDS_Shape is a
    /// compound of triangular faces.  Throws on IO or parse failure.
    ShapeHandle read(const std::string& path, const ObjReadOptions& opts = {});

    /// Fast path that bypasses TopoDS construction and returns the raw
    /// TriangleMesh suitable for direct visualization / regression.
    TriangleMesh readTriangles(const std::string& path,
                               const ObjReadOptions& opts = {});
};

}  // namespace simall::cad
