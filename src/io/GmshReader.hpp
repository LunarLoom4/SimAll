// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/GmshReader.hpp
// Phase  : 23 Pass 1 (v4.1 path); Pass 6 added legacy v2.2 path
//
// Reader for the Gmsh ASCII mesh format.  Supports BOTH the modern
// block-structured v4.x grammar (file-type 0, version 4.0..4.x) and
// the flat-list legacy v2.x grammar (file-type 0, version 2.0..2.x).
// The format version is detected from the `$MeshFormat` header and the
// node/element handlers branch on it.  Optional `$PhysicalNames` and
// (v4 only) `$Entities` sections carry the boundary-condition metadata
// used to derive `BoundaryPatch` records.
//
// Format reference (Gmsh 4.x reference manual, section 9):
//   $MeshFormat
//       version-number  file-type  data-size
//   $EndMeshFormat
//   $PhysicalNames                                (optional)
//       numPhysicalNames
//       dim tag "name"
//       ...
//   $EndPhysicalNames
//   $Entities                                     (optional)
//       numPoints numCurves numSurfaces numVolumes
//       <per-dimension entity records ...>
//   $EndEntities
//   $Nodes
//       numEntityBlocks numNodes minNodeTag maxNodeTag
//       (per block) entityDim entityTag parametric numNodesInBlock
//                   tag\n  ...   (numNodesInBlock lines)
//                   x y z\n ...  (numNodesInBlock lines)
//   $EndNodes
//   $Elements
//       numEntityBlocks numElements minElTag maxElTag
//       (per block) entityDim entityTag elementType numElementsInBlock
//                   tag n1 n2 ... nN\n   (numElementsInBlock lines)
//   $EndElements
//
// Supported Gmsh element type codes (others are reported as an error so the
// caller can refuse to silently drop topology):
//
//     1  -> Bar2   (2-node line)
//     2  -> Tri3   (3-node triangle)
//     3  -> Quad4  (4-node quadrangle)
//     4  -> Tetra4 (4-node tetrahedron)
//     5  -> Hexa8  (8-node hexahedron)
//     6  -> Penta6 (6-node prism / wedge)
//     7  -> Pyra5  (5-node pyramid)
//    15  -> point  (1-node) -- silently skipped (Gmsh emits these for
//                              vertex physical groups; they carry no
//                              topology relevant to a CFD mesh)
//
// Binary `.msh` files (file-type == 1 in the $MeshFormat header) are
// detected and rejected with a clear error message rather than parsed
// incorrectly.
//
// The reader is **single-zone**: every Gmsh file becomes one
// UnstructuredZone whose `sections` vector holds one ElementSection per
// (entityDim, entityTag, elementType) block in v4, or per
// (elementDim, physicalTag, elementType) triple in v2 (where the
// physical tag substitutes for the missing entity hierarchy).  This
// preserves the partition information for downstream zone discovery.
//
// Boundary patches are derived from $PhysicalNames whose dim == dim-1
// where dim is the mesh dimension (so 2 for a 3D mesh, 1 for a 2D mesh).
// Each such physical group becomes one BoundaryPatch whose
// `faceElementIndices` lists the **section indices** (0-based into
// `zone.sections`) of the surface element blocks belonging to entities
// that reference the group. `bcType` defaults to "wall" -- the caller is
// expected to remap on the basis of the patch name.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io
{

struct GmshReadResult
{
    bool ok = false;
    std::string error;
    ImportedMesh mesh;
    std::uint8_t dimension = 0; ///< 1, 2 or 3 (max element dim observed)
};

/// Read a Gmsh ASCII `.msh` file from disk.  Supports v2.x and v4.x;
/// the version is auto-detected from the file's `$MeshFormat` header.
[[nodiscard]] GmshReadResult read_gmsh_msh(const std::string& path);

/// Test hook: parse from an in-memory buffer.
[[nodiscard]] GmshReadResult parse_gmsh_msh_string(const std::string& text,
                                                   std::string sourceHint = "<string>");

/// Translate a Gmsh element-type code (1..7, 15) to the SimAll element
/// type enum. Returns ElementType::Unknown for unsupported codes.
[[nodiscard]] ElementType gmsh_element_type(int gmshCode) noexcept;

} // namespace simall::io
