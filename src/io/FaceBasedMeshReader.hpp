// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/FaceBasedMeshReader.hpp
// Phase  : 23 Pass 4 (renamed in Pass 5 follow-up)
//
// Importer for the "face-based polyhedral mesh directory" layout: a folder
// holding five plain-ASCII files that together describe an unstructured
// mesh through its faces (rather than directly through its cells).
//
// File set
// --------
//   points     -- N (x y z) coordinates wrapped in a `( ... )` list
//   faces      -- N face vertex lists, each `K(v0 v1 ... vK-1)`
//   owner      -- per-face owner-cell index            (size == faces)
//   neighbour  -- per-internal-face neighbour-cell idx (size <= faces;
//                                                        boundary faces
//                                                        start at index
//                                                        neighbour.size())
//   boundary   -- list of patches: `name { type X; nFaces N; startFace S; }`
//
// Each file optionally begins with a `<HeaderWord> { ... }` dictionary
// header which we tokenize and discard for portability with producers
// that emit such a preamble; headerless variants are also accepted.
// Comments `// ...` and `/* ... */` are supported anywhere.
//
// SimAll's in-memory mesh model is element-based (ImportedMesh ->
// UnstructuredZone -> ElementSection), so a face-based input is
// **converted** at import time:
//
//   1. Per-cell face-lists are reconstructed from (owner, neighbour) by
//      inverting them: each cell `c` owns the faces with owner==c and
//      neighbours the faces with neighbour==c.
//
//   2. Each cell is classified by its face-size signature:
//        4 tri-faces                       -> Tetra4
//        6 quad-faces                      -> Hexa8
//        4 tri-faces + 1 quad-face         -> Pyra5
//        2 tri-faces + 3 quad-faces        -> Penta6  (wedge)
//      Anything else is rejected with an actionable error.  Arbitrary
//      polyhedral cells are out of scope for this pass.
//
//   3. Volume connectivity is emitted with "first-appearance" node
//      ordering across the cell's face vertex lists.  Downstream
//      consumers must rebuild canonical face connectivity from
//      meshing::Connectivity; the section name carries the suffix
//      "_unordered" to signal this.
//
// Boundary patches are read verbatim: one BoundaryPatch and one
// ElementSection per patch, the section holding the patch faces.  Each
// patch must be face-type-homogeneous (all tris or all quads); mixed
// patches are rejected.  BoundaryPatch::bcType receives the patch's
// `type` field (defaults to "patch" if absent).
//
// Compatibility note: this 5-file layout is the same one used by some
// open-source CFD pre-processors that emit a face-based polyhedral
// mesh; nothing in SimAll's architecture is built around it.  The
// reader is purely a file-format adapter.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io {

struct FaceBasedMeshReadResult {
    bool          ok        = false;
    std::string   error;
    ImportedMesh  mesh;
    std::uint8_t  dimension = 0;
};

/// Read a face-based mesh directory from disk.  `dirPath` must contain
/// the five ASCII files `points`, `faces`, `owner`, `neighbour`,
/// `boundary` (the latter two may be empty for a single-cell mesh).
[[nodiscard]] FaceBasedMeshReadResult
read_face_based_mesh_dir(const std::string& dirPath);

/// Test hook: parse from in-memory buffers.
[[nodiscard]] FaceBasedMeshReadResult parse_face_based_mesh_strings(
    const std::string& pointsText,
    const std::string& facesText,
    const std::string& ownerText,
    const std::string& neighbourText,
    const std::string& boundaryText,
    std::string        sourceHint = "<string>");

/// Translate a per-cell face-size signature (count of tri faces and
/// count of quad faces) into an ElementType.  Returns
/// ElementType::Unknown for unsupported signatures.  Exposed for tests.
[[nodiscard]] ElementType
classify_face_based_cell(std::uint32_t numTriFaces,
                         std::uint32_t numQuadFaces) noexcept;

}  // namespace simall::io
