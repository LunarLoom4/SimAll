// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/TetgenReader.hpp
// Phase  : 23 Pass 8
//
// Reader for Hang Si's TetGen tetrahedral mesh format.  A TetGen "mesh"
// is a set of sibling ASCII files sharing a common path stem:
//
//     <stem>.node  (vertices, mandatory)
//     <stem>.ele   (tetrahedra, mandatory)
//     <stem>.face  (boundary tri faces, optional)
//
// Each file is a token-stream of lines (`#`-prefixed lines are comments):
//
//     .node file:
//         <nNodes> <dim(=3)> <nAttrs> <hasMarker(0|1)>
//         <id> <x> <y> <z> [<attr...>] [<marker>]
//
//     .ele file:
//         <nTets> <nodesPerTet(4|10)> <nAttrs>
//         <id> <n1> <n2> <n3> <n4> [<n5..n10>] [<attr...>]
//
//     .face file:
//         <nFaces> <hasMarker(0|1)>
//         <id> <n1> <n2> <n3> [<marker>]
//
// Both 0-based and 1-based node numbering schemes are produced by TetGen
// in the wild; the reader infers the first-index from whichever number
// appears as the first row's id and uses it to translate every subsequent
// reference.
//
// Only linear tetrahedra are supported (nodesPerTet == 4).  10-node
// parabolic tets (`tetgen -o2`) are rejected with a clear error.
//
// Element attributes are realised as separate ElementSections keyed by
// the integer attribute value (first attribute column only) so that a
// mesh with material regions emerges with one section per region.  Files
// with zero attributes produce a single section named "tetgen_volume".
//
// Boundary face markers become one BoundaryPatch per unique marker; the
// corresponding face nodes land in a Tri3 ElementSection named
// "marker_<id>".  Files without markers produce one anonymous patch
// "all_faces" if any faces are present.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io {

struct TetgenReadResult {
    bool          ok        = false;
    std::string   error;
    ImportedMesh  mesh;
    std::uint8_t  dimension = 3;   ///< TetGen is always 3D
};

/// Read a TetGen mesh by stem path.  Given `stem = "/path/to/mymesh"`,
/// opens `mymesh.node`, `mymesh.ele`, and (if present) `mymesh.face`.
/// Returns ok=false if either of the mandatory files is missing or
/// unreadable; `.face` is optional and its absence is not an error.
[[nodiscard]] TetgenReadResult read_tetgen(const std::string& stemPath);

/// Test hook: parse from in-memory buffers.  `faceText` may be empty
/// (no boundary information).  `sourceHint` is folded into error messages.
[[nodiscard]] TetgenReadResult parse_tetgen_strings(
    const std::string& nodeText,
    const std::string& eleText,
    const std::string& faceText  = {},
    std::string        sourceHint = "<string>");

}  // namespace simall::io
