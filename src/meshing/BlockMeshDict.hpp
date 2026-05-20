// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/BlockMeshDict.hpp
// Phase  : 23 Pass 14 - OpenFOAM-style blockMeshDict DSL parser.
//
// Parses a textual blockMeshDict subset and produces a populated
// MultiblockHex assembly.  Supported grammar:
//
//   convertToMeters <scalar>;            (optional, default 1.0)
//
//   vertices
//   (
//       (x y z)        // 0
//       (x y z)        // 1
//       ...
//   );
//
//   blocks
//   (
//       hex (v0 v1 v2 v3 v4 v5 v6 v7) <name>? (nx ny nz)
//           simpleGrading (gx gy gz)?   // grading values are PARSED but
//                                       // currently IGNORED (uniform only)
//   );
//
//   boundary
//   (
//       <patchName>
//       {
//           type <patchType>;           // patch/wall/inlet/outlet/empty/...
//           faces
//           (
//               (a b c d)
//               ...
//           );
//       }
//       ...
//   );
//
// Comments:
//   //   line comments
//   /*   block comments  */
//
// Limitations (intentional in this pass):
//   - No edges{} arc/spline support (planar block edges only).
//   - simpleGrading is supported (geometric per-direction expansion);
//     edgeGrading (12 per-edge values) is still rejected.
//   - No mergePatchPairs{}.
//   - No #include / dictionary inheritance.
//
// =============================================================================
#pragma once

#include "meshing/MultiblockHex.hpp"
#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace simall::meshing {

struct BlockMeshBlock {
    std::array<std::uint32_t, 8> vertices{};   // indices into BlockMeshDict::vertices
    std::array<std::uint32_t, 3> divisions{1, 1, 1};
    std::array<double, 3>        grading{1.0, 1.0, 1.0};  // simpleGrading (last/first)
    std::string                  name;                    // optional inline block name
};

struct BlockMeshPatch {
    std::string                                    name;
    std::string                                    type;   // patch/wall/inlet/...
    std::vector<std::array<std::uint32_t, 4>>     faces;
};

struct BlockMeshDict {
    double                       convertToMeters = 1.0;
    std::vector<util::Vec3d>     vertices;
    std::vector<BlockMeshBlock>  blocks;
    std::vector<BlockMeshPatch>  patches;
};

/// Parse a blockMeshDict-format string.  Throws std::runtime_error with a
/// descriptive message (line number, expected token) on syntax error.
BlockMeshDict parse_block_mesh_dict(std::string_view text);

/// Convenience wrapper that reads the file at @p path and parses it.
BlockMeshDict load_block_mesh_dict(const std::string& path);

/// Serialise a BlockMeshDict to a string in the canonical format consumed
/// by parse_block_mesh_dict.  The output is round-trip stable -- feeding
/// it back through the parser yields an equivalent dictionary (same
/// vertices, blocks, patches, grading, divisions, convertToMeters).
std::string write_block_mesh_dict(const BlockMeshDict& dict);

/// Convenience wrapper that writes write_block_mesh_dict() output to @p path.
void save_block_mesh_dict(const BlockMeshDict& dict, const std::string& path);

struct BlockMeshBuildStats {
    std::size_t blocks       = 0;
    std::size_t patches      = 0;
    std::size_t patchFaces   = 0;
    std::size_t orphanFaces  = 0;   // patch faces that did not match any block face
};

/// Populate a MultiblockHex from a parsed dictionary.  Vertices are scaled
/// by convertToMeters, blocks become HexBlock entries (corners taken in the
/// OpenFOAM hex order which matches MultiblockHex's HexBlock convention),
/// and patches are mapped to per-block faceZone tags (1-based, the zone id
/// equals 1 + patchIndex).  Patch faces are matched to block faces by
/// 4-vertex set equality.
BlockMeshBuildStats build_multiblock(const BlockMeshDict& dict,
                                     MultiblockHex&       out);

}  // namespace simall::meshing
