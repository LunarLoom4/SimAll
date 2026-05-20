// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/StarCdReader.hpp
// Phase  : 23 Pass 3
//
// Reader for the legacy Star-CD / Star-CCM+ ProSTAR ASCII mesh format,
// which is delivered as a triple of paired files written by the
// PROSTAR / pro-am exporter:
//
//     <case>.vrt   -- vertex coordinates
//     <case>.cel   -- cell connectivity (always 8 node slots,
//                                        degenerate hexes encode tet /
//                                        pyramid / wedge via repeated nodes)
//     <case>.bnd   -- boundary face connectivity (4 node slots, triangles
//                                                 duplicate the 4th node)
//
// Grammar of each line in each file (whitespace-separated tokens):
//
//   .vrt : <vert_id> <x> <y> <z>
//   .cel : <cell_id> <n1> ... <n8> [<cell_type> [<region_id>]]
//   .bnd : <bnd_id>  <n1> <n2> <n3> <n4> <region_id>
//                       [<type_name> [<patch_name>...]]
//
// Element-type detection from the 8 .cel slots uses the count of distinct
// node tags in order of first appearance -- the canonical PROSTAR
// collapse pattern guarantees:
//
//     4 unique -> Tetra4   (e.g. 1 2 3 4 4 4 4 4)
//     5 unique -> Pyra5    (e.g. 1 2 3 4 5 5 5 5)
//     6 unique -> Penta6   (e.g. 1 2 3 3 4 5 6 6)
//     8 unique -> Hexa8
//
// Faces in .bnd are similarly classified:
//
//     3 unique -> Tri3
//     4 unique -> Quad4
//
// The reader produces a single UnstructuredZone with at most four
// volume ElementSections (one per element type that actually appears)
// plus one surface ElementSection per unique boundary patch.
// faceElementIndices on each BoundaryPatch lists the *section indices*
// of those surface sections.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io {

struct StarCdReadResult {
    bool          ok        = false;
    std::string   error;
    ImportedMesh  mesh;
    std::uint8_t  dimension = 0;
};

/// Read a Star-CD ASCII mesh from disk. The boundary file is optional --
/// pass an empty string to skip boundary patches.
[[nodiscard]] StarCdReadResult read_starcd(const std::string& vrtPath,
                                            const std::string& celPath,
                                            const std::string& bndPath = "");

/// Test hook: parse from in-memory buffers.
[[nodiscard]] StarCdReadResult parse_starcd_strings(
    const std::string& vrtText,
    const std::string& celText,
    const std::string& bndText   = "",
    std::string        sourceHint = "<string>");

/// Classify a 4 or 8 slot connectivity into an ElementType based on the
/// PROSTAR distinct-node-count convention. Returns ElementType::Unknown
/// for unsupported distinct counts. Out parameter `unique` receives the
/// deduplicated node tags in order of first appearance.
[[nodiscard]] ElementType
starcd_classify_cell(const std::vector<long long>& nodes8,
                     std::vector<long long>&       unique) noexcept;

[[nodiscard]] ElementType
starcd_classify_face(const std::vector<long long>& nodes4,
                     std::vector<long long>&       unique) noexcept;

}  // namespace simall::io
