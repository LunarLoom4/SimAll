// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/GambitReader.hpp
// Phase  : 23 Pass 7
//
// Reader for the GAMBIT neutral file format (`.neu`).  GAMBIT was the
// pre-processor that shipped with Fluent through ~v6; the .neu format is
// still the canonical interchange dump for many legacy CFD test cases and
// continues to be emitted by third-party meshers (e.g. snappyToFoam
// converters and certain academic preprocessors).
//
// Section grammar (each section is terminated by the literal `ENDOFSECTION`):
//
//     CONTROL INFO 2.X.Y
//         <title line>
//         PROGRAM: ... VERSION: ...
//         <timestamp line>
//         NUMNP   NELEM   NGRPS   NBSETS   NDFCD   NDFVL
//         <6 integer values>
//     ENDOFSECTION
//
//     NODAL COORDINATES 2.X.Y
//         <id  x  y  [z]>           (NDFCD coords per node; NUMNP rows)
//     ENDOFSECTION
//
//     ELEMENTS/CELLS 2.X.Y
//         <id  type  ndp  n1 ... nN>    (NELEM elements; may wrap lines)
//     ENDOFSECTION
//
//     ELEMENT GROUP 2.X.Y          (repeated NGRPS times)
//         GROUP: <id> ELEMENTS: <m> MATERIAL: <mat> NFLAGS: <nf>
//         <group-name>
//         <nf integer flags>
//         <m element ids>          (10 per line, may wrap)
//     ENDOFSECTION
//
//     BOUNDARY CONDITIONS 2.X.Y    (repeated NBSETS times)
//         <bc-name>  <itype>  <nentry>  <nvalues>  <ibcode1..ibcode5>
//         For itype==0: <node-id> [<nvalues floats>]
//         For itype==1: <cell-id>  <cell-type>  <face-id>
//     ENDOFSECTION
//
// Supported GAMBIT element type codes (mapped via gambit_element_type):
//     1 (Edge, ndp=2)     -> Bar2
//     2 (Quad, ndp=4)     -> Quad4
//     3 (Tri,  ndp=3)     -> Tri3
//     4 (Brick, ndp=8)    -> Hexa8
//     5 (Wedge, ndp=6)    -> Penta6
//     6 (Tet,   ndp=4)    -> Tetra4
//     7 (Pyramid, ndp=5)  -> Pyra5
// High-order (parabolic) ndp values are rejected with a clear error.
//
// Boundary conditions of type 1 (cell+face) are realised by extracting the
// face nodes from the parent volume cell using the GAMBIT face conventions
// and emitting one (or two, for mixed faces) ElementSection per BC patch.
// Boundary conditions of type 0 (node-based) are recorded as patches with
// empty `faceElementIndices` and a synthetic comment in `bcType`.
//
// Element groups are realised as separate volume ElementSections, keyed by
// (groupId, elementType).  Elements not assigned to any group fall into a
// catch-all "ungrouped_<type>" section.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io {

struct GambitReadResult {
    bool          ok        = false;
    std::string   error;
    ImportedMesh  mesh;
    std::uint8_t  dimension = 0;   ///< 2 or 3, copied from CONTROL INFO NDFCD
};

/// Read a GAMBIT neutral file from disk.
[[nodiscard]] GambitReadResult read_gambit_neu(const std::string& path);

/// Test hook: parse from an in-memory buffer.
[[nodiscard]] GambitReadResult parse_gambit_neu_string(const std::string& text,
                                                       std::string        sourceHint
                                                       = "<string>");

/// Translate a GAMBIT element-type code + ndp pair to the SimAll
/// element-type enum.  Returns ElementType::Unknown for unsupported pairs.
[[nodiscard]] ElementType gambit_element_type(int gambitCode, int ndp) noexcept;

}  // namespace simall::io
