// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/UnvReader.hpp
// Phase  : 23 Pass 2
//
// Reader for the SDRC / I-DEAS Universal (".unv") ASCII mesh format. The
// format is the de-facto interchange used by HyperMesh, ANSA, IcemCFD and
// the open-source meshing tools (gmsh exports unv as well), so this reader
// unlocks a large slice of legacy and pre-processor workflows for SimAll.
//
// Grammar (textbook UNV reference, "Universal File Datasets"):
//
//   Every dataset is delimited by a line containing only the sentinel
//   "-1" (right-justified in column 6). Inside the delimiters the first
//   line carries the dataset number, the remaining lines hold the body.
//
//       <whitespace> -1
//       <whitespace> <dataset_number>
//       ... body ...
//       <whitespace> -1
//
// Supported dataset numbers (others are silently skipped, which preserves
// forward-compatibility the same way the Gmsh reader handles unknown
// $Foo sections):
//
//     2411 -- Nodes.   Record-1: label, expCoordSys, dispCoordSys, color.
//                      Record-2: x y z   (Fortran D-format floats).
//     2412 -- Elements. Record-1: label, fe_descriptor, phys_prop_tbl,
//                      mat_prop_tbl, color, num_nodes.
//                      Record-2: beam orientation triple, ONLY for the
//                      beam fe_descriptors 11 / 21 / 22.
//                      Record-3+: connectivity (one or more lines until
//                      num_nodes node labels have been consumed).
//     2467 -- Permanent Groups (legacy).  Header line of 8 ints, then a
//                      name line (40-char field), then ceil(N/2) lines
//                      of 8 ints, where each pair (type tag leafId compId)
//                      describes one group entity.  type==8 means "FE
//                      element"; type==7 means "node".
//     2477 -- Permanent Groups (new, equivalent layout to 2467).
//
// Supported FE descriptor codes (the linear single-element subset which is
// what 99 % of CFD .unv files contain):
//
//     11, 21, 22 -> Bar2    (line / beam, skips the 3-int beam metadata)
//     41, 91     -> Tri3    (plane triangle, thin-shell triangle)
//     44, 94     -> Quad4   (plane quad,      thin-shell quad)
//     111        -> Tetra4
//     112        -> Penta6
//     115        -> Hexa8
//
// Higher-order (parabolic) variants 24, 42, 45, 92, 95, 116, 118 are
// reported as an error rather than silently dropped, matching the
// Phase 23 Pass 1 Gmsh-reader policy.
//
// Boundary patches: any 2467/2477 group composed exclusively of element
// references whose elements all live in surface-dim sections (where
// surfDim == maxElementDim - 1) is promoted to a BoundaryPatch whose
// faceElementIndices lists the **section indices** (0-based into
// `zone.sections`) containing those elements.  Node-only groups are
// stored but not exposed as patches in this pass.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io
{

struct UnvReadResult
{
    bool ok = false;
    std::string error;
    ImportedMesh mesh;
    std::uint8_t dimension = 0;
};

/// Read an Ideas Universal `.unv` ASCII file from disk.
[[nodiscard]] UnvReadResult read_unv(const std::string& path);

/// Test hook: parse from an in-memory buffer.
[[nodiscard]] UnvReadResult parse_unv_string(const std::string& text,
                                             std::string sourceHint = "<string>");

/// Translate a UNV finite-element descriptor code to the SimAll element
/// type enum. Returns ElementType::Unknown for unsupported codes.
[[nodiscard]] ElementType unv_element_type(int feDescriptor) noexcept;

/// Convert a Fortran D-format float token (e.g. "1.500D-02") to an ASCII
/// E-format token the C++ float parser accepts. Exposed for tests.
[[nodiscard]] std::string unv_normalize_float_token(std::string s);

} // namespace simall::io
