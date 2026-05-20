// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/VtuReader.hpp
// Phase  : 23 Pass 9
//
// Reader for the VTK Unstructured Grid XML format (`.vtu`, file extension
// for an `UnstructuredGrid` `<VTKFile>` document).  The VTU format is the
// de-facto interchange dump used by ParaView, VisIt, PyVista, meshio and
// many academic CFD codes when round-tripping unstructured meshes.
//
// Grammar (the subset we accept):
//
//     <?xml version="1.0"?>
//     <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
//       <UnstructuredGrid>
//         <Piece NumberOfPoints="N" NumberOfCells="M">
//           <Points>
//             <DataArray type="Float32|Float64" NumberOfComponents="3"
//                        format="ascii"> x1 y1 z1 x2 y2 z2 ... </DataArray>
//           </Points>
//           <Cells>
//             <DataArray Name="connectivity" type="Int*" format="ascii">
//                 n1 n2 ... </DataArray>
//             <DataArray Name="offsets"      type="Int*" format="ascii">
//                 o1 o2 ... </DataArray>
//             <DataArray Name="types"        type="UInt8" format="ascii">
//                 t1 t2 ... </DataArray>
//           </Cells>
//           [<PointData>... ignored ...</PointData>]
//           [<CellData> ... ignored ...</CellData>]
//         </Piece>
//       </UnstructuredGrid>
//     </VTKFile>
//
// Constraints intentionally enforced:
//   * `format="ascii"` only.  `binary` (base64-inlined) and `appended`
//     blob references are rejected with a clear error so the user is
//     directed to re-export ASCII (`meshio convert --ascii`, or
//     ParaView "Save Data" with the `Use ASCII` checkbox).
//   * `compressor=""` on the <VTKFile> root must be absent or empty;
//     a non-empty compressor is rejected (we do not link zlib here).
//   * Only one <Piece> per file is supported (multi-piece datasets are
//     a parallel-output convenience; meshes saved as a single file
//     always have exactly one piece).
//
// VTK cell-type codes mapped to SimAll element types:
//      3  VTK_LINE       -> Bar2
//      5  VTK_TRIANGLE   -> Tri3
//      9  VTK_QUAD       -> Quad4
//     10  VTK_TETRA      -> Tetra4
//     12  VTK_HEXAHEDRON -> Hexa8
//     13  VTK_WEDGE      -> Penta6
//     14  VTK_PYRAMID    -> Pyra5
//      1  VTK_VERTEX     -> silently skipped (carries no topology)
// All other codes (high-order quadratic / polyhedron / polyline / etc.)
// are rejected with a message naming the offending code.
//
// One ElementSection is emitted per unique cell type encountered, named
// `vtu_<TypeName>`.  No boundary patches are derived from a VTU file -
// the VTU format has no first-class boundary concept; downstream tooling
// either reads a sibling `.vtp` surface or interprets a CellData scalar.
// =============================================================================
#pragma once

#include "io/MeshFormats.hpp"

#include <string>

namespace simall::io {

struct VtuReadResult {
    bool          ok        = false;
    std::string   error;
    ImportedMesh  mesh;
    std::uint8_t  dimension = 3;
};

/// Read a `.vtu` file from disk.
[[nodiscard]] VtuReadResult read_vtu(const std::string& path);

/// Test hook: parse from an in-memory buffer.
[[nodiscard]] VtuReadResult parse_vtu_string(const std::string& text,
                                              std::string        sourceHint
                                              = "<string>");

/// Translate a VTK cell-type code to the SimAll element-type enum.
/// Returns ElementType::Unknown for unsupported codes (including
/// VTK_VERTEX=1, which callers should silently skip).
[[nodiscard]] ElementType vtu_element_type(int vtkCellType) noexcept;

}  // namespace simall::io
