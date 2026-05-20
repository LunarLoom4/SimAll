// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/MeshFormats.hpp
// Week   : 17 (IO completion)
//
// Lightweight intermediate types every reader / writer in the IO subsystem
// targets.  Reading a Fluent .msh, a CGNS file, or a Plot3D grid produces
// one of these objects; the caller (`meshing::Mesh::build_from_imported`)
// then synthesises the SoA face connectivity asynchronously.
//
// We deliberately keep these types isolated from `meshing::Mesh` so the IO
// readers compile standalone (zero meshing dep) and can be unit-tested on
// any CI worker.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::io
{

// -- shared scalar types -----------------------------------------------------
using NodeIdx = std::uint32_t;
using CellIdx = std::uint32_t;
using ZoneIdx = std::uint32_t;

// -- element topology -------------------------------------------------------
// Matches the CGNS element-type ordering for the small subset we support;
// CGNS element type numbers in parentheses for cross-reference.
enum class ElementType : std::uint8_t
{
    Unknown = 0,
    Bar2 = 1,   // (3)  edge
    Tri3 = 2,   // (5)  surface
    Quad4 = 3,  // (7)  surface
    Tetra4 = 4, // (10) volume
    Pyra5 = 5,  // (12)
    Penta6 = 6, // (14) wedge / prism
    Hexa8 = 7,  // (17)
    Poly = 8,   // generic n-gon / polyhedron
};

[[nodiscard]] std::uint8_t vertices_per_element(ElementType t) noexcept;
[[nodiscard]] const char* element_type_name(ElementType t) noexcept;

// -- element section --------------------------------------------------------
// One contiguous block of like-typed elements (CGNS / Plot3D convention).
// For Poly we use the CSR `polyOffsets`/`polyNodes` arrays; for everything
// else we use the dense `nodes` array (size = nElements * verts-per-elt).
struct ElementSection
{
    std::string name;
    ElementType type = ElementType::Unknown;
    std::vector<NodeIdx> nodes;             // 0-based dense connectivity
    std::vector<std::uint32_t> polyOffsets; // populated when type==Poly
    std::vector<NodeIdx> polyNodes;         // populated when type==Poly

    [[nodiscard]] std::size_t element_count() const noexcept;
};

// -- boundary patch / named selection ---------------------------------------
struct BoundaryPatch
{
    std::string name;
    std::string bcType;                            // "wall", "inlet", "interface", ...
    std::vector<std::uint32_t> faceElementIndices; // indices into surface ElementSections
};

// -- a single zone (CGNS terminology) of an unstructured mesh ---------------
struct UnstructuredZone
{
    std::string name;
    std::vector<double> x, y, z;          // node coords, 0-based, size = nNodes
    std::vector<ElementSection> sections; // volume + surface mixed
    std::vector<BoundaryPatch> boundaries;
};

// -- a complete imported mesh (potentially multi-zone, e.g. overset) --------
struct ImportedMesh
{
    std::string sourceFormat; // "fluent_msh", "cgns_native", ...
    std::string sourcePath;
    std::vector<UnstructuredZone> zones;

    [[nodiscard]] std::size_t total_nodes() const noexcept;
    [[nodiscard]] std::size_t total_elements() const noexcept;
};

// -- structured (block-structured) grid --------------------------------------
// Plot3D produces these natively; CGNS supports them as well.
struct StructuredBlock
{
    std::string name;
    std::uint32_t ni = 0, nj = 0, nk = 0;
    std::vector<double> x, y, z; // size = ni*nj*nk, k-fastest
    [[nodiscard]] std::size_t point_count() const noexcept;
    [[nodiscard]] std::size_t cell_count() const noexcept;
};

struct StructuredGrid
{
    std::string sourceFormat;
    std::string sourcePath;
    std::vector<StructuredBlock> blocks;
};

// -- field set (solver output) ----------------------------------------------
enum class FieldLocation : std::uint8_t
{
    Node,
    Cell,
    Face
};

struct FieldComponent
{
    std::string name; // "pressure", "velocity_x", ...
    FieldLocation location = FieldLocation::Cell;
    std::vector<double> values; // size = nNodes / nCells / nFaces
};

struct FieldFrame
{
    double time = 0.0;
    std::uint64_t step = 0;
    std::vector<FieldComponent> components;
};

} // namespace simall::io
