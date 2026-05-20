// =============================================================================
// SimAll Beta - IO Unit Tests
// File   : tests/unit/io/test_vtu_reader.cpp
// Phase  : 23 Pass 9
//
// Coverage for the VTK Unstructured Grid (`.vtu`) ASCII reader.
// All tests use parse_vtu_string() to keep the suite disk-free.
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "io/MeshFormats.hpp"
#include "io/VtuReader.hpp"

#include <string>

using simall::io::ElementType;
using simall::io::parse_vtu_string;
using simall::io::VtuReadResult;
using simall::io::vtu_element_type;

namespace {

// Single tetrahedron, 4 vertices.  VTK cell type 10 = VTK_TETRA, 4 nodes.
constexpr const char* kSingleTet = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
          0.0 0.0 0.0
          1.0 0.0 0.0
          0.0 1.0 0.0
          0.0 0.0 1.0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
          0 1 2 3
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
          4
        </DataArray>
        <DataArray type="UInt8" Name="types" format="ascii">
          10
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";

// Mixed mesh: 1 tet + 1 hex + 1 pyramid + 1 wedge.
constexpr const char* kMixed = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1">
  <UnstructuredGrid>
    <Piece NumberOfPoints="13" NumberOfCells="4">
      <Points>
        <DataArray type="Float32" NumberOfComponents="3" format="ascii">
          0 0 0   1 0 0   0 1 0   0 0 1
          2 0 0   3 0 0   3 1 0   2 1 0
          2 0 1   3 0 1   3 1 1   2 1 1
          2.5 0.5 2
        </DataArray>
      </Points>
      <Cells>
        <DataArray Name="connectivity" type="Int64" format="ascii">
          0 1 2 3
          4 5 6 7 8 9 10 11
          4 5 6 7 12
          0 1 2 4 5 6
        </DataArray>
        <DataArray Name="offsets" type="Int64" format="ascii">
          4 12 17 23
        </DataArray>
        <DataArray Name="types" type="UInt8" format="ascii">
          10 12 14 13
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";

}  // namespace

// =============================================================================
// vtu_element_type translator
// =============================================================================
TEST_CASE("vtu_element_type maps standard VTK codes", "[io][vtu]") {
    REQUIRE(vtu_element_type(3)  == ElementType::Bar2);
    REQUIRE(vtu_element_type(5)  == ElementType::Tri3);
    REQUIRE(vtu_element_type(9)  == ElementType::Quad4);
    REQUIRE(vtu_element_type(10) == ElementType::Tetra4);
    REQUIRE(vtu_element_type(12) == ElementType::Hexa8);
    REQUIRE(vtu_element_type(13) == ElementType::Penta6);
    REQUIRE(vtu_element_type(14) == ElementType::Pyra5);

    REQUIRE(vtu_element_type(1)  == ElementType::Unknown);  // VTK_VERTEX
    REQUIRE(vtu_element_type(21) == ElementType::Unknown);  // quadratic edge
    REQUIRE(vtu_element_type(24) == ElementType::Unknown);  // quadratic tet
    REQUIRE(vtu_element_type(42) == ElementType::Unknown);  // polyhedron
    REQUIRE(vtu_element_type(0)  == ElementType::Unknown);
    REQUIRE(vtu_element_type(99) == ElementType::Unknown);
}

// =============================================================================
// Happy-path: single tetrahedron round-trip
// =============================================================================
TEST_CASE("single tetrahedron round-trip", "[io][vtu]") {
    auto r = parse_vtu_string(kSingleTet, "single_tet.vtu");
    INFO(r.error);
    REQUIRE(r.ok);
    REQUIRE(r.error.empty());
    REQUIRE(r.dimension == 3);
    REQUIRE(r.mesh.zones.size() == 1);

    auto const& z = r.mesh.zones.front();
    REQUIRE(z.x.size() == 4);
    REQUIRE(z.y.size() == 4);
    REQUIRE(z.z.size() == 4);
    REQUIRE(z.x[1] == 1.0);
    REQUIRE(z.y[2] == 1.0);
    REQUIRE(z.z[3] == 1.0);

    REQUIRE(z.sections.size() == 1);
    auto const& sec = z.sections.front();
    REQUIRE(sec.type == ElementType::Tetra4);
    REQUIRE(sec.element_count() == 1);
    REQUIRE(sec.nodes.size() == 4);
    REQUIRE(sec.nodes[0] == 0);
    REQUIRE(sec.nodes[3] == 3);
    REQUIRE(sec.name == std::string("vtu_Tetra4"));

    REQUIRE(z.boundaries.empty());
}

// =============================================================================
// Mixed-type cells generate one section per type, in first-appearance order
// =============================================================================
TEST_CASE("mixed cell types group into per-type sections", "[io][vtu]") {
    auto r = parse_vtu_string(kMixed, "mixed.vtu");
    INFO(r.error);
    REQUIRE(r.ok);
    REQUIRE(r.dimension == 3);

    auto const& z = r.mesh.zones.front();
    REQUIRE(z.x.size() == 13);
    REQUIRE(z.sections.size() == 4);

    // Order matches first-appearance order: tet, hex, pyramid, wedge.
    REQUIRE(z.sections[0].type == ElementType::Tetra4);
    REQUIRE(z.sections[0].name == std::string("vtu_Tetra4"));
    REQUIRE(z.sections[0].element_count() == 1);
    REQUIRE(z.sections[0].nodes.size() == 4);

    REQUIRE(z.sections[1].type == ElementType::Hexa8);
    REQUIRE(z.sections[1].element_count() == 1);
    REQUIRE(z.sections[1].nodes.size() == 8);
    REQUIRE(z.sections[1].nodes.front() == 4);
    REQUIRE(z.sections[1].nodes.back()  == 11);

    REQUIRE(z.sections[2].type == ElementType::Pyra5);
    REQUIRE(z.sections[2].element_count() == 1);
    REQUIRE(z.sections[2].nodes.size() == 5);
    REQUIRE(z.sections[2].nodes.back() == 12);

    REQUIRE(z.sections[3].type == ElementType::Penta6);
    REQUIRE(z.sections[3].element_count() == 1);
    REQUIRE(z.sections[3].nodes.size() == 6);
}

// =============================================================================
// 2D mesh (triangles + quads) yields dimension=2
// =============================================================================
TEST_CASE("2D surface mesh has dimension 2", "[io][vtu]") {
    constexpr const char* kSurface = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="5" NumberOfCells="2">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0   1 0 0   0 1 0   1 1 0   2 0 0
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">
          0 1 2 1 4 3
        </DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii"> 3 6 </DataArray>
        <DataArray Name="types"   type="UInt8" format="ascii"> 5 5 </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kSurface, "surface.vtu");
    INFO(r.error);
    REQUIRE(r.ok);
    REQUIRE(r.dimension == 2);
    REQUIRE(r.mesh.zones.front().sections.size() == 1);
    REQUIRE(r.mesh.zones.front().sections[0].type == ElementType::Tri3);
    REQUIRE(r.mesh.zones.front().sections[0].element_count() == 2);
}

// =============================================================================
// VTK_VERTEX cells are silently skipped (carry no topology of interest)
// =============================================================================
TEST_CASE("VTK_VERTEX cells are silently skipped", "[io][vtu]") {
    constexpr const char* kWithVertex = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="2">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0   1 0 0   0 1 0   0 0 1
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">
          0   0 1 2 3
        </DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii"> 1 5 </DataArray>
        <DataArray Name="types"   type="UInt8" format="ascii"> 1 10 </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kWithVertex, "vert.vtu");
    INFO(r.error);
    REQUIRE(r.ok);
    REQUIRE(r.mesh.zones.front().sections.size() == 1);
    REQUIRE(r.mesh.zones.front().sections[0].type == ElementType::Tetra4);
    REQUIRE(r.mesh.zones.front().sections[0].element_count() == 1);
}

// =============================================================================
// XML niceties: comments, CDATA, single-quoted attributes are tolerated
// =============================================================================
TEST_CASE("XML comments, CDATA, and single-quoted attrs are tolerated",
          "[io][vtu]") {
    constexpr const char* kFancy = R"VTU(<?xml version='1.0'?>
<!-- top-level comment -->
<VTKFile type='UnstructuredGrid' version='0.1'>
  <UnstructuredGrid>
    <!-- piece comment -->
    <Piece NumberOfPoints='4' NumberOfCells='1'>
      <Points>
        <DataArray type='Float64' NumberOfComponents='3' format='ascii'>
          <![CDATA[ 0 0 0  1 0 0  0 1 0  0 0 1 ]]>
        </DataArray>
      </Points>
      <Cells>
        <DataArray Name='connectivity' type='Int32' format='ascii'>0 1 2 3</DataArray>
        <DataArray Name='offsets'      type='Int32' format='ascii'>4</DataArray>
        <DataArray Name='types'        type='UInt8' format='ascii'>10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kFancy, "fancy.vtu");
    INFO(r.error);
    REQUIRE(r.ok);
    REQUIRE(r.mesh.zones.front().sections.size() == 1);
    REQUIRE(r.mesh.zones.front().sections[0].type == ElementType::Tetra4);
}

// =============================================================================
// PointData / CellData siblings are ignored (no boundary derivation)
// =============================================================================
TEST_CASE("PointData and CellData siblings are ignored", "[io][vtu]") {
    constexpr const char* kWithFields = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <PointData Scalars="p">
        <DataArray Name="p" type="Float64" format="ascii">
          1.0 2.0 3.0 4.0
        </DataArray>
      </PointData>
      <CellData Scalars="zoneId">
        <DataArray Name="zoneId" type="Int32" format="ascii"> 7 </DataArray>
      </CellData>
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  1 0 0  0 1 0  0 0 1
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">0 1 2 3</DataArray>
        <DataArray Name="offsets"      type="Int32" format="ascii">4</DataArray>
        <DataArray Name="types"        type="UInt8" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kWithFields, "fields.vtu");
    INFO(r.error);
    REQUIRE(r.ok);
    REQUIRE(r.mesh.zones.front().boundaries.empty());
    REQUIRE(r.mesh.zones.front().sections.size() == 1);
}

// =============================================================================
// Errors: malformed / unsupported documents
// =============================================================================
TEST_CASE("error: wrong root element", "[io][vtu]") {
    auto r = parse_vtu_string("<NotVTKFile/>", "bad.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("expected <VTKFile>") != std::string::npos);
}

TEST_CASE("error: wrong dataset type", "[io][vtu]") {
    constexpr const char* kStructured = R"VTU(<?xml version="1.0"?>
<VTKFile type="StructuredGrid"><StructuredGrid/></VTKFile>)VTU";
    auto r = parse_vtu_string(kStructured, "bad.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("StructuredGrid") != std::string::npos);
    REQUIRE(r.error.find("UnstructuredGrid") != std::string::npos);
}

TEST_CASE("error: binary or appended format rejected", "[io][vtu]") {
    constexpr const char* kBinary = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="binary">
        AAAA====
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">0 1 2 3</DataArray>
        <DataArray Name="offsets"      type="Int32" format="ascii">4</DataArray>
        <DataArray Name="types"        type="UInt8" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kBinary, "bin.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("ascii") != std::string::npos);
}

TEST_CASE("error: AppendedData block rejected", "[io][vtu]") {
    constexpr const char* kAppended = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="0" NumberOfCells="0">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="appended" offset="0"/></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="appended" offset="0"/>
        <DataArray Name="offsets"      type="Int32" format="appended" offset="0"/>
        <DataArray Name="types"        type="UInt8" format="appended" offset="0"/>
      </Cells>
    </Piece>
  </UnstructuredGrid>
  <AppendedData encoding="base64">_</AppendedData>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kAppended, "app.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("AppendedData") != std::string::npos);
}

TEST_CASE("error: compressor attribute rejected", "[io][vtu]") {
    constexpr const char* kComp = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" compressor="vtkZLibDataCompressor">
  <UnstructuredGrid/>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kComp, "comp.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("compressor") != std::string::npos);
}

TEST_CASE("error: multi-piece dataset rejected", "[io][vtu]") {
    constexpr const char* kMulti = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="0" NumberOfCells="0">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii"></DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii"></DataArray>
        <DataArray Name="offsets"      type="Int32" format="ascii"></DataArray>
        <DataArray Name="types"        type="UInt8" format="ascii"></DataArray>
      </Cells>
    </Piece>
    <Piece NumberOfPoints="0" NumberOfCells="0">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii"></DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii"></DataArray>
        <DataArray Name="offsets"      type="Int32" format="ascii"></DataArray>
        <DataArray Name="types"        type="UInt8" format="ascii"></DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kMulti, "multi.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("multi-piece") != std::string::npos);
}

TEST_CASE("error: unsupported VTK cell type", "[io][vtu]") {
    constexpr const char* kBad = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="10" NumberOfCells="1">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  1 0 0  0 1 0  0 0 1
        0.5 0 0  0.5 0.5 0  0 0.5 0
        0.5 0 0.5  0 0.5 0.5  0 0 0.5
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">
          0 1 2 3 4 5 6 7 8 9
        </DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii">10</DataArray>
        <DataArray Name="types"   type="UInt8" format="ascii">24</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kBad, "qtet.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("unsupported VTK cell type 24") != std::string::npos);
}

TEST_CASE("error: cell node-count mismatch", "[io][vtu]") {
    constexpr const char* kBad = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  1 0 0  0 1 0  0 0 1
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">0 1 2</DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii">3</DataArray>
        <DataArray Name="types"   type="UInt8" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kBad, "mismatch.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("expects 4 nodes") != std::string::npos);
}

TEST_CASE("error: out-of-range node index in connectivity", "[io][vtu]") {
    constexpr const char* kBad = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  1 0 0  0 1 0  0 0 1
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">0 1 2 99</DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii">4</DataArray>
        <DataArray Name="types"   type="UInt8" format="ascii">10</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kBad, "badidx.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("out-of-range node index 99") != std::string::npos);
}

TEST_CASE("error: missing required DataArray in <Cells>", "[io][vtu]") {
    constexpr const char* kBad = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="1">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  1 0 0  0 1 0  0 0 1
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii">0 1 2 3</DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii">4</DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kBad, "nokind.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("'types'") != std::string::npos);
}

TEST_CASE("error: malformed XML (unterminated tag)", "[io][vtu]") {
    auto r = parse_vtu_string("<VTKFile type=\"UnstructuredGrid\""
                              , "trunc.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("XML parse error") != std::string::npos);
}

TEST_CASE("error: mismatched close tag", "[io][vtu]") {
    auto r = parse_vtu_string(
        "<VTKFile type=\"UnstructuredGrid\"></WrongFile>", "mm.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("mismatched close tag") != std::string::npos);
}

TEST_CASE("error: non-numeric token in points DataArray", "[io][vtu]") {
    constexpr const char* kBad = R"VTU(<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid">
  <UnstructuredGrid>
    <Piece NumberOfPoints="2" NumberOfCells="0">
      <Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">
        0 0 0  oops 0 0
      </DataArray></Points>
      <Cells>
        <DataArray Name="connectivity" type="Int32" format="ascii"></DataArray>
        <DataArray Name="offsets" type="Int32" format="ascii"></DataArray>
        <DataArray Name="types"   type="UInt8" format="ascii"></DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
)VTU";
    auto r = parse_vtu_string(kBad, "nontok.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("non-numeric token") != std::string::npos);
    REQUIRE(r.error.find("oops") != std::string::npos);
}

TEST_CASE("error: missing <UnstructuredGrid>", "[io][vtu]") {
    auto r = parse_vtu_string(
        "<VTKFile type=\"UnstructuredGrid\"></VTKFile>", "empty.vtu");
    REQUIRE_FALSE(r.ok);
    REQUIRE(r.error.find("missing <UnstructuredGrid>") != std::string::npos);
}
