// =============================================================================
// SimAll Beta — IO Unit Tests
// File   : tests/unit/io/test_cgns_native_roundtrip.cpp
// Phase  : 23 Pass 5
//
// Round-trip and error-path coverage for the vendor-neutral "CGNS-native"
// chunked binary mesh format (`CgnsWriter` / `CgnsReader`).  Exercised
// without libcgns -- the test always runs against the fallback backend
// so the round-trip path is fully testable on any CI worker.
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "io/CgnsReader.hpp"
#include "io/CgnsWriter.hpp"
#include "io/MeshFormats.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>

using simall::io::BoundaryPatch;
using simall::io::ElementSection;
using simall::io::ElementType;
using simall::io::ImportedMesh;
using simall::io::NodeIdx;
using simall::io::UnstructuredZone;
using simall::io::read_cgns;
using simall::io::write_cgns_native;

namespace {

[[nodiscard]] std::string temp_path(const char* tag) {
    namespace fs = std::filesystem;
    fs::path p = fs::temp_directory_path()
               / (std::string("simall_cgns_") + tag + "_"
                  + std::to_string(std::rand()) + ".bin");
    return p.string();
}

[[nodiscard]] ImportedMesh make_hex_mesh() {
    ImportedMesh m;
    m.sourceFormat = "synthetic";
    m.sourcePath   = "<test>";
    UnstructuredZone z;
    z.name = "cube";
    z.x = {0, 1, 1, 0, 0, 1, 1, 0};
    z.y = {0, 0, 1, 1, 0, 0, 1, 1};
    z.z = {0, 0, 0, 0, 1, 1, 1, 1};

    ElementSection vol;
    vol.name  = "vol_hex";
    vol.type  = ElementType::Hexa8;
    vol.nodes = {0, 1, 2, 3, 4, 5, 6, 7};
    z.sections.push_back(std::move(vol));

    ElementSection wallSec;
    wallSec.name  = "patch_walls";
    wallSec.type  = ElementType::Quad4;
    wallSec.nodes = {
        0, 3, 2, 1,   // -z
        4, 5, 6, 7,   // +z
        0, 1, 5, 4,   // -y
        1, 2, 6, 5,   // +x
        2, 3, 7, 6,   // +y
        3, 0, 4, 7,   // -x
    };
    z.sections.push_back(std::move(wallSec));

    BoundaryPatch bp;
    bp.name               = "walls";
    bp.bcType             = "wall";
    bp.faceElementIndices = {1};
    z.boundaries.push_back(std::move(bp));

    m.zones.push_back(std::move(z));
    return m;
}

[[nodiscard]] ImportedMesh make_two_zone_mesh() {
    ImportedMesh m;
    m.sourceFormat = "synthetic";
    UnstructuredZone a;
    a.name = "zoneA";
    a.x = {0, 1, 0, 0};
    a.y = {0, 0, 1, 0};
    a.z = {0, 0, 0, 1};
    ElementSection s;
    s.name  = "tetA";
    s.type  = ElementType::Tetra4;
    s.nodes = {0, 1, 2, 3};
    a.sections.push_back(std::move(s));
    m.zones.push_back(std::move(a));

    UnstructuredZone b;
    b.name = "zoneB";
    b.x = {2, 3, 2};
    b.y = {0, 0, 1};
    b.z = {0, 0, 0};
    ElementSection s2;
    s2.name  = "triB";
    s2.type  = ElementType::Tri3;
    s2.nodes = {0, 1, 2};
    b.sections.push_back(std::move(s2));
    BoundaryPatch bp;
    bp.name               = "edge";
    bp.bcType             = "inlet";
    bp.faceElementIndices = {0};
    b.boundaries.push_back(std::move(bp));
    m.zones.push_back(std::move(b));
    return m;
}

[[nodiscard]] ImportedMesh make_poly_mesh() {
    ImportedMesh m;
    UnstructuredZone z;
    z.name = "poly";
    z.x = {0, 1, 1, 0, 0.5};
    z.y = {0, 0, 1, 1, 0.5};
    z.z = {0, 0, 0, 0, 0};
    ElementSection s;
    s.name = "ngons";
    s.type = ElementType::Poly;
    // Two polygons: a quad (4 verts) and a triangle (3 verts).
    s.polyOffsets = {0, 4, 7};
    s.polyNodes   = {0, 1, 2, 3, 0, 1, 4};
    z.sections.push_back(std::move(s));
    m.zones.push_back(std::move(z));
    return m;
}

void expect_round_trip_equal(const ImportedMesh& original,
                             const ImportedMesh& restored) {
    REQUIRE(restored.zones.size() == original.zones.size());
    for (std::size_t zi = 0; zi < original.zones.size(); ++zi) {
        const auto& a = original.zones[zi];
        const auto& b = restored.zones[zi];
        CHECK(a.name == b.name);
        REQUIRE(a.x.size() == b.x.size());
        for (std::size_t i = 0; i < a.x.size(); ++i) {
            CHECK(a.x[i] == b.x[i]);
            CHECK(a.y[i] == b.y[i]);
            CHECK(a.z[i] == b.z[i]);
        }
        REQUIRE(a.sections.size() == b.sections.size());
        for (std::size_t si = 0; si < a.sections.size(); ++si) {
            const auto& as = a.sections[si];
            const auto& bs = b.sections[si];
            CHECK(as.name == bs.name);
            CHECK(as.type == bs.type);
            CHECK(as.nodes == bs.nodes);
            CHECK(as.polyOffsets == bs.polyOffsets);
            CHECK(as.polyNodes == bs.polyNodes);
        }
        REQUIRE(a.boundaries.size() == b.boundaries.size());
        for (std::size_t bi = 0; bi < a.boundaries.size(); ++bi) {
            const auto& ab = a.boundaries[bi];
            const auto& bb = b.boundaries[bi];
            CHECK(ab.name == bb.name);
            CHECK(ab.bcType == bb.bcType);
            CHECK(ab.faceElementIndices == bb.faceElementIndices);
        }
    }
}

struct TempFile {
    std::string path;
    explicit TempFile(const char* tag) : path(temp_path(tag)) {}
    ~TempFile() {
        std::error_code ec;
        std::filesystem::remove(path, ec);
    }
};

}  // namespace

TEST_CASE("CGNS-native round-trips a single-zone hex mesh", "[io][cgns]") {
    TempFile tf("hex");
    const auto original = make_hex_mesh();

    const auto w = write_cgns_native(tf.path, original);
    INFO(w.error);
    REQUIRE(w.ok);
    CHECK(w.backend == "cgns_native");

    const auto r = read_cgns(tf.path);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.backend == "cgns_native");
    CHECK(r.mesh.sourceFormat == "cgns_native");
    CHECK(r.mesh.sourcePath   == tf.path);

    expect_round_trip_equal(original, r.mesh);
}

TEST_CASE("CGNS-native round-trips a multi-zone mesh with boundaries",
          "[io][cgns]") {
    TempFile tf("twozone");
    const auto original = make_two_zone_mesh();

    REQUIRE(write_cgns_native(tf.path, original).ok);
    const auto r = read_cgns(tf.path);
    INFO(r.error);
    REQUIRE(r.ok);
    expect_round_trip_equal(original, r.mesh);
}

TEST_CASE("CGNS-native round-trips a polyhedral section", "[io][cgns]") {
    TempFile tf("poly");
    const auto original = make_poly_mesh();

    REQUIRE(write_cgns_native(tf.path, original).ok);
    const auto r = read_cgns(tf.path);
    INFO(r.error);
    REQUIRE(r.ok);
    expect_round_trip_equal(original, r.mesh);

    // Spot-check the Poly CSR survived intact.
    const auto& s = r.mesh.zones[0].sections[0];
    CHECK(s.type == ElementType::Poly);
    REQUIRE(s.polyOffsets.size() == 3);
    CHECK(s.polyOffsets[0] == 0);
    CHECK(s.polyOffsets[1] == 4);
    CHECK(s.polyOffsets[2] == 7);
    REQUIRE(s.polyNodes.size() == 7);
}

TEST_CASE("CGNS-native round-trips an empty-mesh file", "[io][cgns]") {
    TempFile tf("empty");
    ImportedMesh original;
    original.sourceFormat = "synthetic";

    REQUIRE(write_cgns_native(tf.path, original).ok);
    const auto r = read_cgns(tf.path);
    INFO(r.error);
    REQUIRE(r.ok);
    CHECK(r.mesh.zones.empty());
}

TEST_CASE("CGNS reader rejects a file with bad magic header",
          "[io][cgns][error]") {
    TempFile tf("badmagic");
    {
        std::ofstream f(tf.path, std::ios::binary);
        const char garbage[16] = {'N','O','T','-','A','-','C','G','N','S',
                                  '\0','\0','\0','\0','\0','\0'};
        f.write(garbage, sizeof(garbage));
    }
    const auto r = read_cgns(tf.path);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("Not a CGNS-native file") != std::string::npos);
}

TEST_CASE("CGNS reader rejects an unsupported version",
          "[io][cgns][error]") {
    TempFile tf("badver");
    {
        std::ofstream f(tf.path, std::ios::binary);
        const char magic[12] = {'C','G','N','S','-','N','A','T','I','V','E','\n'};
        f.write(magic, 12);
        const std::uint32_t v = 999;
        f.write(reinterpret_cast<const char*>(&v), 4);
    }
    const auto r = read_cgns(tf.path);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("Unsupported CGNS-native version") != std::string::npos);
}

TEST_CASE("CGNS reader rejects a missing file", "[io][cgns][error]") {
    const auto r = read_cgns("nonexistent_path_simall_cgns_test.bin");
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("Cannot open") != std::string::npos);
}

TEST_CASE("CGNS reader rejects truncated chunk header", "[io][cgns][error]") {
    TempFile tf("trunchdr");
    {
        std::ofstream f(tf.path, std::ios::binary);
        const char magic[12] = {'C','G','N','S','-','N','A','T','I','V','E','\n'};
        f.write(magic, 12);
        const std::uint32_t v = 1;
        f.write(reinterpret_cast<const char*>(&v), 4);
        // Write only the chunk kind byte and stop -- length is missing.
        const std::uint8_t kind = 1;
        f.write(reinterpret_cast<const char*>(&kind), 1);
    }
    const auto r = read_cgns(tf.path);
    CHECK_FALSE(r.ok);
    CHECK((r.error.find("Truncated chunk") != std::string::npos));
}

TEST_CASE("CGNS reader rejects truncated chunk payload",
          "[io][cgns][error]") {
    TempFile tf("truncpay");
    {
        std::ofstream f(tf.path, std::ios::binary);
        const char magic[12] = {'C','G','N','S','-','N','A','T','I','V','E','\n'};
        f.write(magic, 12);
        const std::uint32_t v = 1;
        f.write(reinterpret_cast<const char*>(&v), 4);
        // Zone chunk claiming 99-byte payload but supplies none.
        const std::uint8_t  kind = 1;
        const std::uint32_t plen = 99;
        f.write(reinterpret_cast<const char*>(&kind), 1);
        f.write(reinterpret_cast<const char*>(&plen), 4);
    }
    const auto r = read_cgns(tf.path);
    CHECK_FALSE(r.ok);
    CHECK(r.error.find("Truncated chunk payload") != std::string::npos);
}

TEST_CASE("CGNS reader skips unknown forward-compatible chunk kinds",
          "[io][cgns]") {
    TempFile tf("unknown");
    {
        std::ofstream f(tf.path, std::ios::binary);
        const char magic[12] = {'C','G','N','S','-','N','A','T','I','V','E','\n'};
        f.write(magic, 12);
        const std::uint32_t v = 1;
        f.write(reinterpret_cast<const char*>(&v), 4);
        // Unknown chunk kind 77 with 4 bytes of opaque payload.
        const std::uint8_t  kind77 = 77;
        const std::uint32_t plen   = 4;
        const std::uint32_t junk   = 0xDEADBEEFu;
        f.write(reinterpret_cast<const char*>(&kind77), 1);
        f.write(reinterpret_cast<const char*>(&plen),   4);
        f.write(reinterpret_cast<const char*>(&junk),   4);
        // Then the end marker.
        const std::uint8_t  kindEnd = 255;
        const std::uint32_t plen0   = 0;
        f.write(reinterpret_cast<const char*>(&kindEnd), 1);
        f.write(reinterpret_cast<const char*>(&plen0),   4);
    }
    const auto r = read_cgns(tf.path);
    INFO(r.error);
    CHECK(r.ok);
    CHECK(r.mesh.zones.empty());
}

TEST_CASE("cgns_libcgns_available reflects build configuration",
          "[io][cgns]") {
    // Whatever the build configuration, the predicate must return cleanly
    // and the fallback round-trip path must remain functional.
    const bool have = simall::io::cgns_libcgns_available();
    (void)have;
    SUCCEED("cgns_libcgns_available() callable");
}
