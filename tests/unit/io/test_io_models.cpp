// =============================================================================
// SimAll Beta - tests/unit/io/test_io_models.cpp
// Week 17 - IO completion + project graph.
//
// Covers FluentMshReader / Plot3dReader / Cgns native round-trip / Tecplot
// ASCII writer / Hdf5ResultStore native SRS / ProjectGraph commit & topo /
// ProjectSerializer round-trip / MovieEncoder APNG header.
// =============================================================================
#include "io/CgnsReader.hpp"
#include "io/CgnsWriter.hpp"
#include "io/FluentMshReader.hpp"
#include "io/Hdf5ResultStore.hpp"
#include "io/MeshFormats.hpp"
#include "io/MovieEncoder.hpp"
#include "io/Plot3dReader.hpp"
#include "io/ProjectGraph.hpp"
#include "io/ProjectSerializer.hpp"
#include "io/TecplotWriter.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;
using namespace simall::io;

namespace
{
std::string tmp_path(const char* tag)
{
    auto p = fs::temp_directory_path() / (std::string("simall_w17_") + tag);
    return p.string();
}
std::string slurp(const std::string& path)
{
    std::ifstream f(path, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}
} // namespace

TEST_CASE("FluentMshReader parses synthetic ASCII msh", "[io][fluent]")
{
    // Minimal synthetic Fluent ASCII file:
    //   (0 "comment")
    //   (2 3)               -- 3D
    //   (10 (0 1 4 1 3))    -- nodes header: zone 0, ids 1..4, type 1, dim 3
    //   (10 (1 1 4 1 3)(
    //     0.0 0.0 0.0
    //     1.0 0.0 0.0
    //     0.0 1.0 0.0
    //     0.0 0.0 1.0))
    //   (12 (0 1 1 1 4))    -- cells header (1 tet)
    //   (12 (2 1 1 1 4))    -- 1 cell of type Tetra4 (cell-type 4)
    const char* text = "(0 \"hello\")\n"
                       "(2 3)\n"
                       "(10 (0 1 4 1 3))\n"
                       "(10 (1 1 4 1 3)(\n"
                       "  0.0 0.0 0.0\n"
                       "  1.0 0.0 0.0\n"
                       "  0.0 1.0 0.0\n"
                       "  0.0 0.0 1.0))\n"
                       "(12 (0 1 1 1 4))\n"
                       "(12 (2 1 1 1 4))\n";
    auto r = parse_fluent_msh_string(text, "<test>");
    REQUIRE(r.ok);
    REQUIRE(r.dimension == 3);
    REQUIRE(!r.mesh.zones.empty());
    REQUIRE(r.mesh.zones[0].x.size() == 4);
    REQUIRE(r.mesh.zones[0].x[1] == 1.0);
}

TEST_CASE("Plot3D formatted multi-block parse", "[io][plot3d]")
{
    const std::string path = tmp_path("grid.p3d");
    {
        std::ofstream f(path);
        f << "2\n";
        f << "2 2 1\n";
        f << "2 2 1\n";
        // block 0 x,y,z
        for (int i = 0; i < 4; ++i)
            f << double(i) << " ";
        f << "\n";
        for (int i = 0; i < 4; ++i)
            f << double(i * 2) << " ";
        f << "\n";
        for (int i = 0; i < 4; ++i)
            f << 0.0 << " ";
        f << "\n";
        // block 1
        for (int i = 0; i < 4; ++i)
            f << double(10 + i) << " ";
        f << "\n";
        for (int i = 0; i < 4; ++i)
            f << double(20 + i) << " ";
        f << "\n";
        for (int i = 0; i < 4; ++i)
            f << 1.0 << " ";
        f << "\n";
    }
    Plot3dOptions opt;
    opt.forceFormatted = true;
    auto r = read_plot3d(path, opt);
    REQUIRE(r.ok);
    REQUIRE(r.grid.blocks.size() == 2);
    REQUIRE(r.grid.blocks[0].ni == 2);
    REQUIRE(r.grid.blocks[0].nj == 2);
    REQUIRE(r.grid.blocks[0].nk == 1);
    REQUIRE(r.grid.blocks[1].z[3] == 1.0);
    std::remove(path.c_str());
}

TEST_CASE("CGNS native writer/reader round-trip", "[io][cgns]")
{
    ImportedMesh m;
    m.sourceFormat = "synthetic";
    UnstructuredZone z;
    z.name = "fluid";
    z.x = {0.0, 1.0, 0.0, 0.0};
    z.y = {0.0, 0.0, 1.0, 0.0};
    z.z = {0.0, 0.0, 0.0, 1.0};
    ElementSection s;
    s.name = "Tetras";
    s.type = ElementType::Tetra4;
    s.nodes = {0, 1, 2, 3};
    z.sections.push_back(s);
    BoundaryPatch bp;
    bp.name = "wall";
    bp.bcType = "wall";
    bp.faceElementIndices = {0};
    z.boundaries.push_back(bp);
    m.zones.push_back(z);

    const std::string path = tmp_path("native.cgns");
    auto w = write_cgns(path, m);
    REQUIRE(w.ok);

    auto r = read_cgns(path);
    REQUIRE(r.ok);
    REQUIRE(r.backend == "cgns_native");
    REQUIRE(r.mesh.zones.size() == 1);
    REQUIRE(r.mesh.zones[0].x.size() == 4);
    REQUIRE(r.mesh.zones[0].sections.size() == 1);
    REQUIRE(r.mesh.zones[0].sections[0].type == ElementType::Tetra4);
    REQUIRE(r.mesh.zones[0].boundaries.size() == 1);
    REQUIRE(r.mesh.zones[0].boundaries[0].bcType == "wall");
    std::remove(path.c_str());
}

TEST_CASE("Tecplot ASCII writer emits expected header tokens", "[io][tecplot]")
{
    ImportedMesh m;
    UnstructuredZone z;
    z.name = "blk";
    z.x = {0, 1, 0, 0};
    z.y = {0, 0, 1, 0};
    z.z = {0, 0, 0, 1};
    ElementSection s;
    s.type = ElementType::Tetra4;
    s.nodes = {0, 1, 2, 3};
    z.sections.push_back(s);
    m.zones.push_back(z);

    const std::string path = tmp_path("out.dat");
    auto r = write_tecplot_ascii(path, m);
    REQUIRE(r.ok);
    auto body = slurp(path);
    REQUIRE(body.find("ZONE T=\"blk\"") != std::string::npos);
    REQUIRE(body.find("FETETRAHEDRON") != std::string::npos);
    std::remove(path.c_str());
}

TEST_CASE("Hdf5ResultStore native SRS round-trip", "[io][srs]")
{
    Hdf5ResultStore s;
    s.write_field("/run/iter_000001",
                  "pressure",
                  std::vector<double>{1.0, 2.0, 3.0, 4.0},
                  FieldLocation::Cell);
    s.root().set_attribute("/run", "solver", "SIMPLE");
    const std::string path = tmp_path("results.srs");
    REQUIRE(s.save_native(path));

    Hdf5ResultStore t;
    REQUIRE(t.load_native(path));
    auto v = t.read_field_f64("/run/iter_000001", "pressure");
    REQUIRE(v.size() == 4);
    REQUIRE(v[2] == 3.0);
    auto* g = t.root().find_group("/run");
    REQUIRE(g != nullptr);
    REQUIRE(g->attributes["solver"] == "SIMPLE");
    std::remove(path.c_str());
}

TEST_CASE("ProjectGraph hash commit invalidates downstream", "[io][graph]")
{
    ProjectGraph g;
    auto* geom = g.add_node(NodeKind::Geometry, "geom");
    auto* mesh = g.add_node(NodeKind::Mesh, "mesh");
    auto* sol = g.add_node(NodeKind::Solution, "sol");
    g.link(geom, mesh);
    g.link(mesh, sol);

    const std::string a = "geomBlobA";
    g.commit(geom, a.data(), a.size());
    g.commit(mesh, "meshBlob1", 9);
    g.commit(sol, "solBlob1", 8);

    REQUIRE(geom->state() == NodeState::Valid);
    REQUIRE(mesh->state() == NodeState::Valid);
    REQUIRE(sol->state() == NodeState::Valid);

    // Mutate geometry → both mesh and sol must be Stale.
    const std::string b = "geomBlobB";
    g.commit(geom, b.data(), b.size());
    REQUIRE(geom->state() == NodeState::Valid);
    REQUIRE(mesh->state() == NodeState::Stale);
    REQUIRE(sol->state() == NodeState::Stale);

    REQUIRE_FALSE(g.uuid_of(geom).empty());
    REQUIRE(g.find_by_uuid(g.uuid_of(geom)) == geom);
    REQUIRE(g.topological_order().size() == 3);
}

TEST_CASE("ProjectSerializer chunked .simall round-trip", "[io][simall]")
{
    ProjectSerializer w;
    ProjectChunk c;
    c.mime = "application/vnd.simall.mesh.cgns";
    c.bytes = {0xDE, 0xAD, 0xBE, 0xEF};
    w.put("mesh.cgns", c);
    ProjectChunk c2;
    c2.mime = "text/plain";
    c2.bytes = {'h', 'i'};
    w.put("readme.txt", c2);

    const std::string path = tmp_path("proj.simall");
    REQUIRE(w.save(path));

    ProjectSerializer r;
    REQUIRE(r.load(path));
    REQUIRE(r.size() == 2);
    auto* g = r.get("mesh.cgns");
    REQUIRE(g != nullptr);
    REQUIRE(g->mime == "application/vnd.simall.mesh.cgns");
    REQUIRE(g->bytes.size() == 4);
    REQUIRE(g->bytes[0] == 0xDE);
    std::remove(path.c_str());
}

TEST_CASE("APNG encoder emits PNG signature + acTL", "[io][movie]")
{
    const std::uint32_t w = 4, h = 4;
    std::vector<std::uint8_t> rgb(w * h * 3, 0x80);
    std::vector<std::vector<std::uint8_t>> frames = {rgb, rgb};
    auto bytes = apng_encode(w, h, 30, true, frames);
    REQUIRE(bytes.size() > 100);
    // PNG signature
    REQUIRE(bytes[0] == 0x89);
    REQUIRE(bytes[1] == 'P');
    REQUIRE(bytes[2] == 'N');
    REQUIRE(bytes[3] == 'G');
    // acTL appears somewhere in the stream
    auto haystack = std::string(bytes.begin(), bytes.end());
    REQUIRE(haystack.find("acTL") != std::string::npos);
    REQUIRE(haystack.find("IHDR") != std::string::npos);
    REQUIRE(haystack.find("IDAT") != std::string::npos);
    REQUIRE(haystack.find("IEND") != std::string::npos);
}

TEST_CASE("MovieEncoder APNG end-to-end", "[io][movie]")
{
    const std::uint32_t w = 8, h = 8;
    const std::string path = tmp_path("movie.apng");
    MovieEncoder enc;
    MovieOptions opt;
    opt.codec = MovieCodec::Apng;
    opt.fps = 24;
    REQUIRE(enc.begin(path, w, h, opt));
    std::vector<std::uint8_t> frame(w * h * 3, 0);
    for (int k = 0; k < 3; ++k) {
        std::fill(frame.begin(), frame.end(), std::uint8_t(k * 64));
        REQUIRE(enc.add_frame_rgb8(frame.data(), frame.size()));
    }
    REQUIRE(enc.end());
    REQUIRE(enc.frame_count() == 3);
    REQUIRE(enc.active_codec() == MovieCodec::Apng);
    REQUIRE(fs::exists(path));
    REQUIRE(fs::file_size(path) > 50);
    std::remove(path.c_str());
}
