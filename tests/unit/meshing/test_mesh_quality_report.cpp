// =============================================================================
// SimAll Beta - Meshing Unit Tests
// File   : tests/unit/meshing/test_mesh_quality_report.cpp
// Phase  : 23 Pass 17
//
// Coverage for the extended MeshQuality reporting surface:
//   - format_histogram() emits one labelled line per bin with bar chars
//   - extended format() now contains all three histogram sections
//   - to_csv_faces() row count == nFaces + 1 (header)
//   - to_csv_cells() row count == nCells + 1 (header)
//   - empty-mesh histogram path emits "(no samples)"
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "meshing/MeshQuality.hpp"
#include "meshing/MeshStorage.hpp"
#include "meshing/CartesianMesher.hpp"

#include <algorithm>
#include <numeric>
#include <string>

using namespace simall;

namespace {

std::size_t count_lines(const std::string& s) {
    return static_cast<std::size_t>(std::count(s.begin(), s.end(), '\n'));
}

meshing::Mesh build_brick(int Nx, int Ny, int Nz) {
    meshing::CartesianGridSpec spec;
    spec.origin = {0, 0, 0};
    spec.extent = {1, 1, 1};
    spec.Nx = Nx; spec.Ny = Ny; spec.Nz = Nz;
    meshing::Mesh m;
    meshing::CartesianMesher{spec}.generate(m);
    return m;
}

}  // namespace

// =============================================================================
TEST_CASE("MeshQuality::format includes the three histogram sections",
          "[meshing][quality][histogram]") {
    auto m = build_brick(3, 3, 3);
    const auto r = meshing::MeshQuality::evaluate(m);
    const auto text = meshing::MeshQuality::format(r);

    REQUIRE(text.find("Mesh quality report")        != std::string::npos);
    REQUIRE(text.find("histogram: skewness")        != std::string::npos);
    REQUIRE(text.find("histogram: non-orthogonality (deg)") != std::string::npos);
    REQUIRE(text.find("histogram: aspect ratio")    != std::string::npos);
}

// =============================================================================
TEST_CASE("MeshQuality::format_histogram emits one line per bin plus header",
          "[meshing][quality][histogram]") {
    auto m = build_brick(2, 2, 2);
    const auto r = meshing::MeshQuality::evaluate(m);
    const auto h = meshing::MeshQuality::format_histogram(r.histSkewness, "skew");
    // Expect: 1 "histogram: skew" header line + 10 bin rows  OR
    //         1 header + 1 "(no samples)" line.
    REQUIRE((count_lines(h) == 11 || count_lines(h) == 2));
    REQUIRE(h.find("histogram: skew") != std::string::npos);
}

// =============================================================================
TEST_CASE("MeshQuality::format_histogram on empty bins reports (no samples)",
          "[meshing][quality][histogram]") {
    meshing::QualityHistogram empty;        // all bins zero
    const auto h = meshing::MeshQuality::format_histogram(empty, "x");
    REQUIRE(h.find("(no samples)") != std::string::npos);
}

// =============================================================================
TEST_CASE("MeshQuality::to_csv_faces yields nFaces+1 lines with correct header",
          "[meshing][quality][csv]") {
    auto m = build_brick(2, 2, 2);
    const auto r   = meshing::MeshQuality::evaluate(m);
    const auto csv = meshing::MeshQuality::to_csv_faces(r);

    REQUIRE(csv.rfind("faceId,skewness,nonOrthoDeg", 0) == 0);
    REQUIRE(count_lines(csv) == r.nFaces + 1);
}

// =============================================================================
TEST_CASE("MeshQuality::to_csv_cells yields nCells+1 lines with correct header",
          "[meshing][quality][csv]") {
    auto m = build_brick(2, 2, 2);
    const auto r   = meshing::MeshQuality::evaluate(m);
    const auto csv = meshing::MeshQuality::to_csv_cells(r);

    REQUIRE(csv.rfind("cellId,aspectRatio,negativeVolume", 0) == 0);
    REQUIRE(count_lines(csv) == r.nCells + 1);
}

// =============================================================================
TEST_CASE("MeshQuality::per_zone_stats groups boundary faces by zone",
          "[meshing][quality][zones]") {
    auto m       = build_brick(2, 2, 2);
    const auto r = meshing::MeshQuality::evaluate(m);
    const auto z = meshing::MeshQuality::per_zone_stats(m, r);

    // CartesianMesher creates xMin/xMax/yMin/yMax/zMin/zMax — 6 zones,
    // 4 boundary faces each on a 2x2x2 brick.
    REQUIRE(z.size() == 6);
    std::size_t totalFaces = 0;
    for (const auto& s : z) {
        REQUIRE(s.nFaces == 4);
        REQUIRE(s.meanSkewness <= s.maxSkewness);
        REQUIRE(s.meanNonOrtho <= s.maxNonOrtho);
        REQUIRE_FALSE(s.name.empty());
        totalFaces += s.nFaces;
    }
    REQUIRE(totalFaces == 24);
}

// =============================================================================
TEST_CASE("MeshQuality::format_per_zone emits header row + one row per zone",
          "[meshing][quality][zones]") {
    auto m       = build_brick(2, 2, 2);
    const auto r = meshing::MeshQuality::evaluate(m);
    const auto z = meshing::MeshQuality::per_zone_stats(m, r);
    const auto t = meshing::MeshQuality::format_per_zone(z);

    REQUIRE(t.find("Per-zone quality") != std::string::npos);
    REQUIRE(t.find("xMin") != std::string::npos);
    REQUIRE(t.find("zMax") != std::string::npos);
    // 1 banner + 1 header + 6 data rows = 8 newlines.
    REQUIRE(count_lines(t) == 8);
}

// =============================================================================
TEST_CASE("MeshQuality::format_per_zone empty list reports no zones",
          "[meshing][quality][zones]") {
    std::vector<meshing::ZoneQualityStats> empty;
    const auto t = meshing::MeshQuality::format_per_zone(empty);
    REQUIRE(t.find("(no boundary zones)") != std::string::npos);
}

// =============================================================================
TEST_CASE("MeshQuality::evaluate populates volume stats on a unit brick",
          "[meshing][quality][volume]") {
    auto m       = build_brick(2, 2, 2);   // 8 cells, each volume = 0.125
    const auto r = meshing::MeshQuality::evaluate(m);

    REQUIRE(r.nCells == 8);
    REQUIRE(r.totalVolume == Catch::Approx(1.0).margin(1e-12));
    REQUIRE(r.minVolume   == Catch::Approx(0.125).margin(1e-12));
    REQUIRE(r.maxVolume   == Catch::Approx(0.125).margin(1e-12));
    REQUIRE(r.meanVolume  == Catch::Approx(0.125).margin(1e-12));
    // Uniform cells → all samples collapse into bin 0 (vmin == vmax).
    const std::size_t total =
        std::accumulate(r.histVolume.bins.begin(),
                        r.histVolume.bins.end(), std::size_t{0});
    REQUIRE(total == r.nCells);
}

// =============================================================================
TEST_CASE("MeshQuality::format includes cell volume summary and histogram",
          "[meshing][quality][volume]") {
    auto m       = build_brick(2, 2, 2);
    const auto r = meshing::MeshQuality::evaluate(m);
    const auto t = meshing::MeshQuality::format(r);

    REQUIRE(t.find("cell volume: min=")     != std::string::npos);
    REQUIRE(t.find("histogram: cell volume") != std::string::npos);
}
