// =============================================================================
// SimAll Beta — Solver Unit Tests
// File   : tests/unit/solver/test_wall_distance_service.cpp
// Phase  : 22 Pass 3
//
// Validates the new core::IWallDistanceService interface and its
// solver-side adapter WallDistanceService. Uses the Poisson-Eikonal
// back-end on a 1-D channel with walls at both ends, where the analytic
// distance is  d(x) = min(x, L-x).
//
// Poisson derivation for the 1-D channel of length L with d=0 at x=0,L:
//   -d²φ/dx² = 1, φ(0) = φ(L) = 0   ⇒   φ(x) = x(L-x)/2
//   |dφ/dx| = |L/2 - x|
//   d(x) = -|∇φ| + √(|∇φ|² + 2φ)
//        = -|L/2-x| + √((L/2-x)² + x(L-x))
//        = -|L/2-x| + L/2
//        = min(x, L-x)                        ✓
// =============================================================================
#include "core/IWallDistanceService.hpp"
#include "meshing/MeshStorage.hpp"
#include "solver/LinearSolvers.hpp"
#include "solver/Solver.hpp"
#include "solver/WallDistanceService.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace simall;

namespace
{

// 1-D row of N hex cells (dx × 1 × 1) with WALL zones on BOTH ends.
meshing::Mesh build_channel_walls(int N, double dx)
{
    meshing::Mesh m;
    auto& C = m.cells();
    auto& F = m.faces();

    C.volume.assign(N, dx);
    C.centroidX.resize(N);
    C.centroidY.assign(N, 0.5);
    C.centroidZ.assign(N, 0.5);
    for (int i = 0; i < N; ++i)
        C.centroidX[i] = (i + 0.5) * dx;

    const int Fint = N - 1;
    const int Ftot = Fint + 2;
    F.owner.resize(Ftot);
    F.neighbor.resize(Ftot);
    F.areaX.assign(Ftot, 1.0);
    F.areaY.assign(Ftot, 0.0);
    F.areaZ.assign(Ftot, 0.0);
    F.centroidX.resize(Ftot);
    F.centroidY.assign(Ftot, 0.5);
    F.centroidZ.assign(Ftot, 0.5);
    F.boundaryZone.assign(Ftot, 0);

    for (int i = 0; i < Fint; ++i) {
        F.owner[i] = static_cast<meshing::CellId>(i);
        F.neighbor[i] = static_cast<meshing::CellId>(i + 1);
        F.centroidX[i] = (i + 1) * dx;
    }
    // Left wall (zone 1)
    F.owner[Fint] = 0;
    F.neighbor[Fint] = meshing::kBoundaryCell;
    F.areaX[Fint] = -1.0;
    F.centroidX[Fint] = 0.0;
    F.boundaryZone[Fint] = 1;
    // Right wall (zone 2)
    F.owner[Fint + 1] = static_cast<meshing::CellId>(N - 1);
    F.neighbor[Fint + 1] = meshing::kBoundaryCell;
    F.areaX[Fint + 1] = 1.0;
    F.centroidX[Fint + 1] = N * dx;
    F.boundaryZone[Fint + 1] = 2;

    C.faceOffsets.assign(N + 1, 0);
    for (int i = 0; i < N; ++i) {
        int cnt = 0;
        if (i > 0)
            ++cnt;
        if (i < N - 1)
            ++cnt;
        if (i == 0)
            ++cnt;
        if (i == N - 1)
            ++cnt;
        C.faceOffsets[i + 1] = C.faceOffsets[i] + cnt;
    }
    C.faceIndices.resize(C.faceOffsets[N]);
    std::vector<int> ptr(N, 0);
    for (int i = 0; i < N; ++i) {
        if (i > 0)
            C.faceIndices[C.faceOffsets[i] + ptr[i]++] = i - 1;
        if (i < N - 1)
            C.faceIndices[C.faceOffsets[i] + ptr[i]++] = i;
        if (i == 0)
            C.faceIndices[C.faceOffsets[i] + ptr[i]++] = Fint;
        if (i == N - 1)
            C.faceIndices[C.faceOffsets[i] + ptr[i]++] = Fint + 1;
    }
    return m;
}

std::vector<solver::BoundarySpec> make_two_wall_bcs()
{
    std::vector<solver::BoundarySpec> bcs;
    {
        solver::BoundarySpec b;
        b.zone = 1;
        b.type = solver::BCType::NoSlipWall;
        bcs.push_back(b);
    }
    {
        solver::BoundarySpec b;
        b.zone = 2;
        b.type = solver::BCType::NoSlipWall;
        bcs.push_back(b);
    }
    return bcs;
}

solver::LinearSolverConfig make_lin_cfg()
{
    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::CG;
    c.preconditioner = solver::PreconditionerKind::Jacobi;
    c.tolerance = 1.0e-12;
    c.maxIterations = 1000;
    return c;
}

} // namespace

TEST_CASE("WallDistanceService (Poisson) reproduces 1-D analytic d(x)=min(x,L-x)",
          "[solver][wall-distance][service]")
{
    const int N = 40;
    const double L = 1.0;
    const double dx = L / N;

    auto mesh = build_channel_walls(N, dx);
    auto bcs = make_two_wall_bcs();
    auto lin = solver::make_linear_solver(make_lin_cfg());

    auto svc = solver::make_wall_distance_service_poisson(mesh, bcs, *lin);
    REQUIRE(svc != nullptr);

    // Polymorphic handle — exercise the core interface, not the concrete.
    core::IWallDistanceService& iface = *svc;
    REQUIRE(iface.size() == 0); // before prepare()

    const std::size_t n = iface.prepare();
    REQUIRE(n == static_cast<std::size_t>(N));
    REQUIRE(iface.size() == n);
    REQUIRE(iface.data() != nullptr);
    REQUIRE(std::string(iface.implementation_name()) == "poisson-eikonal");

    // Compare against analytic d(x) = min(x, L-x). Coarse-grid Poisson has
    // 5 % L∞ error against the exact Eikonal solution — tighten if the
    // discretisation improves.
    double linf = 0.0;
    for (int i = 0; i < N; ++i) {
        const double xc = (i + 0.5) * dx;
        const double exact = std::min(xc, L - xc);
        const double got = iface[i];
        const double err = std::fabs(got - exact);
        if (err > linf)
            linf = err;
        CHECK(got >= 0.0);
    }
    CHECK(linf < 0.05 * L);
}

TEST_CASE("WallDistanceService prepare() is idempotent", "[solver][wall-distance][service]")
{
    const int N = 20;
    const double dx = 1.0 / N;

    auto mesh = build_channel_walls(N, dx);
    auto bcs = make_two_wall_bcs();
    auto lin = solver::make_linear_solver(make_lin_cfg());
    auto svc = solver::make_wall_distance_service_poisson(mesh, bcs, *lin);

    const auto n1 = svc->prepare();
    std::vector<double> snap1(svc->data(), svc->data() + n1);

    const auto n2 = svc->prepare();
    std::vector<double> snap2(svc->data(), svc->data() + n2);

    REQUIRE(n1 == n2);
    for (std::size_t i = 0; i < n1; ++i) {
        CHECK(snap1[i] == Catch::Approx(snap2[i]).margin(1.0e-12));
    }
}

TEST_CASE("WallDistanceService exact backend returns zero when no walls present",
          "[solver][wall-distance][service]")
{
    const int N = 10;
    const double dx = 0.1;

    auto mesh = build_channel_walls(N, dx);
    std::vector<solver::BoundarySpec> noWallBcs; // no wall zones at all
    auto svc = solver::make_wall_distance_service_exact(mesh, noWallBcs);

    const auto n = svc->prepare();
    REQUIRE(n == static_cast<std::size_t>(N));
    REQUIRE(std::string(svc->implementation_name()) == "exact-nearest");

    // ExactNearest fallback: with no wall triangles it fills the field with
    // zeros (cannot establish a distance reference, so reports d=0).
    for (std::size_t i = 0; i < n; ++i) {
        CHECK((*svc)[i] == Catch::Approx(0.0).margin(1.0e-12));
    }
}
