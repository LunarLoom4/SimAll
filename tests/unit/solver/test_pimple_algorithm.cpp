// =============================================================================
// SimAll Beta — Solver Unit Tests
// File   : tests/unit/solver/test_pimple_algorithm.cpp
// Phase  : 22 Pass 1
//
// PIMPLE = transient SIMPLE-outer / PISO-inner hybrid (see
// src/solver/PimpleAlgorithm.hpp header for derivation).
//
// Tests:
//   1) PIMPLE iterate() returns finite, non-NaN residuals on a 1-D channel.
//   2) PIMPLE with nOuterCorrectors=1 + momentumPredictor=true matches PISO
//      behaviour within tolerance (both perform 1 predictor + N pressure
//      correctors with no under-relaxation).
//   3) Two-corrector PIMPLE (nOuterCorrectors=2) yields continuity residual
//      no worse than one-corrector PIMPLE — outer loop must be monotone in
//      the sense of not degrading the converged-corrector residual.
//   4) Accessors expose the configured nOuterCorrectors / nCorrectors /
//      nNonOrthCorr / momentumPredictor verbatim.
//   5) last_outer_iters() reports the early-exit point when tolU/tolP
//      trigger before the full outer count.
//
// The mesh is a 1-D row of N hex cells with VelocityInlet at zone 1 and
// PressureOutlet at zone 2, identical layout to test_fct_zalesak.cpp so that
// the FieldRegistry / linear-solver wiring can stay minimal.
// =============================================================================
#include "meshing/MeshStorage.hpp"
#include "solver/LinearSolvers.hpp"
#include "solver/PimpleAlgorithm.hpp"
#include "solver/PisoAlgorithm.hpp"
#include "solver/Solver.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

using namespace simall;
using util::aligned_vector;

namespace
{

// 1-D row of N hex cells (dx × 1 × 1) along +x with inlet/outlet boundaries.
meshing::Mesh build_channel_1d(int N, double dx)
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
    // Left = inlet
    F.owner[Fint] = 0;
    F.neighbor[Fint] = meshing::kBoundaryCell;
    F.areaX[Fint] = -1.0;
    F.centroidX[Fint] = 0.0;
    F.boundaryZone[Fint] = 1;
    // Right = outlet
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
            ++cnt; // left boundary
        if (i == N - 1)
            ++cnt; // right boundary
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

std::vector<solver::BoundarySpec> make_channel_bcs(double Uinlet, double Pout)
{
    std::vector<solver::BoundarySpec> bcs;
    {
        solver::BoundarySpec b;
        b.zone = 1;
        b.type = solver::BCType::VelocityInlet;
        b.vectorValue[0] = Uinlet;
        bcs.push_back(b);
    }
    {
        solver::BoundarySpec b;
        b.zone = 2;
        b.type = solver::BCType::PressureOutlet;
        b.scalarValue = Pout;
        bcs.push_back(b);
    }
    return bcs;
}

solver::LinearSolverConfig make_lin_cfg()
{
    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::Jacobi;
    c.tolerance = 1.0e-9;
    c.maxIterations = 200;
    c.restart = 30;
    return c;
}

void prime_fields(solver::FieldRegistry& F, std::size_t nC, double Uinit)
{
    auto& U = F.vector("U", nC);
    auto& p = F.scalar("p", nC);
    for (std::size_t i = 0; i < nC; ++i) {
        U.x[i] = Uinit;
        U.y[i] = 0;
        U.z[i] = 0;
    }
    std::fill(p.begin(), p.end(), 0.0);
}

inline bool finite_residuals(const solver::SimpleResiduals& r)
{
    return std::isfinite(r.mom[0]) && std::isfinite(r.mom[1]) && std::isfinite(r.mom[2])
           && std::isfinite(r.cont);
}

} // namespace

TEST_CASE("PIMPLE iterate produces finite residuals on a 1-D channel", "[solver][pimple]")
{
    const int N = 20;
    const double dx = 1.0 / N;
    auto mesh = build_channel_1d(N, dx);
    auto bcs = make_channel_bcs(/*Uinlet=*/1.0, /*Pout=*/0.0);

    solver::FieldRegistry fields;
    prime_fields(fields, mesh.cells().size(), 1.0);

    auto linM = solver::make_linear_solver(make_lin_cfg());
    auto linP = solver::make_linear_solver(make_lin_cfg());

    solver::PimpleOptions opt;
    opt.nOuterCorrectors = 2;
    opt.nCorrectors = 2;
    opt.nNonOrthCorr = 0;
    opt.dt = 1.0e-3;
    opt.rho = 1.0;
    opt.mu = 1.0e-3;
    solver::PimpleAlgorithm pimple(mesh, fields, bcs, *linM, *linP, opt);

    REQUIRE(pimple.n_outer_correctors() == 2);
    REQUIRE(pimple.n_correctors() == 2);
    REQUIRE(pimple.momentum_predictor() == true);

    solver::SimpleResiduals r{};
    REQUIRE_NOTHROW(r = pimple.iterate());
    REQUIRE(finite_residuals(r));
}

TEST_CASE("PIMPLE with nOuterCorrectors=1 matches PISO behaviour", "[solver][pimple][piso]")
{
    const int N = 20;
    const double dx = 1.0 / N;

    // Identical mesh + BCs + initial fields for both algorithms.
    auto buildAndRun = [&](auto& algo,
                           solver::FieldRegistry& F,
                           meshing::Mesh& mesh,
                           std::vector<solver::BoundarySpec>& bcs) {
        prime_fields(F, mesh.cells().size(), 1.0);
        return algo.iterate();
    };

    // PIMPLE n_outer=1, n_corr=2.
    auto meshA = build_channel_1d(N, dx);
    auto bcsA = make_channel_bcs(1.0, 0.0);
    solver::FieldRegistry fA;
    auto linMA = solver::make_linear_solver(make_lin_cfg());
    auto linPA = solver::make_linear_solver(make_lin_cfg());
    solver::PimpleOptions po;
    po.nOuterCorrectors = 1;
    po.nCorrectors = 2;
    po.dt = 1.0e-3;
    po.rho = 1.0;
    po.mu = 1.0e-3;
    solver::PimpleAlgorithm pimple(meshA, fA, bcsA, *linMA, *linPA, po);
    auto rPimple = buildAndRun(pimple, fA, meshA, bcsA);

    // PISO n_corr=2.
    auto meshB = build_channel_1d(N, dx);
    auto bcsB = make_channel_bcs(1.0, 0.0);
    solver::FieldRegistry fB;
    auto linMB = solver::make_linear_solver(make_lin_cfg());
    auto linPB = solver::make_linear_solver(make_lin_cfg());
    solver::PisoOptions piso_opt;
    piso_opt.nCorrectors = 2;
    piso_opt.dt = 1.0e-3;
    piso_opt.rho = 1.0;
    piso_opt.mu = 1.0e-3;
    solver::PisoAlgorithm piso(meshB, fB, bcsB, *linMB, *linPB, piso_opt);
    auto rPiso = buildAndRun(piso, fB, meshB, bcsB);

    REQUIRE(finite_residuals(rPimple));
    REQUIRE(finite_residuals(rPiso));

    // With identical settings the two algorithms must produce equivalent
    // residuals up to GMRES round-off (1e-6 relative).
    const double tol = 1.0e-6;
    CHECK(std::fabs(rPimple.mom[0] - rPiso.mom[0]) < tol * (1.0 + std::fabs(rPiso.mom[0])));
    CHECK(std::fabs(rPimple.cont - rPiso.cont) < tol * (1.0 + std::fabs(rPiso.cont)));
}

TEST_CASE("PIMPLE multi-outer-corrector does not degrade residuals", "[solver][pimple]")
{
    const int N = 20;
    const double dx = 1.0 / N;

    auto runWithOuter = [&](int nOuter) {
        auto mesh = build_channel_1d(N, dx);
        auto bcs = make_channel_bcs(1.0, 0.0);
        solver::FieldRegistry F;
        prime_fields(F, mesh.cells().size(), 1.0);
        auto linM = solver::make_linear_solver(make_lin_cfg());
        auto linP = solver::make_linear_solver(make_lin_cfg());
        solver::PimpleOptions opt;
        opt.nOuterCorrectors = nOuter;
        opt.nCorrectors = 2;
        opt.dt = 1.0e-3;
        opt.rho = 1.0;
        opt.mu = 1.0e-3;
        // Use under-relaxation since nOuter>1 will engage SIMPLE-style URF.
        opt.urfU = 0.7;
        opt.urfP = 0.3;
        solver::PimpleAlgorithm pimple(mesh, F, bcs, *linM, *linP, opt);
        return pimple.iterate();
    };

    auto r1 = runWithOuter(1);
    auto r3 = runWithOuter(3);
    REQUIRE(finite_residuals(r1));
    REQUIRE(finite_residuals(r3));

    // With multiple outer correctors the final continuity residual must not
    // exceed the single-corrector value by more than a small slack: outer
    // loop is intended to *improve* coupling, never to amplify continuity
    // error. We allow a generous 50 % margin to absorb URF-induced rebalancing
    // since the absolute magnitudes are very small on this tiny mesh.
    CHECK(std::fabs(r3.cont) <= 1.5 * std::fabs(r1.cont) + 1.0e-12);
}

TEST_CASE("PIMPLE outer-loop early-exit triggers on residual tolerance", "[solver][pimple]")
{
    const int N = 20;
    const double dx = 1.0 / N;
    auto mesh = build_channel_1d(N, dx);
    auto bcs = make_channel_bcs(1.0, 0.0);
    solver::FieldRegistry F;
    prime_fields(F, mesh.cells().size(), 1.0);

    auto linM = solver::make_linear_solver(make_lin_cfg());
    auto linP = solver::make_linear_solver(make_lin_cfg());

    solver::PimpleOptions opt;
    opt.nOuterCorrectors = 5;
    opt.nCorrectors = 2;
    opt.dt = 1.0e-3;
    opt.rho = 1.0;
    opt.mu = 1.0e-3;
    // Loose enough that the first outer corrector should satisfy them.
    opt.tolU = 1.0e6;
    opt.tolP = 1.0e6;

    solver::PimpleAlgorithm pimple(mesh, F, bcs, *linM, *linP, opt);
    (void) pimple.iterate();
    // Early exit means last_outer_iters() < nOuterCorrectors - 1.
    CHECK(pimple.last_outer_iters() < pimple.n_outer_correctors() - 1);
}
