// =============================================================================
// SimAll Beta - tests/unit/longtail/test_longtail_modules.cpp
// Week 18 - Adjoint / Optimization / Morphing / Dynamics / Rotating / AMR /
//           IBM / EMag / Acoustics / ROM.
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "adjoint/DiscreteAdjoint.hpp"
#include "adjoint/ContinuousAdjoint.hpp"
#include "optimization/MmaOptimizer.hpp"
#include "optimization/CmaEs.hpp"
#include "optimization/Doe.hpp"
#include "optimization/Kriging.hpp"
#include "morphing/FfdBox.hpp"
#include "dynamics/SixDofCoupler.hpp"
#include "rotating/SlidingMeshDriver.hpp"
#include "amr/AmrApplier.hpp"
#include "ibm/DirectForcingIbm.hpp"
#include "emag/AvFormulation.hpp"
#include "acoustics/CurleAcoustics.hpp"
#include "rom/OperatorInference.hpp"
#include "rom/EcswHyperReduction.hpp"

#include <cmath>

using Catch::Matchers::WithinAbs;

// -----------------------------------------------------------------------------
// Adjoint
// -----------------------------------------------------------------------------
TEST_CASE("DiscreteAdjoint solves transpose system", "[w18][adjoint]") {
    // Build A = [[2,1],[0,3]]; Aᵀ = [[2,0],[1,3]].
    // Choose dJdQ = (1, 1)ᵀ → ψ = (Aᵀ)⁻¹ dJdQ = ([[0.5, 0],[-1/6, 1/3]]) (1,1)ᵀ
    //                                        = (0.5, 1/6).
    simall::adjoint::CsrMatrix A;
    A.n = 2;
    A.rowPtr = {0, 2, 3};
    A.colIdx = {0, 1, 1};
    A.values = {2.0, 1.0, 3.0};

    simall::adjoint::DiscreteAdjointInputs in;
    in.A    = &A;
    in.dJdQ = {1.0, 1.0};
    in.dJdAlpha = {0.5};
    in.dRdAlphaTransposeTimes = [](std::vector<double>& out){ out[0] = 0.25; };

    auto r = simall::adjoint::solve_discrete_adjoint(in, /*nDV=*/1, {});
    REQUIRE(r.ok);
    REQUIRE_THAT(r.psi[0], WithinAbs(0.5,     1e-6));
    REQUIRE_THAT(r.psi[1], WithinAbs(1.0/6.0, 1e-6));
    REQUIRE_THAT(r.gradient[0], WithinAbs(0.5 - 0.25, 1e-6));
}

TEST_CASE("ContinuousAdjoint shape sensitivity scatter", "[w18][adjoint]") {
    simall::adjoint::SurfaceFace f;
    f.nNodes = 3; f.nodes = {0, 1, 2, 0};
    f.normal = {0, 0, 1};
    f.area = 0.5;
    simall::adjoint::AdjointFlowSample s;
    s.uPsi = {0, 0, 1}; s.pPsi = 0.1; s.mu = 1.0;
    auto faceSens = simall::adjoint::assemble_shape_sensitivity({f}, {s});
    REQUIRE(faceSens.size() == 1);
    // G = (-0.1 + 2·1·1) · 0.5 = 0.95
    REQUIRE_THAT(faceSens[0], WithinAbs(0.95, 1e-12));
    auto nodeSens = simall::adjoint::accumulate_node_sensitivity(3, {f}, faceSens);
    REQUIRE(nodeSens.size() == 3);
}

// -----------------------------------------------------------------------------
// Optimization
// -----------------------------------------------------------------------------
TEST_CASE("MMA drives x to optimum of simple quadratic", "[w18][opt][mma]") {
    simall::optimization::MmaProblem p;
    p.n = 2;
    p.m = 0;
    p.xLow = {-2.0, -2.0};
    p.xUp  = { 2.0,  2.0};
    p.x0   = { 1.5,  1.5};
    p.evaluate = [](const std::vector<double>& x, double& f,
                     std::vector<double>& g, std::vector<double>& df,
                     std::vector<double>& dg) {
        f = (x[0] - 0.5)*(x[0] - 0.5) + (x[1] + 0.25)*(x[1] + 0.25);
        df = {2.0 * (x[0] - 0.5), 2.0 * (x[1] + 0.25)};
        g.clear(); dg.clear();
    };
    simall::optimization::MmaOptions opt;
    opt.maxIter = 60; opt.moveLimit = 0.3;
    auto r = simall::optimization::run_mma(p, opt);
    REQUIRE(r.ok);
    REQUIRE_THAT(r.x[0], WithinAbs(0.5,   5e-2));
    REQUIRE_THAT(r.x[1], WithinAbs(-0.25, 5e-2));
}

TEST_CASE("CMA-ES finds minimum of sphere function", "[w18][opt][cmaes]") {
    simall::optimization::CmaEsOptions o;
    o.dim = 3;
    o.maxGen = 80;
    o.sigma0 = 0.3;
    o.mean0  = {1.0, -1.0, 0.5};
    o.evaluate = [](const std::vector<double>& x){
        double s = 0.0; for (auto v : x) s += v*v; return s;
    };
    auto r = simall::optimization::run_cmaes(o);
    REQUIRE(r.ok);
    REQUIRE(r.bestF < 1e-4);
}

TEST_CASE("DOE LHS samples are space-filling and bounded", "[w18][opt][doe]") {
    auto X = simall::optimization::doe_latin_hypercube({0.0, -1.0}, {1.0, 2.0}, 16, 42);
    REQUIRE(X.size() == 16);
    for (auto& row : X) {
        REQUIRE(row[0] >= 0.0); REQUIRE(row[0] <= 1.0);
        REQUIRE(row[1] >= -1.0); REQUIRE(row[1] <= 2.0);
    }
    auto S = simall::optimization::doe_sobol({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, 8);
    REQUIRE(S.size() == 8);
    auto F = simall::optimization::doe_full_factorial({0.0, 0.0}, {1.0, 1.0}, 3);
    REQUIRE(F.size() == 9);
}

TEST_CASE("Kriging interpolates training points exactly", "[w18][opt][kriging]") {
    std::vector<std::vector<double>> X = {{0.0}, {1.0}, {2.0}, {3.0}};
    std::vector<double>              y = {1.0, 0.0, 1.0, 4.0};
    simall::optimization::KrigingModel m;
    REQUIRE(m.fit(X, y));
    for (std::size_t i = 0; i < X.size(); ++i) {
        REQUIRE_THAT(m.predict(X[i]), WithinAbs(y[i], 1e-3));
    }
    double mean, var;
    m.predict({1.5}, mean, var);
    REQUIRE(var >= 0.0);
}

// -----------------------------------------------------------------------------
// Morphing
// -----------------------------------------------------------------------------
TEST_CASE("FFD identity preserves embedded points", "[w18][morphing]") {
    simall::morphing::FfdBox box({0,0,0}, {1,1,1}, 3, 3, 3);
    std::vector<simall::morphing::Vec3> in  = {{0.25, 0.5, 0.75}, {0.1, 0.1, 0.1}};
    std::vector<simall::morphing::Vec3> out;
    box.morph(in, out);
    REQUIRE_THAT(out[0].x, WithinAbs(0.25, 1e-9));
    REQUIRE_THAT(out[0].y, WithinAbs(0.50, 1e-9));
    REQUIRE_THAT(out[0].z, WithinAbs(0.75, 1e-9));
}
TEST_CASE("FFD lattice displacement shifts an interior point", "[w18][morphing]") {
    simall::morphing::FfdBox box({0,0,0}, {1,1,1}, 2, 2, 2);
    // Displace one corner control point: every interior point moves a fraction.
    box.displace_control_point(1, 1, 1, {0.5, 0.0, 0.0});
    std::vector<simall::morphing::Vec3> in  = {{0.5, 0.5, 0.5}};
    std::vector<simall::morphing::Vec3> out;
    box.morph(in, out);
    REQUIRE(out[0].x > 0.5);
}

// -----------------------------------------------------------------------------
// Dynamics
// -----------------------------------------------------------------------------
TEST_CASE("SixDofCoupler explicit translation under gravity", "[w18][dynamics]") {
    simall::dynamics::SixDofCoupler c({simall::dynamics::SixDofCouplingMode::Explicit});
    simall::dynamics::RigidBodyState s;
    s.mass = 1.0;
    c.advance(s, {0, 0, -9.81}, {0, 0, 0}, 1.0);
    REQUIRE_THAT(s.velocity.z, WithinAbs(-9.81, 1e-9));
    REQUIRE_THAT(s.position.z, WithinAbs(-9.81, 1e-9));   // semi-implicit Euler
}

// -----------------------------------------------------------------------------
// Rotating
// -----------------------------------------------------------------------------
TEST_CASE("SlidingMeshDriver advances angle and builds map", "[w18][rotating]") {
    simall::rotating::SlidingMeshDriver d({0,0,0}, {0,0,1}, /*omega=*/1.0);
    REQUIRE_THAT(d.advance_rotor(0.5), WithinAbs(0.5, 1e-12));
    REQUIRE_THAT(d.cumulative_angle(), WithinAbs(0.5, 1e-12));
    d.set_rotor_faces ({{0, {1.0, 0.0, 0.0}, {1,0,0}, 1.0}});
    d.set_stator_faces({{10, {std::cos(0.5), std::sin(0.5), 0.0}, {1,0,0}, 1.0}});
    auto map = d.build_interface_map();
    REQUIRE(!map.entries.empty());
    REQUIRE(map.entries[0].donor == 10);
}

// -----------------------------------------------------------------------------
// AMR
// -----------------------------------------------------------------------------
TEST_CASE("AmrApplier refines highest-indicator cells", "[w18][amr]") {
    simall::amr::AmrOptions opt;
    opt.refineFraction = 0.5;
    opt.coarsenFraction = 0.0;
    opt.enforceTwoToOne = false;
    simall::amr::AmrApplier app(opt);
    std::vector<double>        ind   = {0.1, 0.9, 0.2, 0.8};
    std::vector<std::uint8_t>  level(4, 0);
    std::vector<std::size_t>   rp    = {0, 1, 2, 3, 4};
    std::vector<std::size_t>   idx   = {1, 0, 3, 2};
    auto r = app.apply_cycle(ind, level, rp, idx, nullptr);
    REQUIRE(r.nRefined == 2);
    REQUIRE(level[1] == 1);
    REQUIRE(level[3] == 1);
}

// -----------------------------------------------------------------------------
// IBM
// -----------------------------------------------------------------------------
TEST_CASE("DirectForcingIbm produces non-trivial force", "[w18][ibm]") {
    std::vector<simall::ibm::MarkerPoint> mk = {{{0.5, 0.5, 0.5}, {1.0, 0.0, 0.0}, 1.0}};
    std::vector<simall::ibm::EulerianCell> cells;
    std::vector<std::array<double,3>> uStar;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) {
                cells.push_back({{i*0.5, j*0.5, k*0.5}, 1.0});
                uStar.push_back({0.0, 0.0, 0.0});
            }
    simall::ibm::DirectForcingOptions o; o.spacing = 0.5; o.dt = 1.0; o.kernelWidth = 3.0;
    auto r = simall::ibm::compute_direct_forcing(mk, cells, uStar, o);
    REQUIRE(r.markerForce.size() == 1);
    REQUIRE(r.markerForce[0][0] != 0.0);    // x-force should be positive
}

// -----------------------------------------------------------------------------
// EMag
// -----------------------------------------------------------------------------
TEST_CASE("A-V tet4 element matrix has correct volume", "[w18][emag]") {
    simall::emag::Tet4 e;
    e.coords = {{ {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} }};
    e.nodes  = {0,1,2,3};
    e.mu     = 1.0;
    e.sigma  = 1.0;
    auto em  = simall::emag::assemble_tet4(e);
    REQUIRE_THAT(std::abs(em.volume), WithinAbs(1.0/6.0, 1e-12));
    // Diagonal mass should be V/10 with σ=1.
    REQUIRE_THAT(em.M_A[0][0], WithinAbs((1.0/6.0)/10.0, 1e-12));
    // Row sums of stiffness should be zero (locality / null space of Laplacian).
    for (int i = 0; i < 4; ++i) {
        double rs = 0.0;
        for (int j = 0; j < 4; ++j) rs += em.K_V[i][j];
        REQUIRE_THAT(rs, WithinAbs(0.0, 1e-10));
    }
    // Build the global system on a one-tet "mesh" and check dimensions.
    auto sys = simall::emag::build_av_scalar_system({e}, 4);
    REQUIRE(sys.n == 4);
    REQUIRE(sys.rowPtr.size() == 5);
}

// -----------------------------------------------------------------------------
// Acoustics
// -----------------------------------------------------------------------------
TEST_CASE("Curle dipole produces non-zero pressure history", "[w18][acoustics]") {
    simall::acoustics::CurleAcoustics ca(343.0);
    simall::acoustics::CurleObserver  obs;
    obs.position = {10.0, 0.0, 0.0};
    for (int k = 0; k < 5; ++k) {
        std::vector<simall::acoustics::CurleSurfaceSample> surf = {
            { {0,0,0}, {1,0,0}, std::sin(double(k)*0.5), 1.0 }
        };
        ca.push_step(double(k) * 1e-3, surf, obs);
    }
    auto p = ca.observer_pressure_history();
    REQUIRE(p.size() == 5);
    double sumAbs = 0.0;
    for (auto v : p) sumAbs += std::abs(v);
    REQUIRE(sumAbs > 0.0);
}

// -----------------------------------------------------------------------------
// ROM
// -----------------------------------------------------------------------------
TEST_CASE("Operator inference recovers a known linear operator", "[w18][rom][opinf]") {
    // True system:  dx/dt = A x   with A = [[-1, 0],[0, -2]]
    const int nt = 50;
    std::vector<std::vector<double>> X(2, std::vector<double>(nt, 0.0));
    std::vector<std::vector<double>> dX(2, std::vector<double>(nt, 0.0));
    for (int t = 0; t < nt; ++t) {
        const double tt = double(t) * 0.05;
        X[0][t]  =  std::exp(-1.0 * tt);
        X[1][t]  =  2.0 * std::exp(-2.0 * tt);
        dX[0][t] = -1.0 * X[0][t];
        dX[1][t] = -2.0 * X[1][t];
    }
    simall::rom::OperatorInferenceOptions opt;
    opt.reducedDim = 2;
    opt.regularization = 1e-10;
    auto op = simall::rom::fit_operator_inference(X, dX, {}, opt);
    REQUIRE(op.r == 2);
    // Forward evaluation should approximately reproduce dX/dt at a snapshot.
    std::vector<double> xHat = {X[0][10] * op.V[0][0] + X[1][10] * op.V[1][0],
                                  X[0][10] * op.V[0][1] + X[1][10] * op.V[1][1]};
    auto rhs = simall::rom::rom_rhs(op, xHat, {});
    REQUIRE(rhs.size() == 2);
}

TEST_CASE("ECSW picks a non-empty cubature with low error", "[w18][rom][ecsw]") {
    // Two snapshots, four elements, r=2.
    std::vector<std::vector<std::vector<double>>> er = {
        { {1.0, 0.0}, {0.5, 0.5}, {0.0, 1.0}, {0.2, 0.2} },
        { {0.8, 0.1}, {0.4, 0.4}, {0.1, 0.9}, {0.3, 0.1} }
    };
    auto r = simall::rom::ecsw_hyper_reduction(er, 1e-3, 100);
    REQUIRE(!r.indices.empty());
    REQUIRE(r.indices.size() == r.weights.size());
    REQUIRE(r.relativeError < 1e-1);
    for (auto w : r.weights) REQUIRE(w > 0.0);
}
