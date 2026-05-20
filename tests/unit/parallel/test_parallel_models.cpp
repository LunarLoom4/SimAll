// =============================================================================
// SimAll Beta — Parallel/HPC Unit Tests (Week 12)
// File   : tests/unit/parallel/test_parallel_models.cpp
//
// Sanity tests for the six W12 parallel/HPC modules.  All tests run on
// a single rank in the no-MPI build (collective ops are pass-through).
//   1. MpiContext            (rank/size, serial allreduce identity)
//   2. GhostExchange         (empty-neighbour no-op, layout query)
//   3. DomainPartition       (single-rank → all cells local, 0 ghosts)
//   4. FieldReducer          (sum / L2 / min / max / mean)
//   5. LoadBalancer          (balanced samples → no repart)
//   6. NumaPinning           (discover non-zero core count)
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "parallel/MpiContext.hpp"
#include "parallel/GhostExchange.hpp"
#include "parallel/DomainPartition.hpp"
#include "parallel/FieldReducer.hpp"
#include "parallel/LoadBalancer.hpp"
#include "parallel/NumaPinning.hpp"
#include "meshing/CartesianMesher.hpp"
#include "meshing/MeshStorage.hpp"

using namespace simall;
using Catch::Approx;

TEST_CASE("MpiContext::world has valid rank/size and serial collectives", "[parallel][mpi]") {
    auto& ctx = parallel::MpiContext::world();
    REQUIRE(ctx.size() >= 1);
    REQUIRE(ctx.rank() >= 0);
    REQUIRE(ctx.allreduce_sum(3.5) == Approx(3.5 * ctx.size()));
    REQUIRE(ctx.allreduce_max(7.0) == Approx(7.0));
    ctx.barrier();
}

TEST_CASE("GhostExchange with no neighbours is a no-op", "[parallel][ghost]") {
    auto& ctx = parallel::MpiContext::world();
    parallel::GhostExchange gx(ctx);
    gx.set_layout(/*nOwned=*/10, /*nGhost=*/4);
    REQUIRE(gx.owned_count() == 10);
    REQUIRE(gx.ghost_count() == 4);
    REQUIRE(gx.neighbour_count() == 0);
    util::aligned_vector<double> v(14, 1.0);
    gx.sync_scalar(v);
    REQUIRE(v.size() == 14);
    REQUIRE(v[0] == Approx(1.0));
}

TEST_CASE("DomainPartition on one rank assigns all cells locally", "[parallel][partition]") {
    auto& ctx = parallel::MpiContext::world();
    meshing::CartesianGridSpec spec{}; spec.Nx = 4; spec.Ny = 2; spec.Nz = 2;
    spec.origin = {0,0,0}; spec.extent = {4,2,2};
    meshing::CartesianMesher cm(spec);
    meshing::Mesh mesh;
    cm.generate(mesh);

    parallel::DomainPartition dp(ctx);
    parallel::DomainPartitionProps props;
    props.nParts = ctx.size();
    dp.initialize(props);
    parallel::DomainPlan plan;
    const auto nLocal = dp.partition(mesh, plan);
    REQUIRE(nLocal > 0);
    REQUIRE(plan.cellRank.size() == mesh.cells().size());
    if (ctx.size() == 1) {
        REQUIRE(nLocal == mesh.cells().size());
        REQUIRE(plan.ghostCells.empty());
        REQUIRE(plan.neighbours.empty());
    }
    parallel::GhostExchange gx(ctx);
    dp.install_into(gx, plan);
    REQUIRE(gx.owned_count() == plan.localCells.size());
}

TEST_CASE("DomainPartition::extract_local_mesh wraps ops::extract_subdomain",
          "[parallel][partition][decomposition]") {
    auto& ctx = parallel::MpiContext::world();
    meshing::CartesianGridSpec spec{}; spec.Nx = 4; spec.Ny = 2; spec.Nz = 2;
    spec.origin = {0,0,0}; spec.extent = {4,2,2};
    meshing::CartesianMesher cm(spec);
    meshing::Mesh mesh;
    cm.generate(mesh);

    parallel::DomainPartition dp(ctx);
    parallel::DomainPartitionProps props;
    props.nParts = ctx.size();
    dp.initialize(props);
    parallel::DomainPlan plan;
    dp.partition(mesh, plan);

    meshing::Mesh local;
    auto stats = dp.extract_local_mesh(mesh, plan, local);

    REQUIRE(stats.ownedCells == plan.localCells.size());
    if (ctx.size() == 1) {
        REQUIRE(stats.ghostCells     == 0);
        REQUIRE(stats.interfaceFaces == 0);
        REQUIRE(local.cells().size() == mesh.cells().size());
    } else {
        REQUIRE(local.cells().size() ==
                stats.ownedCells + stats.ghostCells);
    }
}

TEST_CASE("FieldReducer scalar reductions match serial values", "[parallel][reducer]") {
    auto& ctx = parallel::MpiContext::world();
    parallel::FieldReducer red(ctx);
    util::aligned_vector<double> v{1.0, 2.0, 3.0, 4.0};
    const std::size_t n = v.size();
    REQUIRE(red.sum (v, n) == Approx(10.0 * ctx.size()));
    REQUIRE(red.l2  (v, n) == Approx(std::sqrt(30.0 * ctx.size())));
    REQUIRE(red.min (v, n) == Approx(1.0));
    REQUIRE(red.max (v, n) == Approx(4.0));
    REQUIRE(red.mean(v, n) == Approx(2.5));
}

TEST_CASE("LoadBalancer with balanced samples does not request repartition", "[parallel][balance]") {
    auto& ctx = parallel::MpiContext::world();
    parallel::LoadBalancer lb(ctx);
    parallel::LoadBalancerProps p; p.repartThreshold = 0.50;
    lb.initialize(p);
    for (int i = 0; i < 3; ++i) lb.record_iteration(1.0, 100);
    const auto stats = lb.evaluate();
    REQUIRE(stats.imbalance == Approx(0.0).margin(1e-9));
    REQUIRE_FALSE(stats.shouldRepart);
}

TEST_CASE("NumaPinning discovers at least one core", "[parallel][numa]") {
    const auto topo = parallel::NumaPinning::discover();
    REQUIRE(topo.totalCores >= 1);
    REQUIRE(topo.numaNodes >= 1);
    const auto map = parallel::NumaPinning::compute_thread_to_core_map(4);
    REQUIRE(map.size() == 4);
    for (auto c : map) REQUIRE(c >= 0);
    // pin_thread_to_core may legitimately return false on unsupported
    // platforms; just verify the API can be called without crashing.
    parallel::NumaPinning::pin_thread_to_core(0);
}
