// =============================================================================
// SimAll Beta - GPU subsystem unit tests (Week 13).
//
// Exercises the SIMALL_NO_CUDA fall-back paths of every W13 module so the
// suite can run on CI nodes without a GPU.  Where CUDA is present, the
// same tests still pass (kernels delegate to cuBLAS / cuSPARSE / on-device
// implementations).
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "gpu/CudaContext.hpp"
#include "gpu/DeviceMemoryPool.hpp"
#include "gpu/StreamScheduler.hpp"
#include "gpu/HostDeviceMirror.hpp"
#include "gpu/GpuKernels.hpp"
#include "gpu/CgGpu.hpp"
#include "gpu/BicgstabGpu.hpp"
#include "gpu/RenderingBridge.hpp"

#include <cmath>
#include <numeric>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace sg = simall::gpu;

// -----------------------------------------------------------------------------
TEST_CASE("CudaContext singleton is always available", "[gpu][context]") {
    auto& ctx = sg::CudaContext::current();
    // Either a real device or the serial fall-back.
    REQUIRE(ctx.device_count() >= 0);
    REQUIRE_NOTHROW(ctx.synchronize());
    REQUIRE(!ctx.last_error().empty());
    auto p = ctx.properties();
    REQUIRE(p.warp_size >= 1);
}

// -----------------------------------------------------------------------------
TEST_CASE("DeviceMemoryPool bucket math + cache hit/miss", "[gpu][pool]") {
    sg::DeviceMemoryPool pool;
    // 100 B  →  bucket k=8 (256 B), 300 B → k=9 (512 B), 1025 B → k=11 (2 KiB)
    REQUIRE(sg::DeviceMemoryPool::bucket_for(1)      == 8);
    REQUIRE(sg::DeviceMemoryPool::bucket_for(100)    == 8);
    REQUIRE(sg::DeviceMemoryPool::bucket_for(300)    == 9);
    REQUIRE(sg::DeviceMemoryPool::bucket_for(1025)   == 11);

    void* a = pool.allocate(300);
    REQUIRE(a != nullptr);
    auto s1 = pool.stats();
    REQUIRE(s1.allocation_calls == 1);
    REQUIRE(s1.freelist_misses  == 1);
    REQUIRE(s1.total_in_use_bytes == 512);

    pool.release(a, 300);
    auto s2 = pool.stats();
    REQUIRE(s2.total_in_use_bytes      == 0);
    REQUIRE(s2.total_in_freelist_bytes == 512);

    void* b = pool.allocate(300);
    auto s3 = pool.stats();
    REQUIRE(s3.freelist_hits == 1);
    REQUIRE(b == a);   // LIFO cache hit
    pool.release(b, 300);
    pool.purge();
    auto s4 = pool.stats();
    REQUIRE(s4.total_in_freelist_bytes == 0);
    REQUIRE(s4.total_reserved_bytes    == 0);
}

// -----------------------------------------------------------------------------
TEST_CASE("StreamScheduler acquires and releases streams", "[gpu][stream]") {
    sg::StreamScheduler sched(3);
    REQUIRE(sched.stream_count() == 3);

    void* s1 = sched.acquire();
    void* s2 = sched.acquire();
    REQUIRE(s1 != nullptr);
    REQUIRE(s2 != nullptr);
    // record/wait/release form a closed cycle; must not crash.
    void* e = sched.record_event(s1);
    REQUIRE(e != nullptr);
    sched.wait_event(s2, e);
    sched.release_event(e);
    sched.release(s1);
    sched.release(s2);
    sched.synchronize_all();
}

// -----------------------------------------------------------------------------
TEST_CASE("HostDeviceMirror RAII + sync", "[gpu][mirror]") {
    sg::HostDeviceMirror<double> m(8);
    REQUIRE(m.size() == 8);
    REQUIRE(m.size_in_bytes() == 8 * sizeof(double));
    REQUIRE(m.host_ptr() != nullptr);
    REQUIRE(m.device_ptr() != nullptr);

    std::iota(m.host_begin(), m.host_end(), 1.0);
    m.to_device();
    m.zero_device();    // wipes device side
    m.to_host();
    for (double v : m.host_view()) REQUIRE_THAT(v, WithinAbs(0.0, 1e-15));

    // Move constructor must not double-free.
    sg::HostDeviceMirror<double> m2 = std::move(m);
    REQUIRE(m2.size() == 8);
}

// -----------------------------------------------------------------------------
TEST_CASE("GpuKernels: axpy / dot / nrm2 / scal", "[gpu][kernels]") {
    constexpr std::size_t n = 5;
    sg::HostDeviceMirror<double> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
        x.host_ptr()[i] = double(i + 1);     // 1..5
        y.host_ptr()[i] = 1.0;
    }
    x.to_device(); y.to_device();

    sg::axpy(2.0, x.device_ptr(), y.device_ptr(), n);
    y.to_host();
    for (std::size_t i = 0; i < n; ++i)
        REQUIRE_THAT(y.host_ptr()[i], WithinAbs(1.0 + 2.0 * (i + 1), 1e-12));

    sg::scal(0.5, y.device_ptr(), n);
    y.to_host();
    REQUIRE_THAT(y.host_ptr()[0], WithinAbs(1.5, 1e-12));    // (1+2*1)/2

    const double d = sg::dot(x.device_ptr(), x.device_ptr(), n);
    REQUIRE_THAT(d, WithinAbs(55.0, 1e-12));                  // 1+4+9+16+25
    const double r = sg::nrm2(x.device_ptr(), n);
    REQUIRE_THAT(r, WithinRel(std::sqrt(55.0), 1e-12));
}

// -----------------------------------------------------------------------------
TEST_CASE("GpuKernels: CSR SpMV reference", "[gpu][kernels]") {
    // 4x4 diagonal-dominant SPD matrix
    //   4 -1  0  0
    //  -1  4 -1  0
    //   0 -1  4 -1
    //   0  0 -1  4
    std::vector<int>    rp = {0, 2, 5, 8, 10};
    std::vector<int>    ci = {0,1, 0,1,2, 1,2,3, 2,3};
    std::vector<double> vs = {4,-1, -1,4,-1, -1,4,-1, -1,4};
    sg::HostDeviceMirror<int>    dRow(rp.size()), dCol(ci.size());
    sg::HostDeviceMirror<double> dVal(vs.size()), dX(4), dY(4);
    std::copy(rp.begin(), rp.end(), dRow.host_ptr());
    std::copy(ci.begin(), ci.end(), dCol.host_ptr());
    std::copy(vs.begin(), vs.end(), dVal.host_ptr());
    for (int i = 0; i < 4; ++i) dX.host_ptr()[i] = double(i + 1);  // 1..4
    dRow.to_device(); dCol.to_device(); dVal.to_device(); dX.to_device();

    sg::spmv_csr(4, dRow.device_ptr(), dCol.device_ptr(), dVal.device_ptr(),
                 dX.device_ptr(), dY.device_ptr());
    dY.to_host();
    // y = A x  →  [4-2, -1+8-3, -2+12-4, -3+16] = [2,4,6,13]
    REQUIRE_THAT(dY.host_ptr()[0], WithinAbs(2.0,  1e-12));
    REQUIRE_THAT(dY.host_ptr()[1], WithinAbs(4.0,  1e-12));
    REQUIRE_THAT(dY.host_ptr()[2], WithinAbs(6.0,  1e-12));
    REQUIRE_THAT(dY.host_ptr()[3], WithinAbs(13.0, 1e-12));
}

// -----------------------------------------------------------------------------
TEST_CASE("CgGpu solves tridiagonal SPD system", "[gpu][cg]") {
    // n=10, A = tridiag(-1, 4, -1), b = e_1+e_n.
    constexpr int n = 10;
    std::vector<int>    rp;  rp.reserve(n + 1);
    std::vector<int>    ci;
    std::vector<double> vs;
    rp.push_back(0);
    for (int i = 0; i < n; ++i) {
        if (i > 0)     { ci.push_back(i - 1); vs.push_back(-1.0); }
        ci.push_back(i);                        vs.push_back( 4.0);
        if (i + 1 < n) { ci.push_back(i + 1); vs.push_back(-1.0); }
        rp.push_back(static_cast<int>(ci.size()));
    }
    std::vector<double> b(n, 0.0), x(n, 0.0);
    b.front() = 1.0; b.back() = 1.0;

    sg::HostCsrView view{static_cast<std::size_t>(n),
                         static_cast<std::size_t>(ci.size()),
                         rp.data(), ci.data(), vs.data()};

    sg::GpuSolverConfig cfg;
    cfg.preconditioner = sg::GpuPreconditioner::Jacobi;
    cfg.tolerance      = 1.0e-10;
    cfg.maxIterations  = 200;
    sg::CgGpu cg(cfg);
    cg.set_matrix(view);
    auto st = cg.solve(b, x);
    REQUIRE(st.converged);
    REQUIRE(st.iterations < 50);

    // Verify residual ||b - A x|| ≤ tol*||b||.
    std::vector<double> r(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = rp[i]; k < rp[i + 1]; ++k) s += vs[k] * x[ci[k]];
        r[i] = b[i] - s;
    }
    double rn = 0.0;
    for (double v : r) rn += v*v;
    REQUIRE(std::sqrt(rn) < 1e-8);
}

// -----------------------------------------------------------------------------
TEST_CASE("BicgstabGpu solves non-symmetric system", "[gpu][bicg]") {
    // n=8, A = upper-bidiagonal with diag 3 and super-diag 1, lower -0.5
    constexpr int n = 8;
    std::vector<int>    rp{0};
    std::vector<int>    ci;
    std::vector<double> vs;
    for (int i = 0; i < n; ++i) {
        if (i > 0)     { ci.push_back(i - 1); vs.push_back(-0.5); }
        ci.push_back(i);                        vs.push_back( 3.0);
        if (i + 1 < n) { ci.push_back(i + 1); vs.push_back( 1.0); }
        rp.push_back(static_cast<int>(ci.size()));
    }
    std::vector<double> b(n, 1.0), x(n, 0.0);

    sg::HostCsrView view{static_cast<std::size_t>(n),
                         static_cast<std::size_t>(ci.size()),
                         rp.data(), ci.data(), vs.data()};

    sg::GpuSolverConfig cfg;
    cfg.preconditioner = sg::GpuPreconditioner::ILU0;
    cfg.tolerance      = 1.0e-10;
    cfg.maxIterations  = 200;
    sg::BicgstabGpu bicg(cfg);
    bicg.set_matrix(view);
    auto st = bicg.solve(b, x);
    REQUIRE(st.converged);

    // Residual check.
    std::vector<double> r(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = rp[i]; k < rp[i + 1]; ++k) s += vs[k] * x[ci[k]];
        r[i] = b[i] - s;
    }
    double rn = 0.0;
    for (double v : r) rn += v*v;
    REQUIRE(std::sqrt(rn) < 1e-7);
}

// -----------------------------------------------------------------------------
TEST_CASE("RenderingBridge registers / acquires / releases", "[gpu][bridge]") {
    auto& br = sg::RenderingBridge::instance();
    sg::HostDeviceMirror<float> velocityMagnitude(64);
    for (std::size_t i = 0; i < velocityMagnitude.size(); ++i)
        velocityMagnitude.host_ptr()[i] = float(i);
    velocityMagnitude.to_device();

    auto id = br.register_scalar("|U|", &velocityMagnitude);
    REQUIRE(id != sg::kInvalidBufferId);
    void* p = br.acquire(id);
    REQUIRE(p != nullptr);
    br.release(id);

    auto listing = br.list_buffers();
    bool found = false;
    for (const auto& d : listing) if (d.id == id) { found = true; break; }
    REQUIRE(found);

    br.unregister(id);
}
