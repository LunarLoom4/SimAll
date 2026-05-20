// =============================================================================
// SimAll Beta — Solver Unit Tests
// File   : tests/unit/solver/test_tfqmr_solver.cpp
//
// Validates transpose-free QMR on:
//   1) SPD 1-D Laplacian (TFQMR should still converge).
//   2) Non-symmetric convection-diffusion (this is TFQMR's home turf).
// =============================================================================
#include "solver/LinearSolvers.hpp"
#include "solver/Solver.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

using namespace simall::solver;
using util::aligned_vector;

namespace
{
CSRMatrix make_laplacian(int n)
{
    CSRMatrix A;
    A.rowPtr.push_back(0);
    for (int i = 0; i < n; ++i) {
        if (i > 0) {
            A.colIdx.push_back(i - 1);
            A.values.push_back(-1.0);
        }
        A.colIdx.push_back(i);
        A.values.push_back(2.0);
        if (i < n - 1) {
            A.colIdx.push_back(i + 1);
            A.values.push_back(-1.0);
        }
        A.rowPtr.push_back(static_cast<int>(A.colIdx.size()));
    }
    return A;
}

// Non-symmetric 1-D convection-diffusion: -u_xx + Pe u_x = b
// Discretised with upwind on the convective term, Pe = 5.
CSRMatrix make_conv_diff(int n)
{
    CSRMatrix A;
    const double Pe = 5.0;
    A.rowPtr.push_back(0);
    for (int i = 0; i < n; ++i) {
        if (i > 0) {
            A.colIdx.push_back(i - 1);
            A.values.push_back(-1.0 - Pe);
        }
        A.colIdx.push_back(i);
        A.values.push_back(2.0 + Pe);
        if (i < n - 1) {
            A.colIdx.push_back(i + 1);
            A.values.push_back(-1.0);
        }
        A.rowPtr.push_back(static_cast<int>(A.colIdx.size()));
    }
    return A;
}
} // namespace

TEST_CASE("TFQMR converges on SPD Laplacian", "[solver][tfqmr]")
{
    const int n = 30;
    auto A = make_laplacian(n);
    aligned_vector<double> b(n, 1.0), x(n, 0.0);

    LinearSolverConfig cfg;
    cfg.kind = LinearSolverKind::TFQMR;
    cfg.preconditioner = PreconditionerKind::Jacobi;
    cfg.maxIterations = 500;
    cfg.tolerance = 1.0e-8;
    auto s = make_linear_solver(cfg);
    s->solve(A, b, x);

    aligned_vector<double> Ax(n, 0);
    A.spmv(x, Ax);
    double res = 0;
    for (int i = 0; i < n; ++i) {
        const double d = Ax[i] - b[i];
        res += d * d;
    }
    REQUIRE(std::sqrt(res) < 1.0e-5);
}

TEST_CASE("TFQMR converges on non-symmetric convection-diffusion", "[solver][tfqmr][nonsym]")
{
    const int n = 40;
    auto A = make_conv_diff(n);
    aligned_vector<double> b(n, 1.0), x(n, 0.0);

    LinearSolverConfig cfg;
    cfg.kind = LinearSolverKind::TFQMR;
    cfg.preconditioner = PreconditionerKind::ILU;
    cfg.maxIterations = 500;
    cfg.tolerance = 1.0e-8;
    auto s = make_linear_solver(cfg);
    s->solve(A, b, x);

    aligned_vector<double> Ax(n, 0);
    A.spmv(x, Ax);
    double res = 0;
    for (int i = 0; i < n; ++i) {
        const double d = Ax[i] - b[i];
        res += d * d;
    }
    REQUIRE(std::sqrt(res) < 1.0e-4);
}
