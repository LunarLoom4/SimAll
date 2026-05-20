// =============================================================================
// SimAll Beta — Solver Unit Tests
// File   : tests/unit/solver/test_ssor_preconditioner.cpp
//
// Validates SSOR symmetric Gauss-Seidel preconditioner:
//   1) Applied to identity matrix returns input.
//   2) Reduces residual when fed into CG on an SPD 1-D Laplacian.
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
// Build 1-D Laplacian: tridiag(-1, 2, -1) of size n  (SPD).
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
} // namespace

TEST_CASE("SSOR preconditioner sanity on identity", "[solver][ssor]")
{
    CSRMatrix I;
    const int n = 5;
    I.rowPtr.push_back(0);
    for (int i = 0; i < n; ++i) {
        I.colIdx.push_back(i);
        I.values.push_back(1.0);
        I.rowPtr.push_back(static_cast<int>(I.colIdx.size()));
    }
    auto M = make_preconditioner(PreconditionerKind::SSOR);
    M->setup(I);
    aligned_vector<double> r{1, 2, 3, 4, 5}, z;
    M->apply(r, z);
    REQUIRE(z.size() == r.size());
    // For identity, SSOR(ω=1) returns (2-ω) · r = 1 · r  (the (2-ω) outer
    // scale cancels with the back-sweep). Within numerical tolerance:
    for (int i = 0; i < n; ++i)
        REQUIRE(z[i] == Catch::Approx(r[i]).margin(1e-12));
}

TEST_CASE("CG with SSOR converges on 1-D Laplacian", "[solver][ssor][cg]")
{
    const int n = 50;
    auto A = make_laplacian(n);
    aligned_vector<double> b(n, 1.0), x(n, 0.0);

    LinearSolverConfig cfg;
    cfg.kind = LinearSolverKind::CG;
    cfg.preconditioner = PreconditionerKind::SSOR;
    cfg.maxIterations = 200;
    cfg.tolerance = 1.0e-10;
    auto solver = make_linear_solver(cfg);
    const int iters = solver->solve(A, b, x);

    aligned_vector<double> Ax(n, 0);
    A.spmv(x, Ax);
    double res = 0;
    for (int i = 0; i < n; ++i) {
        const double d = Ax[i] - b[i];
        res += d * d;
    }
    REQUIRE(std::sqrt(res) < 1.0e-7);
    REQUIRE(iters > 0);
}
