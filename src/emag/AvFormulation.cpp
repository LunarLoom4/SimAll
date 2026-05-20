// =============================================================================
// SimAll Beta - Electromagnetics Subsystem
// File   : src/emag/AvFormulation.cpp
// =============================================================================
#include "emag/AvFormulation.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace simall::emag
{

namespace
{

double det3(const double m[9])
{
    return m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
           + m[2] * (m[3] * m[7] - m[4] * m[6]);
}

// Compute the 4 shape-function gradients ∇N_i (constant on a linear tet) and
// the tet volume.  ∇N_i = b_i / (6V) where b_i are co-factors of the [1 x y z]
// matrix.  We return ∇N as a 4×3 array and the signed volume.
void tet_grads(const std::array<std::array<double, 3>, 4>& X,
               std::array<std::array<double, 3>, 4>& grads,
               double& volume)
{
    const double M[16] = {1.0,
                          X[0][0],
                          X[0][1],
                          X[0][2],
                          1.0,
                          X[1][0],
                          X[1][1],
                          X[1][2],
                          1.0,
                          X[2][0],
                          X[2][1],
                          X[2][2],
                          1.0,
                          X[3][0],
                          X[3][1],
                          X[3][2]};
    // Volume = (1/6) det
    auto det4 = [&](const double m[16]) {
        // Expand along first column (all ones).
        const double a[9] = {m[5], m[6], m[7], m[9], m[10], m[11], m[13], m[14], m[15]};
        const double b[9] = {m[1], m[2], m[3], m[9], m[10], m[11], m[13], m[14], m[15]};
        const double c[9] = {m[1], m[2], m[3], m[5], m[6], m[7], m[13], m[14], m[15]};
        const double d[9] = {m[1], m[2], m[3], m[5], m[6], m[7], m[9], m[10], m[11]};
        return det3(a) - det3(b) + det3(c) - det3(d);
    };
    const double D = det4(M);
    volume = D / 6.0;
    // ∂N_i/∂x_j = co-factor C_{i+1, j+1} of M scaled by 1/D.
    auto cofactor = [&](int row, int col) {
        double sub[9];
        int k = 0;
        for (int r = 0; r < 4; ++r) {
            if (r == row)
                continue;
            for (int c = 0; c < 4; ++c) {
                if (c == col)
                    continue;
                sub[k++] = M[r * 4 + c];
            }
        }
        const double sign = ((row + col) & 1) ? -1.0 : 1.0;
        return sign * det3(sub);
    };
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            grads[i][j] = cofactor(i, j + 1) / D;
        }
    }
}

} // namespace

ElementMatrices assemble_tet4(const Tet4& e)
{
    ElementMatrices em;
    std::array<std::array<double, 3>, 4> grads;
    tet_grads(e.coords, grads, em.volume);
    const double V = std::abs(em.volume);
    const double invMu = 1.0 / std::max(e.mu, 1e-30);

    // K_A_ij = (1/μ) ∇N_i · ∇N_j · V
    // K_V_ij = σ      ∇N_i · ∇N_j · V
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const double gg =
                grads[i][0] * grads[j][0] + grads[i][1] * grads[j][1] + grads[i][2] * grads[j][2];
            em.K_A[i][j] = invMu * gg * V;
            em.K_V[i][j] = e.sigma * gg * V;
        }
    }
    // M_A: σ ∫ N_i N_j dV; on linear tet, ∫ N_i N_j dV = V/20 (i==j) or V/20·(off-diag = V/20 too?
    // actually V/20 for i==j and V/20·... ) Correct closed form: ∫ N_i N_j = V/20 · (1+δ_ij), so
    // diag = V/10, off-diag = V/20.
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            em.M_A[i][j] = e.sigma * V * ((i == j) ? 1.0 / 10.0 : 1.0 / 20.0);
    return em;
}

CsrFromTets build_av_scalar_system(const std::vector<Tet4>& mesh, std::size_t n)
{
    CsrFromTets sys;
    sys.n = n;
    // First pass: row-wise std::map of col→value.
    std::vector<std::map<std::size_t, double>> rows(n);
    for (const auto& e : mesh) {
        auto em = assemble_tet4(e);
        for (int i = 0; i < 4; ++i) {
            const std::size_t r = e.nodes[i];
            if (r >= n)
                continue;
            for (int j = 0; j < 4; ++j) {
                const std::size_t c = e.nodes[j];
                if (c >= n)
                    continue;
                rows[r][c] += em.K_V[i][j];
            }
        }
    }
    sys.rowPtr.assign(n + 1, 0);
    for (std::size_t i = 0; i < n; ++i)
        sys.rowPtr[i + 1] = sys.rowPtr[i] + rows[i].size();
    sys.colIdx.resize(sys.rowPtr.back());
    sys.values.resize(sys.rowPtr.back());
    std::size_t k = 0;
    for (std::size_t i = 0; i < n; ++i) {
        for (const auto& [c, v] : rows[i]) {
            sys.colIdx[k] = c;
            sys.values[k] = v;
            ++k;
        }
    }
    return sys;
}

} // namespace simall::emag
