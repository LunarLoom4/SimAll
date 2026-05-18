// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/AndersonAcceleration.cpp
// =============================================================================
#include "solver/AndersonAcceleration.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace simall::solver {

namespace {
// Solve the m×m symmetric positive-definite system M γ = b in-place
// using Cholesky factorisation. Returns false if M is not SPD (rank
// deficient) — caller should drop oldest column and retry.
bool spd_solve(std::vector<std::vector<double>>& M, std::vector<double>& b) {
    const std::size_t n = M.size();
    if (n == 0) return true;
    // Cholesky: M = L L^T.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double s = M[i][j];
            for (std::size_t k = 0; k < j; ++k) s -= M[i][k] * M[j][k];
            if (i == j) {
                if (s <= 1e-18) return false;
                M[i][i] = std::sqrt(s);
            } else {
                M[i][j] = s / M[j][j];
            }
        }
    }
    // Forward L y = b.
    for (std::size_t i = 0; i < n; ++i) {
        double s = b[i];
        for (std::size_t k = 0; k < i; ++k) s -= M[i][k] * b[k];
        b[i] = s / M[i][i];
    }
    // Backward L^T x = y.
    for (std::size_t i = n; i-- > 0; ) {
        double s = b[i];
        for (std::size_t k = i+1; k < n; ++k) s -= M[k][i] * b[k];
        b[i] = s / M[i][i];
    }
    return true;
}
}  // namespace

AndersonAcceleration::AndersonAcceleration(std::size_t depth, double beta)
    : depth_(std::max<std::size_t>(1, depth)), beta_(beta) {}

void AndersonAcceleration::update(std::vector<double>& x,
                                  const std::vector<double>& g)
{
    const std::size_t n = x.size();
    if (g.size() != n)
        throw std::runtime_error("Anderson: x and g size mismatch");

    std::vector<double> f(n);
    for (std::size_t i = 0; i < n; ++i) f[i] = g[i] - x[i];

    if (have_prev_) {
        std::vector<double> dx(n), df(n);
        for (std::size_t i = 0; i < n; ++i) {
            dx[i] = x[i] - x_prev_[i];
            df[i] = f[i] - f_prev_[i];
        }
        dx_hist_.push_back(std::move(dx));
        df_hist_.push_back(std::move(df));
        while (dx_hist_.size() > depth_) { dx_hist_.pop_front(); df_hist_.pop_front(); }
    }
    x_prev_ = x; f_prev_ = f; have_prev_ = true;

    const std::size_t m = df_hist_.size();
    std::vector<double> gamma(m, 0.0);
    if (m > 0) {
        // Solve (F^T F) γ = F^T f via Cholesky on the m×m Gram matrix.
        std::vector<std::vector<double>> Gm(m, std::vector<double>(m, 0.0));
        std::vector<double> rhs(m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = i; j < m; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < n; ++k) s += df_hist_[i][k] * df_hist_[j][k];
                Gm[i][j] = Gm[j][i] = s;
            }
            double s = 0.0;
            for (std::size_t k = 0; k < n; ++k) s += df_hist_[i][k] * f[k];
            rhs[i] = s;
        }
        // Tikhonov regularisation to handle near-rank-deficient history.
        double trace = 0.0; for (std::size_t i = 0; i < m; ++i) trace += Gm[i][i];
        const double reg = 1e-12 * std::max(1.0, trace / static_cast<double>(m));
        for (std::size_t i = 0; i < m; ++i) Gm[i][i] += reg;
        if (!spd_solve(Gm, rhs)) {
            // Drop oldest and clear γ — fall back to plain Picard step.
            dx_hist_.pop_front(); df_hist_.pop_front();
            std::fill(gamma.begin(), gamma.end(), 0.0);
        } else {
            gamma = rhs;
        }
    }

    // x_{k+1} = x_k + β f_k - (X + β F) γ.
    for (std::size_t i = 0; i < n; ++i) {
        double acc = x[i] + beta_ * f[i];
        for (std::size_t j = 0; j < gamma.size(); ++j)
            acc -= (dx_hist_[j][i] + beta_ * df_hist_[j][i]) * gamma[j];
        x[i] = acc;
    }
}

}  // namespace simall::solver
