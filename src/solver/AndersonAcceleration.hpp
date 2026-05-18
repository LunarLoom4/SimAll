// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/AndersonAcceleration.hpp
// Phase  : 7.7 — Anderson mixing (Type-II) for fixed-point convergence
// acceleration of the outer SIMPLE / segregated pressure-velocity loop.
//
// Given a fixed-point iteration x_{k+1} = g(x_k), Anderson(m) replaces
// it with a least-squares correction over the last m residuals:
//
//   x_{k+1} = x_k + β f_k - (X_k + β F_k) γ
//   γ = arg min ‖f_k + F_k γ‖₂        (unconstrained QR via normal eqs)
//
// where f_k = g(x_k) - x_k is the residual, X_k = [Δx_{k-m+1} … Δx_{k-1}],
// F_k = [Δf_{k-m+1} … Δf_{k-1}]. The mixing factor β ∈ (0, 1] (typ. 1).
// =============================================================================
#pragma once

#include <cstddef>
#include <deque>
#include <vector>

namespace simall::solver {

class AndersonAcceleration {
public:
    /// `depth` is the Anderson history length m (use 3–10 in practice).
    /// `beta`  is the relaxation factor.
    AndersonAcceleration(std::size_t depth = 5, double beta = 1.0);

    /// Push the most recent fixed-point pair (x_k, g(x_k)) and produce the
    /// next iterate x_{k+1}. Replaces `x` in place.
    void update(std::vector<double>& x, const std::vector<double>& g);

    /// Clear the history (call when restarting the solver).
    void reset() { dx_hist_.clear(); df_hist_.clear(); have_prev_ = false; }

    std::size_t depth()  const noexcept { return depth_; }
    std::size_t window() const noexcept { return dx_hist_.size(); }

private:
    std::size_t depth_;
    double      beta_;
    bool        have_prev_ = false;
    std::vector<double>            x_prev_;
    std::vector<double>            f_prev_;
    std::deque<std::vector<double>> dx_hist_;
    std::deque<std::vector<double>> df_hist_;
};

}  // namespace simall::solver
