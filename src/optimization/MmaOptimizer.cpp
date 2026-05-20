// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/MmaOptimizer.cpp
// =============================================================================
#include "optimization/MmaOptimizer.hpp"

#include <algorithm>
#include <cmath>

namespace simall::optimization
{

MmaResult run_mma(MmaProblem p, MmaOptions opt)
{
    MmaResult r;
    const std::size_t n = p.n;
    const std::size_t m = p.m;
    if (n == 0 || p.x0.size() != n || p.xLow.size() != n || p.xUp.size() != n || !p.evaluate) {
        r.error = "Invalid MMA problem";
        return r;
    }
    r.x = p.x0;
    std::vector<double> xPrev = r.x;
    std::vector<double> xPP = r.x;
    std::vector<double> low(n), up(n);
    std::vector<double> df(n), dg(m * n);
    std::vector<double> g(m);
    double f = 0.0;

    for (std::size_t k = 0; k < opt.maxIter; ++k) {
        p.evaluate(r.x, f, g, df, dg);
        r.history.push_back(f);

        // ---- Asymptote update ---------------------------------------------
        for (std::size_t i = 0; i < n; ++i) {
            const double range = std::max(p.xUp[i] - p.xLow[i], 1e-30);
            if (k < 2) {
                low[i] = r.x[i] - opt.asyInit * range;
                up[i] = r.x[i] + opt.asyInit * range;
            } else {
                const double s = (r.x[i] - xPrev[i]) * (xPrev[i] - xPP[i]);
                double gamma = 1.0;
                if (s > 0)
                    gamma = opt.asyIncr;
                else if (s < 0)
                    gamma = opt.asyDecr;
                low[i] = r.x[i] - gamma * (xPrev[i] - low[i]);
                up[i] = r.x[i] + gamma * (up[i] - xPrev[i]);
            }
            // Tighten by move limit.
            const double ml = opt.moveLimit * range;
            low[i] = std::max(low[i], r.x[i] - ml);
            up[i] = std::min(up[i], r.x[i] + ml);
            low[i] = std::max(low[i], p.xLow[i]);
            up[i] = std::min(up[i], p.xUp[i]);
        }

        // ---- MMA convex subproblem: minimise separable approximation -----
        // Build p_i / q_i for objective and constraints (Svanberg form).
        // For a single objective with grad df: p_i⁰ = (up-x)²·max(df,0),
        //                                       q_i⁰ = (x-low)²·max(-df,0).
        std::vector<double> p0(n), q0(n);
        std::vector<double> pj(m * n), qj(m * n);
        for (std::size_t i = 0; i < n; ++i) {
            const double up_x = std::max(up[i] - r.x[i], 1e-12);
            const double x_lo = std::max(r.x[i] - low[i], 1e-12);
            p0[i] = up_x * up_x * std::max(df[i], 0.0);
            q0[i] = x_lo * x_lo * std::max(-df[i], 0.0);
            for (std::size_t j = 0; j < m; ++j) {
                pj[j * n + i] = up_x * up_x * std::max(dg[j * n + i], 0.0);
                qj[j * n + i] = x_lo * x_lo * std::max(-dg[j * n + i], 0.0);
            }
        }

        // Use projected gradient on x with a fixed penalty for constraint
        // violations (sufficient when m is small and we restart at each k).
        const double penalty = 1e3;
        std::vector<double> xs = r.x;
        for (std::size_t it = 0; it < opt.subIters; ++it) {
            std::vector<double> grad(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                const double up_x = std::max(up[i] - xs[i], 1e-12);
                const double x_lo = std::max(xs[i] - low[i], 1e-12);
                grad[i] = p0[i] / (up_x * up_x) - q0[i] / (x_lo * x_lo);
            }
            // Recompute approx g via the MMA expansion:
            for (std::size_t j = 0; j < m; ++j) {
                double gApprox = g[j];
                for (std::size_t i = 0; i < n; ++i) {
                    const double up_x = std::max(up[i] - xs[i], 1e-12);
                    const double x_lo = std::max(xs[i] - low[i], 1e-12);
                    gApprox += pj[j * n + i] / up_x + qj[j * n + i] / x_lo
                               - pj[j * n + i] / std::max(up[i] - r.x[i], 1e-12)
                               - qj[j * n + i] / std::max(r.x[i] - low[i], 1e-12);
                }
                if (gApprox > 0.0) {
                    for (std::size_t i = 0; i < n; ++i) {
                        const double up_x = std::max(up[i] - xs[i], 1e-12);
                        const double x_lo = std::max(xs[i] - low[i], 1e-12);
                        grad[i] +=
                            penalty * gApprox
                            * (pj[j * n + i] / (up_x * up_x) - qj[j * n + i] / (x_lo * x_lo));
                    }
                }
            }
            const double step = 1e-3;
            for (std::size_t i = 0; i < n; ++i) {
                xs[i] -= step * grad[i];
                xs[i] = std::min(std::max(xs[i], low[i]), up[i]);
            }
        }

        // ---- Convergence test --------------------------------------------
        double dx = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            dx = std::max(dx, std::abs(xs[i] - r.x[i]));
        }
        xPP = xPrev;
        xPrev = r.x;
        r.x = xs;
        r.iterations = k + 1;
        if (dx < opt.tolX)
            break;
    }
    p.evaluate(r.x, r.f, r.g, df, dg);
    r.ok = true;
    return r;
}

} // namespace simall::optimization
