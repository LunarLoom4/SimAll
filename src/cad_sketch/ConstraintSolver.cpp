// =============================================================================
// SimAll Beta -- Parametric CAD / Sketcher
// File   : src/cad_sketch/ConstraintSolver.cpp
// Phase  : 23 Pass 23.1
//
// Damped Gauss-Newton (Levenberg-Marquardt) solver with rank diagnostics.
//
// Residuals are evaluated by `emit_constraint_residuals` which dispatches on
// constraint kind.  The Jacobian is built by central finite differences --
// at sketcher scale (n, m typically << 1000) this is fast, exact-up-to-
// rounding, and avoids the maintenance burden of analytic per-constraint
// partials.  A later pass can add analytic Jacobians for the hot cases.
//
// The damped normal-equations system  (J^T J + lambda I) dx = -J^T r  is
// solved by Eigen's column-pivoting Householder QR, which also yields a
// numerical rank estimate (with our `rank_tol`).  Pivot columns rejected
// by QR are reported as `degenerate_parameters` so the UI can highlight
// the under-/over-constrained subset.
// =============================================================================
#include "cad_sketch/ConstraintSolver.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <unordered_map>

namespace simall::cad::sketch {

namespace {

// ---------------------------------------------------------------------------
// Unknowns map: column j of x corresponds to non-fixed ParameterId col_to_pid[j]
// ---------------------------------------------------------------------------
struct UnknownsMap {
    std::vector<ParameterId>             col_to_pid;
    std::unordered_map<ParameterId, int> pid_to_col;
};

UnknownsMap build_unknowns(const Sketch& s) {
    UnknownsMap u;
    for (const auto& p : s.parameters()) {
        if (!p.fixed) {
            u.pid_to_col.emplace(p.id, static_cast<int>(u.col_to_pid.size()));
            u.col_to_pid.push_back(p.id);
        }
    }
    return u;
}

// ---------------------------------------------------------------------------
// Helpers for fetching entity geometry by reading current parameter values.
// All return current values from the sketch (not cached) so the solver can
// re-evaluate after each parameter perturbation.
// ---------------------------------------------------------------------------
inline Vec2 pt(const Sketch& s, EntityId e) {
    const auto& ent = s.entity(e);
    return { s.get_parameter(ent.params[0]), s.get_parameter(ent.params[1]) };
}
inline Vec2 lineA(const Sketch& s, EntityId e) {
    const auto& ent = s.entity(e);
    return { s.get_parameter(ent.params[0]), s.get_parameter(ent.params[1]) };
}
inline Vec2 lineB(const Sketch& s, EntityId e) {
    const auto& ent = s.entity(e);
    return { s.get_parameter(ent.params[2]), s.get_parameter(ent.params[3]) };
}
inline Vec2 circCenter(const Sketch& s, EntityId e) {
    const auto& ent = s.entity(e);
    return { s.get_parameter(ent.params[0]), s.get_parameter(ent.params[1]) };
}
inline double circRadius(const Sketch& s, EntityId e) {
    const auto& ent = s.entity(e);
    return s.get_parameter(ent.params[2]);
}

// ---------------------------------------------------------------------------
// Per-constraint residual emission.  Each kind appends 1 or 2 doubles.
// Constraints that don't contribute residuals (e.g. Fix, which is encoded
// by marking parameters as fixed at construction time) emit nothing.
// ---------------------------------------------------------------------------
void emit_constraint_residuals(const Sketch& s, const SketchConstraint& c,
                               std::vector<double>& out) {
    switch (c.kind) {
        case ConstraintKind::Coincident: {
            const Vec2 p1 = pt(s, c.a), p2 = pt(s, c.b);
            out.push_back(p1.x - p2.x);
            out.push_back(p1.y - p2.y);
            return;
        }
        case ConstraintKind::Horizontal: {
            const Vec2 a = lineA(s, c.a), b = lineB(s, c.a);
            out.push_back(a.y - b.y);
            return;
        }
        case ConstraintKind::Vertical: {
            const Vec2 a = lineA(s, c.a), b = lineB(s, c.a);
            out.push_back(a.x - b.x);
            return;
        }
        case ConstraintKind::Parallel: {
            const Vec2 d1 = lineB(s, c.a) - lineA(s, c.a);
            const Vec2 d2 = lineB(s, c.b) - lineA(s, c.b);
            out.push_back(d1.cross(d2));            // == 0  <=>  parallel
            return;
        }
        case ConstraintKind::Perpendicular: {
            const Vec2 d1 = lineB(s, c.a) - lineA(s, c.a);
            const Vec2 d2 = lineB(s, c.b) - lineA(s, c.b);
            out.push_back(d1.dot(d2));              // == 0  <=>  perpendicular
            return;
        }
        case ConstraintKind::Distance: {
            const Vec2 p1 = pt(s, c.a), p2 = pt(s, c.b);
            out.push_back((p1 - p2).norm() - c.value.value_or(0.0));
            return;
        }
        case ConstraintKind::Angle: {
            // dot(d1,d2) - cos(theta) * |d1| * |d2|  = 0
            const Vec2 d1 = lineB(s, c.a) - lineA(s, c.a);
            const Vec2 d2 = lineB(s, c.b) - lineA(s, c.b);
            const double n1 = d1.norm(), n2 = d2.norm();
            const double target = std::cos(c.value.value_or(0.0)) * n1 * n2;
            out.push_back(d1.dot(d2) - target);
            return;
        }
        case ConstraintKind::Radius: {
            out.push_back(circRadius(s, c.a) - c.value.value_or(0.0));
            return;
        }
        case ConstraintKind::PointOnLine: {
            const Vec2 p = pt(s, c.a);
            const Vec2 a = lineA(s, c.b), b = lineB(s, c.b);
            out.push_back((b - a).cross(p - a));    // signed area; 0 iff collinear
            return;
        }
        case ConstraintKind::PointOnCircle: {
            const Vec2 p = pt(s, c.a);
            const Vec2 cc = circCenter(s, c.b);
            out.push_back((p - cc).norm() - circRadius(s, c.b));
            return;
        }
        case ConstraintKind::EqualLength: {
            const Vec2 d1 = lineB(s, c.a) - lineA(s, c.a);
            const Vec2 d2 = lineB(s, c.b) - lineA(s, c.b);
            out.push_back(d1.norm2() - d2.norm2());
            return;
        }
        case ConstraintKind::Fix:
            // Encoded by pinning the underlying parameters (fixed=true) when
            // the constraint is created; contributes nothing to residual.
            return;
    }
}

std::vector<double> compute_residuals_impl(const Sketch& s) {
    std::vector<double> r;
    r.reserve(s.constraints().size() * 2);
    for (const auto& c : s.constraints()) emit_constraint_residuals(s, c, r);
    return r;
}

// ---------------------------------------------------------------------------
// Finite-difference Jacobian: J[i, j] = d r_i / d x_j, evaluated by central
// differences with step `eps`.  This perturbs Sketch parameters in place;
// callers must restore the baseline value (we save+restore around each col).
// ---------------------------------------------------------------------------
Eigen::MatrixXd jacobian_fd(Sketch& s, const UnknownsMap& u,
                            const std::vector<double>& /*r0_unused*/, double eps) {
    const int m = static_cast<int>(compute_residuals_impl(s).size());
    const int n = static_cast<int>(u.col_to_pid.size());
    Eigen::MatrixXd J(m, n);

    for (int j = 0; j < n; ++j) {
        const ParameterId pid = u.col_to_pid[j];
        const double      x0  = s.get_parameter(pid);
        // Scale eps by parameter magnitude for robustness.
        const double h = eps * (1.0 + std::abs(x0));

        s.set_parameter(pid, x0 + h);
        const std::vector<double> rp = compute_residuals_impl(s);

        s.set_parameter(pid, x0 - h);
        const std::vector<double> rm = compute_residuals_impl(s);

        s.set_parameter(pid, x0);   // restore

        for (int i = 0; i < m; ++i) {
            J(i, j) = (rp[i] - rm[i]) / (2.0 * h);
        }
    }
    return J;
}

inline double l2(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x * x;
    return std::sqrt(s);
}

inline bool any_nonfinite(const std::vector<double>& v) {
    for (double x : v) if (!std::isfinite(x)) return true;
    return false;
}

}  // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

std::vector<double> evaluate_residuals(const Sketch& s) {
    return compute_residuals_impl(s);
}

SolveReport solve(Sketch& s, const SolverConfig& cfg) {
    SolveReport rep;

    const UnknownsMap unk = build_unknowns(s);
    const int n = static_cast<int>(unk.col_to_pid.size());

    std::vector<double> r = compute_residuals_impl(s);
    const int m = static_cast<int>(r.size());

    rep.unknowns       = static_cast<std::size_t>(n);
    rep.residual_count = static_cast<std::size_t>(m);
    rep.residual_norm  = l2(r);
    rep.residual_history.push_back(rep.residual_norm);

    if (n == 0 || m == 0) {
        rep.status  = SolveStatus::Empty;
        rep.message = "No unknowns or no constraints.";
        return rep;
    }

    double      lambda = cfg.initial_damping;
    std::size_t iter   = 0;

    for (; iter < cfg.max_iterations; ++iter) {
        if (rep.residual_norm < cfg.residual_tol) break;

        // --- Jacobian ----------------------------------------------------------
        const Eigen::MatrixXd J = jacobian_fd(s, unk, r, cfg.fd_eps);
        Eigen::VectorXd       rv(m);
        for (int i = 0; i < m; ++i) rv[i] = r[i];

        // --- Damped normal equations  (J^T J + lambda I) dx = -J^T r ----------
        Eigen::MatrixXd H = J.transpose() * J;
        H.diagonal().array() += lambda;
        const Eigen::VectorXd g = -J.transpose() * rv;

        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(H);
        qr.setThreshold(cfg.rank_tol);
        const Eigen::VectorXd dx = qr.solve(g);

        if (!dx.allFinite()) {
            rep.status  = SolveStatus::Singular;
            rep.message = "Non-finite step computed (Jacobian breakdown).";
            break;
        }

        const double step_norm = dx.norm();

        // --- Trial step --------------------------------------------------------
        std::vector<double> x0(n);
        for (int j = 0; j < n; ++j) x0[j] = s.get_parameter(unk.col_to_pid[j]);
        for (int j = 0; j < n; ++j) s.set_parameter(unk.col_to_pid[j], x0[j] + dx[j]);

        std::vector<double> r_trial = compute_residuals_impl(s);
        const double        nrm_trial = any_nonfinite(r_trial)
                                          ? std::numeric_limits<double>::infinity()
                                          : l2(r_trial);

        if (nrm_trial < rep.residual_norm) {
            // Accept step -- shrink damping.
            r                  = std::move(r_trial);
            rep.residual_norm  = nrm_trial;
            rep.residual_history.push_back(rep.residual_norm);
            lambda             = std::max(lambda * cfg.damping_shrink, 1e-15);

            if (step_norm < cfg.step_tol) {
                rep.status = (rep.residual_norm < cfg.residual_tol)
                             ? SolveStatus::Converged
                             : SolveStatus::StalledStepTooSmall;
                ++iter;
                break;
            }
        } else {
            // Reject -- restore x, grow damping.
            for (int j = 0; j < n; ++j) s.set_parameter(unk.col_to_pid[j], x0[j]);
            lambda *= cfg.damping_grow;
            if (lambda > 1e20) {
                rep.status  = SolveStatus::Singular;
                rep.message = "Damping diverged; cannot reduce residual.";
                ++iter;
                break;
            }
        }
    }

    rep.iterations = iter;
    if (rep.status == SolveStatus::Empty) {
        // Status was not set during the loop -- decide based on residual norm.
        if (rep.residual_norm < cfg.residual_tol) rep.status = SolveStatus::Converged;
        else if (iter >= cfg.max_iterations)      rep.status = SolveStatus::StalledMaxIters;
    }

    // --- Final rank / degeneracy diagnostics --------------------------------
    {
        const Eigen::MatrixXd Jf = jacobian_fd(s, unk, r, cfg.fd_eps);
        if (Jf.size() > 0) {
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrf(Jf);
            qrf.setThreshold(cfg.rank_tol);
            rep.rank = static_cast<std::size_t>(qrf.rank());
            rep.dof  = (rep.unknowns > rep.rank) ? rep.unknowns - rep.rank : 0;

            // Columns of J ranked by QR pivoting: the last (n - rank) pivots
            // correspond to rank-deficient parameters.
            const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& P =
                qrf.colsPermutation();
            const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& perm = P.indices();
            for (Eigen::Index k = static_cast<Eigen::Index>(rep.rank); k < perm.size(); ++k) {
                rep.degenerate_parameters.push_back(unk.col_to_pid[static_cast<std::size_t>(perm[k])]);
            }
        }
    }

    // Refine status using rank info.
    if (rep.status == SolveStatus::StalledMaxIters || rep.status == SolveStatus::StalledStepTooSmall) {
        if (rep.rank < rep.unknowns && rep.residual_norm < cfg.residual_tol * 1e3) {
            rep.status  = SolveStatus::Underconstrained;
            rep.message = "Sketch is under-constrained: multiple solutions exist.";
        } else if (rep.rank == rep.unknowns) {
            rep.status  = SolveStatus::Overconstrained;
            rep.message = "Sketch is over-constrained: residual cannot reach tolerance.";
        }
    }
    if (rep.status == SolveStatus::Converged && rep.rank < rep.unknowns) {
        rep.message = "Converged, but sketch retains " +
                      std::to_string(rep.dof) + " free DOF.";
    }

    return rep;
}

}  // namespace simall::cad::sketch
