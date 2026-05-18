// =============================================================================
// SimAll Beta - Dynamics Subsystem
// File   : src/dynamics/SixDofCoupler.hpp
// Week   : 18
//
// Six-degree-of-freedom rigid-body coupler that drives a SixDof body from
// fluid-side surface forces & moments.  Acts as the bridge between the
// solver subsystem (which reports per-time-step integrated wall forces
// and moments) and the SixDof integrator (which advances the body state).
//
// Coupling modes:
//
//   * Explicit         — F_n+1 ← solver, then body marches with F_n+1.
//   * Loose-implicit   — sub-iterate (force, move-mesh, F again) up to
//                        `nMaxSub` times or relative-position tolerance.
//   * Two-way w/ relax — under-relaxed loose-implicit; safe for bodies
//                        whose added-mass is comparable to inertial mass.
//
// External fluid forces and moments are expressed in the global (world)
// frame; the SixDof integrator handles the world ↔ body transformation.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace simall::dynamics {

struct Vec3 { double x = 0.0, y = 0.0, z = 0.0; };

enum class SixDofCouplingMode : std::uint8_t { Explicit, LooseImplicit, UnderRelaxed };

struct RigidBodyState {
    double  mass        = 1.0;
    Vec3    inertiaDiag { 1.0, 1.0, 1.0 };
    Vec3    position;
    Vec3    velocity;
    Vec3    angularVelocity;             // body frame, rad/s
    std::array<double, 9> rotation{ 1,0,0,  0,1,0,  0,0,1 };  // row-major R
};

struct CouplingOptions {
    SixDofCouplingMode mode      = SixDofCouplingMode::Explicit;
    std::size_t        nMaxSub   = 5;
    double             relax     = 0.7;     // 0..1
    double             tolPos    = 1e-6;
};

class SixDofCoupler {
public:
    explicit SixDofCoupler(CouplingOptions opt = {}) : opt_(opt) {}

    /// Advance the body by `dt` under fluid force `F` and moment `M`.
    /// `solverReevaluate` is invoked once per sub-iteration in loose-implicit
    /// and under-relaxed modes; it returns the *new* fluid (F,M) sampled
    /// at the *predicted* body pose.  Skipped entirely in Explicit mode.
    void advance(RigidBodyState&       state,
                  Vec3                   F,
                  Vec3                   M,
                  double                 dt,
                  void (*solverReevaluate)(const RigidBodyState&, Vec3&, Vec3&) = nullptr);

    [[nodiscard]] std::size_t last_subiterations() const noexcept { return lastSub_; }

private:
    CouplingOptions opt_;
    std::size_t     lastSub_ = 0;

    static RigidBodyState explicit_step(const RigidBodyState& s, const Vec3& F, const Vec3& M, double dt);
};

}  // namespace simall::dynamics
