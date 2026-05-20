// =============================================================================
// SimAll Beta - Dynamics Subsystem
// File   : src/dynamics/SixDof.hpp
// Phase  : 14.5 — Six Degree-of-Freedom rigid-body integrator. Couples
// to the fluid solver through io::ForcesReport: each fluid step computes
// the resultant force F and moment M on the body's wetted surface, and
// SixDof.advance(F, M, dt) updates linear & angular state using a
// semi-implicit Euler (symplectic Newmark-β for orientation via
// quaternion exponential map).
//
//   m  dv/dt = F + m g                  (linear in world frame)
//   I  dω/dt = M - ω × (I ω)            (Euler eqns in body frame)
//   dx/dt    = v
//   dq/dt    = ½  Ω(ω_b)  q
//
// where Ω(ω) is the 4×4 quaternion rate matrix. Used together with the
// RBF morpher (mesh follows the body) and the IBM / sliding interface
// systems (no mesh deformation, body acts as a moving immersed surface).
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"

#include <array>

namespace simall::dynamics
{

struct Quat
{
    double w = 1, x = 0, y = 0, z = 0;
};

struct SixDofState
{
    util::Vec3d position{0, 0, 0};
    util::Vec3d velocity{0, 0, 0};
    Quat orientation{};                  // unit quaternion (w, x, y, z)
    util::Vec3d angularVelBody{0, 0, 0}; // ω expressed in body frame
};

struct SixDofParams
{
    double mass = 1.0;
    /// Body-frame inertia tensor (assumed principal-axis aligned, so only
    /// diagonal entries are stored). If you have a full 3×3 tensor, rotate
    /// the body axes onto its eigenvectors first.
    util::Vec3d inertiaDiag{1, 1, 1};
    util::Vec3d gravity{0, 0, -9.81};
    /// Per-DOF lock mask (1 = free, 0 = constrained to zero motion).
    std::array<int, 6> dofFree{1, 1, 1, 1, 1, 1};
};

class SixDof
{
public:
    void initialize(SixDofParams p, SixDofState s0);

    /// Advance one fluid-coupling step under external world-frame force F
    /// and world-frame moment M (both about the body centre of mass).
    void advance(const util::Vec3d& Fworld, const util::Vec3d& Mworld, double dt);

    const SixDofState& state() const noexcept { return state_; }
    const SixDofParams& params() const noexcept { return params_; }

    /// Rotate a body-frame vector into the world frame using current
    /// orientation (used by RbfMorpher to transform surface displacements).
    util::Vec3d body_to_world(const util::Vec3d& v) const;
    util::Vec3d world_to_body(const util::Vec3d& v) const;

private:
    SixDofParams params_{};
    SixDofState state_{};
};

} // namespace simall::dynamics
