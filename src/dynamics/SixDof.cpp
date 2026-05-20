// =============================================================================
// SimAll Beta - Dynamics Subsystem
// File   : src/dynamics/SixDof.cpp
// =============================================================================
#include "dynamics/SixDof.hpp"

#include <cmath>

namespace simall::dynamics
{

namespace
{

Quat quat_normalize(Quat q)
{
    const double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    if (n < 1e-30)
        return {1, 0, 0, 0};
    return {q.w / n, q.x / n, q.y / n, q.z / n};
}

Quat quat_mul(const Quat& a, const Quat& b)
{
    return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

/// Quaternion exponential map of the half-angle vector θ/2 (body frame).
Quat quat_exp_half(const util::Vec3d& omega, double dt)
{
    const double phi =
        0.5 * dt * std::sqrt(omega.x * omega.x + omega.y * omega.y + omega.z * omega.z);
    if (phi < 1e-12)
        return {1, 0.5 * dt * omega.x, 0.5 * dt * omega.y, 0.5 * dt * omega.z};
    const double s =
        std::sin(phi) / std::sqrt(omega.x * omega.x + omega.y * omega.y + omega.z * omega.z);
    return {std::cos(phi), omega.x * s, omega.y * s, omega.z * s};
}

util::Vec3d rotate_by_quat(const util::Vec3d& v, const Quat& q)
{
    // v' = q · (0, v) · q*
    const Quat qv{0, v.x, v.y, v.z};
    const Quat qc{q.w, -q.x, -q.y, -q.z};
    const Quat r = quat_mul(quat_mul(q, qv), qc);
    return {r.x, r.y, r.z};
}

} // namespace

void SixDof::initialize(SixDofParams p, SixDofState s0)
{
    params_ = p;
    state_ = s0;
    state_.orientation = quat_normalize(state_.orientation);
}

util::Vec3d SixDof::body_to_world(const util::Vec3d& v) const
{
    return rotate_by_quat(v, state_.orientation);
}
util::Vec3d SixDof::world_to_body(const util::Vec3d& v) const
{
    const Quat qc{
        state_.orientation.w, -state_.orientation.x, -state_.orientation.y, -state_.orientation.z};
    return rotate_by_quat(v, qc);
}

void SixDof::advance(const util::Vec3d& Fworld, const util::Vec3d& Mworld, double dt)
{
    const auto& F = params_;
    // --- Linear: semi-implicit Euler in world frame ---
    util::Vec3d a{(Fworld.x + F.mass * F.gravity.x) / F.mass,
                  (Fworld.y + F.mass * F.gravity.y) / F.mass,
                  (Fworld.z + F.mass * F.gravity.z) / F.mass};
    if (F.dofFree[0])
        state_.velocity.x += dt * a.x;
    else
        state_.velocity.x = 0;
    if (F.dofFree[1])
        state_.velocity.y += dt * a.y;
    else
        state_.velocity.y = 0;
    if (F.dofFree[2])
        state_.velocity.z += dt * a.z;
    else
        state_.velocity.z = 0;
    state_.position.x += dt * state_.velocity.x;
    state_.position.y += dt * state_.velocity.y;
    state_.position.z += dt * state_.velocity.z;

    // --- Angular: Euler's equations in body frame ---
    const util::Vec3d Mbody = world_to_body(Mworld);
    const util::Vec3d& w = state_.angularVelBody;
    const util::Vec3d& I = F.inertiaDiag;
    // dω_b/dt = I^{-1} (M_b - ω × (I ω))
    const util::Vec3d Iw{I.x * w.x, I.y * w.y, I.z * w.z};
    const util::Vec3d wxIw{
        w.y * Iw.z - w.z * Iw.y, w.z * Iw.x - w.x * Iw.z, w.x * Iw.y - w.y * Iw.x};
    util::Vec3d aw{(Mbody.x - wxIw.x) / std::max(I.x, 1e-30),
                   (Mbody.y - wxIw.y) / std::max(I.y, 1e-30),
                   (Mbody.z - wxIw.z) / std::max(I.z, 1e-30)};
    if (F.dofFree[3])
        state_.angularVelBody.x += dt * aw.x;
    else
        state_.angularVelBody.x = 0;
    if (F.dofFree[4])
        state_.angularVelBody.y += dt * aw.y;
    else
        state_.angularVelBody.y = 0;
    if (F.dofFree[5])
        state_.angularVelBody.z += dt * aw.z;
    else
        state_.angularVelBody.z = 0;

    // Quaternion update via exponential map (preserves unit norm).
    const Quat dq = quat_exp_half(state_.angularVelBody, dt);
    state_.orientation = quat_normalize(quat_mul(state_.orientation, dq));
}

} // namespace simall::dynamics
