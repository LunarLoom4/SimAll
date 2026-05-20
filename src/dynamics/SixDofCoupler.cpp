// =============================================================================
// SimAll Beta - Dynamics Subsystem
// File   : src/dynamics/SixDofCoupler.cpp
// =============================================================================
#include "dynamics/SixDofCoupler.hpp"

#include <algorithm>
#include <cmath>

namespace simall::dynamics
{

namespace
{

Vec3 mat3_vec(const std::array<double, 9>& R, const Vec3& v)
{
    return {R[0] * v.x + R[1] * v.y + R[2] * v.z,
            R[3] * v.x + R[4] * v.y + R[5] * v.z,
            R[6] * v.x + R[7] * v.y + R[8] * v.z};
}
Vec3 mat3T_vec(const std::array<double, 9>& R, const Vec3& v)
{
    return {R[0] * v.x + R[3] * v.y + R[6] * v.z,
            R[1] * v.x + R[4] * v.y + R[7] * v.z,
            R[2] * v.x + R[5] * v.y + R[8] * v.z};
}

// Rodrigues-formula incremental rotation about a (small) world-frame
// angular-velocity vector × dt.  Composes onto the existing rotation.
std::array<double, 9> integrate_rotation(const std::array<double, 9>& R,
                                         const Vec3& omegaWorld,
                                         double dt)
{
    const double theta = std::sqrt(omegaWorld.x * omegaWorld.x + omegaWorld.y * omegaWorld.y
                                   + omegaWorld.z * omegaWorld.z)
                         * dt;
    if (theta < 1e-14)
        return R;
    const Vec3 axis{
        omegaWorld.x / (theta / dt), omegaWorld.y / (theta / dt), omegaWorld.z / (theta / dt)};
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    const double C = 1.0 - c;
    const std::array<double, 9> dR = {c + axis.x * axis.x * C,
                                      axis.x * axis.y * C - axis.z * s,
                                      axis.x * axis.z * C + axis.y * s,
                                      axis.y * axis.x * C + axis.z * s,
                                      c + axis.y * axis.y * C,
                                      axis.y * axis.z * C - axis.x * s,
                                      axis.z * axis.x * C - axis.y * s,
                                      axis.z * axis.y * C + axis.x * s,
                                      c + axis.z * axis.z * C};
    std::array<double, 9> out{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s2 = 0.0;
            for (int k = 0; k < 3; ++k)
                s2 += dR[i * 3 + k] * R[k * 3 + j];
            out[i * 3 + j] = s2;
        }
    return out;
}

} // namespace

RigidBodyState SixDofCoupler::explicit_step(const RigidBodyState& s0,
                                            const Vec3& F,
                                            const Vec3& M,
                                            double dt)
{
    RigidBodyState s = s0;
    // Translational: a = F/m, semi-implicit Euler.
    s.velocity.x += (F.x / s.mass) * dt;
    s.velocity.y += (F.y / s.mass) * dt;
    s.velocity.z += (F.z / s.mass) * dt;
    s.position.x += s.velocity.x * dt;
    s.position.y += s.velocity.y * dt;
    s.position.z += s.velocity.z * dt;
    // Rotational: I_body · ω̇ = M_body - ω × (I_body · ω).  Express M in body frame.
    const Vec3 Mbody = mat3T_vec(s0.rotation, M);
    const Vec3 Iw{s.inertiaDiag.x * s.angularVelocity.x,
                  s.inertiaDiag.y * s.angularVelocity.y,
                  s.inertiaDiag.z * s.angularVelocity.z};
    const Vec3 wxIw{s.angularVelocity.y * Iw.z - s.angularVelocity.z * Iw.y,
                    s.angularVelocity.z * Iw.x - s.angularVelocity.x * Iw.z,
                    s.angularVelocity.x * Iw.y - s.angularVelocity.y * Iw.x};
    s.angularVelocity.x += dt * (Mbody.x - wxIw.x) / std::max(s.inertiaDiag.x, 1e-30);
    s.angularVelocity.y += dt * (Mbody.y - wxIw.y) / std::max(s.inertiaDiag.y, 1e-30);
    s.angularVelocity.z += dt * (Mbody.z - wxIw.z) / std::max(s.inertiaDiag.z, 1e-30);
    const Vec3 omegaWorld = mat3_vec(s.rotation, s.angularVelocity);
    s.rotation = integrate_rotation(s.rotation, omegaWorld, dt);
    return s;
}

void SixDofCoupler::advance(RigidBodyState& state,
                            Vec3 F,
                            Vec3 M,
                            double dt,
                            void (*solverReevaluate)(const RigidBodyState&, Vec3&, Vec3&))
{
    lastSub_ = 0;
    if (opt_.mode == SixDofCouplingMode::Explicit || !solverReevaluate) {
        state = explicit_step(state, F, M, dt);
        return;
    }
    RigidBodyState prev = state;
    for (std::size_t k = 0; k < opt_.nMaxSub; ++k) {
        RigidBodyState candidate = explicit_step(state, F, M, dt);
        const double dp = std::sqrt(std::pow(candidate.position.x - prev.position.x, 2.0)
                                    + std::pow(candidate.position.y - prev.position.y, 2.0)
                                    + std::pow(candidate.position.z - prev.position.z, 2.0));
        prev = candidate;
        // Apply optional under-relaxation in mixed mode.
        if (opt_.mode == SixDofCouplingMode::UnderRelaxed) {
            candidate.velocity.x =
                (1.0 - opt_.relax) * state.velocity.x + opt_.relax * candidate.velocity.x;
            candidate.velocity.y =
                (1.0 - opt_.relax) * state.velocity.y + opt_.relax * candidate.velocity.y;
            candidate.velocity.z =
                (1.0 - opt_.relax) * state.velocity.z + opt_.relax * candidate.velocity.z;
        }
        ++lastSub_;
        Vec3 Fnew = F, Mnew = M;
        solverReevaluate(candidate, Fnew, Mnew);
        F = Fnew;
        M = Mnew;
        if (dp < opt_.tolPos) {
            state = candidate;
            return;
        }
    }
    state = prev;
}

} // namespace simall::dynamics
