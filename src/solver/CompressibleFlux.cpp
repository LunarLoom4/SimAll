// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/CompressibleFlux.cpp
// =============================================================================
#include "solver/CompressibleFlux.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver
{

namespace
{

/// Roe-averaged sound speed estimates used to bound wave speeds (Einfeldt).
inline void wave_speeds(double rhoL,
                        double uL,
                        double pL,
                        double aL,
                        double rhoR,
                        double uR,
                        double pR,
                        double aR,
                        double& SL,
                        double& SR)
{
    // Roe averages
    const double sqrL = std::sqrt(rhoL), sqrR = std::sqrt(rhoR);
    const double uTilde = (sqrL * uL + sqrR * uR) / (sqrL + sqrR);
    const double HL = (pL / (rhoL)) + 0.5 * uL * uL; // partial enthalpy term
    const double HR = (pR / (rhoR)) + 0.5 * uR * uR;
    (void) HL;
    (void) HR;
    const double aTilde = 0.5 * (aL + aR);
    SL = std::min(uL - aL, uTilde - aTilde);
    SR = std::max(uR + aR, uTilde + aTilde);
}

} // namespace

void CompressibleFlux::rotate_to_face(const util::Vec3d& n,
                                      double& u,
                                      double& v,
                                      double& w,
                                      double& tx_x,
                                      double& tx_y,
                                      double& tx_z,
                                      double& ty_x,
                                      double& ty_y,
                                      double& ty_z)
{
    // Build orthonormal basis (n, tx, ty)
    util::Vec3d ref{1, 0, 0};
    if (std::abs(n.x) > 0.9)
        ref = {0, 1, 0};
    util::Vec3d tx{n.y * ref.z - n.z * ref.y, n.z * ref.x - n.x * ref.z, n.x * ref.y - n.y * ref.x};
    const double txn = std::sqrt(tx.x * tx.x + tx.y * tx.y + tx.z * tx.z);
    tx.x /= txn;
    tx.y /= txn;
    tx.z /= txn;
    util::Vec3d ty{n.y * tx.z - n.z * tx.y, n.z * tx.x - n.x * tx.z, n.x * tx.y - n.y * tx.x};
    tx_x = tx.x;
    tx_y = tx.y;
    tx_z = tx.z;
    ty_x = ty.x;
    ty_y = ty.y;
    ty_z = ty.z;
    const double un = u * n.x + v * n.y + w * n.z;
    const double ut = u * tx.x + v * tx.y + w * tx.z;
    const double us = u * ty.x + v * ty.y + w * ty.z;
    u = un;
    v = ut;
    w = us;
}

ConsFlux CompressibleFlux::hllc(const PrimState& Ls,
                                const PrimState& Rs,
                                const util::Vec3d& nIn,
                                double A) const
{
    // Normalise normal (caller may pass area vector already split).
    const double nn = std::sqrt(nIn.x * nIn.x + nIn.y * nIn.y + nIn.z * nIn.z);
    util::Vec3d n =
        (nn > 1e-30) ? util::Vec3d{nIn.x / nn, nIn.y / nn, nIn.z / nn} : util::Vec3d{1, 0, 0};

    PrimState L = Ls, R = Rs;
    double txx, txy, txz, tyx, tyy, tyz;
    rotate_to_face(n, L.u, L.v, L.w, txx, txy, txz, tyx, tyy, tyz);
    {
        double dummy_txx, dummy_txy, dummy_txz, dummy_tyx, dummy_tyy, dummy_tyz;
        rotate_to_face(
            n, R.u, R.v, R.w, dummy_txx, dummy_txy, dummy_txz, dummy_tyx, dummy_tyy, dummy_tyz);
    }

    const double g = gamma_;
    const double aL = std::sqrt(g * L.p / std::max(L.rho, 1e-30));
    const double aR = std::sqrt(g * R.p / std::max(R.rho, 1e-30));
    const double EL = L.p / ((g - 1.0) * L.rho) + 0.5 * (L.u * L.u + L.v * L.v + L.w * L.w);
    const double ER = R.p / ((g - 1.0) * R.rho) + 0.5 * (R.u * R.u + R.v * R.v + R.w * R.w);

    double SL, SR;
    wave_speeds(L.rho, L.u, L.p, aL, R.rho, R.u, R.p, aR, SL, SR);
    // Contact wave (Toro eq. 10.37)
    const double Sstar = (R.p - L.p + L.rho * L.u * (SL - L.u) - R.rho * R.u * (SR - R.u))
                         / (L.rho * (SL - L.u) - R.rho * (SR - R.u) + 1e-30);

    auto flux_state = [&](const PrimState& S, double E) {
        return ConsFlux{S.rho * S.u,
                        S.rho * S.u * S.u + S.p,
                        S.rho * S.u * S.v,
                        S.rho * S.u * S.w,
                        S.u * (S.rho * E + S.p)};
    };
    auto cons_state = [&](const PrimState& S, double E) {
        return ConsFlux{S.rho, S.rho * S.u, S.rho * S.v, S.rho * S.w, S.rho * E};
    };

    ConsFlux F;
    if (SL >= 0.0) {
        F = flux_state(L, EL);
    } else if (SR <= 0.0) {
        F = flux_state(R, ER);
    } else {
        // Star region (HLLC)
        auto star = [&](const PrimState& K, double EK, double SK) {
            const double pre = K.rho * (SK - K.u) / (SK - Sstar + 1e-30);
            const double Estar = EK + (Sstar - K.u) * (Sstar + K.p / (K.rho * (SK - K.u + 1e-30)));
            return ConsFlux{pre, pre * Sstar, pre * K.v, pre * K.w, pre * Estar};
        };
        if (Sstar >= 0.0) {
            const ConsFlux FL = flux_state(L, EL);
            const ConsFlux UL = cons_state(L, EL);
            const ConsFlux US = star(L, EL, SL);
            for (int i = 0; i < 5; ++i)
                F[i] = FL[i] + SL * (US[i] - UL[i]);
        } else {
            const ConsFlux FR = flux_state(R, ER);
            const ConsFlux UR = cons_state(R, ER);
            const ConsFlux US = star(R, ER, SR);
            for (int i = 0; i < 5; ++i)
                F[i] = FR[i] + SR * (US[i] - UR[i]);
        }
    }

    // Rotate momentum components back to global frame, scale by A.
    const double Fx = F[1] * n.x + F[2] * txx + F[3] * tyx;
    const double Fy = F[1] * n.y + F[2] * txy + F[3] * tyy;
    const double Fz = F[1] * n.z + F[2] * txz + F[3] * tyz;
    F[1] = Fx;
    F[2] = Fy;
    F[3] = Fz;
    for (int i = 0; i < 5; ++i)
        F[i] *= A;
    return F;
}

} // namespace simall::solver
