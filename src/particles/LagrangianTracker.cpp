// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/LagrangianTracker.cpp
// =============================================================================
#include "particles/LagrangianTracker.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace simall::particles {

namespace {

inline double dot3(const util::Vec3d& a, const util::Vec3d& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

/// Returns true if cell centroid-to-point line stays on the owner side of
/// every face. (Robust enough for face-walking; not used as the canonical
/// inside test — see locate() below.)
bool point_in_cell(const meshing::Mesh& m, meshing::CellId c, const util::Vec3d& p) {
    const auto& C = m.cells();
    const auto& F = m.faces();
    const util::Vec3d ctr{ C.centroidX[c], C.centroidY[c], C.centroidZ[c] };
    const int s = C.faceOffsets[c];
    const int e = C.faceOffsets[c+1];
    for (int k = s; k < e; ++k) {
        const meshing::FaceId fid = C.faceIndices[k];
        util::Vec3d n{ F.areaX[fid], F.areaY[fid], F.areaZ[fid] };
        if (F.owner[fid] != c) { n.x = -n.x; n.y = -n.y; n.z = -n.z; }
        const util::Vec3d fc{ F.centroidX[fid], F.centroidY[fid], F.centroidZ[fid] };
        // Outward-pointing normal n; point p inside if (p-fc)·n ≤ 0.
        const util::Vec3d dp{ p.x - fc.x, p.y - fc.y, p.z - fc.z };
        if (dot3(dp, n) > 1e-12 * (std::abs(dot3({ctr.x-fc.x,ctr.y-fc.y,ctr.z-fc.z}, n)) + 1e-30))
            return false;
    }
    return true;
}

}  // namespace

void LagrangianTracker::initialize(const meshing::Mesh& m,
                                   double rhoF, double muF,
                                   const ParticleSpec& s) {
    mesh_ = &m; rhoF_ = rhoF; muF_ = muF; spec_ = s;
    const std::size_t nC = m.cells().size();
    srcX_.assign(nC, 0.0);
    srcY_.assign(nC, 0.0);
    srcZ_.assign(nC, 0.0);
}

void LagrangianTracker::inject_point(const util::Vec3d& loc,
                                     const util::Vec3d& vel,
                                     std::size_t count) {
    if (!mesh_) return;
    // Find host cell by brute scan (one-time per injection).
    meshing::CellId host = static_cast<meshing::CellId>(-1);
    const std::size_t nC = mesh_->cells().size();
    for (std::size_t c = 0; c < nC; ++c) {
        if (point_in_cell(*mesh_, c, loc)) { host = c; break; }
    }
    if (host == static_cast<meshing::CellId>(-1)) {
        SIMALL_LOG_WARN("Particles", "inject_point: location outside domain");
        return;
    }
    const double V = (4.0/3.0) * M_PI * std::pow(0.5 * spec_.diameter, 3.0);
    const double m = spec_.density * V;
    const double I = 0.1 * m * spec_.diameter * spec_.diameter;  // (1/10) m d_p²
    parts_.reserve(parts_.size() + count);
    for (std::size_t i = 0; i < count; ++i) {
        ParticleState ps;
        ps.x = loc; ps.v = vel; ps.mass = m; ps.inertia = I;
        ps.cell = host; ps.active = true;
        parts_.push_back(ps);
    }
}

meshing::CellId LagrangianTracker::locate(const util::Vec3d& p,
                                          meshing::CellId hint) const {
    if (hint == static_cast<meshing::CellId>(-1)) return hint;
    if (point_in_cell(*mesh_, hint, p)) return hint;
    // Face-walk: move to neighbour through the face that point p violates
    // most strongly. Cap walk length to avoid infinite loops on bad meshes.
    meshing::CellId c = hint;
    const auto& C = mesh_->cells();
    const auto& F = mesh_->faces();
    for (int step = 0; step < 64; ++step) {
        const int s = C.faceOffsets[c];
        const int e = C.faceOffsets[c+1];
        double worst = 0.0;
        meshing::FaceId fexit = static_cast<meshing::FaceId>(-1);
        for (int k = s; k < e; ++k) {
            const meshing::FaceId fid = C.faceIndices[k];
            util::Vec3d n{ F.areaX[fid], F.areaY[fid], F.areaZ[fid] };
            if (F.owner[fid] != c) { n.x = -n.x; n.y = -n.y; n.z = -n.z; }
            const util::Vec3d fc{ F.centroidX[fid], F.centroidY[fid], F.centroidZ[fid] };
            const double d = dot3({p.x-fc.x, p.y-fc.y, p.z-fc.z}, n);
            if (d > worst) { worst = d; fexit = fid; }
        }
        if (fexit == static_cast<meshing::FaceId>(-1)) return c;
        const meshing::CellId nb = (F.owner[fexit] == c) ? F.neighbor[fexit]
                                                         : F.owner[fexit];
        if (nb == meshing::kBoundaryCell) return static_cast<meshing::CellId>(-1);
        c = nb;
        if (point_in_cell(*mesh_, c, p)) return c;
    }
    return c;
}

void LagrangianTracker::advance(double dt, const solver::FieldRegistry& fields) {
    if (!mesh_) return;
    auto* Uf = const_cast<solver::VectorField*>(
        const_cast<solver::FieldRegistry&>(fields).find_vector("U"));
    if (!Uf) return;

    std::fill(srcX_.begin(), srcX_.end(), 0.0);
    std::fill(srcY_.begin(), srcY_.end(), 0.0);
    std::fill(srcZ_.begin(), srcZ_.end(), 0.0);

    const auto& C = mesh_->cells();
    const std::size_t nC = C.size();

    // Vorticity ω_f = ∇×U per cell (used for Saffman & Magnus). Built only
    // when at least one lift mechanism is enabled to keep the hot path lean.
    std::vector<util::Vec3d> vortF;
    if (spec_.saffmanLift || spec_.magnusLift) {
        solver::LeastSquaresGradient G(*mesh_);
        solver::ScalarField ux(nC), uy(nC), uz(nC);
        std::copy(Uf->x.begin(), Uf->x.end(), ux.begin());
        std::copy(Uf->y.begin(), Uf->y.end(), uy.begin());
        std::copy(Uf->z.begin(), Uf->z.end(), uz.begin());
        solver::VectorField gx, gy, gz;
        G.evaluate(ux, gx); G.evaluate(uy, gy); G.evaluate(uz, gz);
        vortF.resize(nC);
        for (std::size_t c = 0; c < nC; ++c) {
            vortF[c] = { gz.y[c] - gy.z[c],
                         gx.z[c] - gz.x[c],
                         gy.x[c] - gx.y[c] };
        }
    }

    for (ParticleState& P : parts_) {
        if (!P.active) continue;
        if (P.cell == static_cast<meshing::CellId>(-1)) { P.active = false; continue; }

        // Cell-centred fluid velocity at particle position (zero-order).
        const util::Vec3d uF{ Uf->x[P.cell], Uf->y[P.cell], Uf->z[P.cell] };
        const util::Vec3d uRel{ uF.x - P.v.x, uF.y - P.v.y, uF.z - P.v.z };
        const double uRelMag = std::sqrt(dot3(uRel, uRel));
        const double Rep = rhoF_ * uRelMag * spec_.diameter / std::max(muF_, 1e-30);
        double Cd = 0.44;
        if (Rep < 1000.0 && Rep > 1e-12)
            Cd = 24.0 / Rep * (1.0 + 0.15 * std::pow(Rep, 0.687));

        // F_drag on particle (N): F = ½ ρ |u_rel| u_rel C_d A_p
        const double Ap = M_PI * 0.25 * spec_.diameter * spec_.diameter;
        const util::Vec3d Fd{
            0.5 * rhoF_ * uRelMag * uRel.x * Cd * Ap,
            0.5 * rhoF_ * uRelMag * uRel.y * Cd * Ap,
            0.5 * rhoF_ * uRelMag * uRel.z * Cd * Ap
        };

        const util::Vec3d Fg{ P.mass * spec_.gravity.x,
                              P.mass * spec_.gravity.y,
                              P.mass * spec_.gravity.z };

        // --- Saffman shear-induced lift (Mei 1992 correction of Saffman 1965) ---
        util::Vec3d Fs{0,0,0};
        if (spec_.saffmanLift && !vortF.empty()) {
            const util::Vec3d w = vortF[P.cell];
            const double wMag = std::sqrt(dot3(w, w));
            if (wMag > 1e-30) {
                // K_sl = 1.615 d_p² (ρ_f μ_f)^½ |ω|^½  (Saffman 1965)
                const double K = 1.615 * spec_.diameter * spec_.diameter
                                * std::sqrt(rhoF_ * muF_) / std::sqrt(wMag);
                // F_S = K (ω × u_rel)
                Fs = { K * (w.y * uRel.z - w.z * uRel.y),
                       K * (w.z * uRel.x - w.x * uRel.z),
                       K * (w.x * uRel.y - w.y * uRel.x) };
                // Mei (1992) Reynolds-number correction f(Re_p, Re_s)
                if (Rep > 0.1) {
                    const double Res = rhoF_ * wMag * spec_.diameter * spec_.diameter
                                     / std::max(muF_, 1e-30);
                    const double beta = 0.5 * Res / std::max(Rep, 1e-12);
                    const double f = (Rep <= 40.0)
                        ? (1.0 - 0.3314 * std::sqrt(beta)) * std::exp(-Rep * 0.1)
                          + 0.3314 * std::sqrt(beta)
                        : 0.0524 * std::sqrt(beta * Rep);
                    Fs.x *= f; Fs.y *= f; Fs.z *= f;
                }
            }
        }

        // --- Magnus rotation-induced lift ---
        util::Vec3d Fm{0,0,0};
        if (spec_.magnusLift && spec_.trackRotation) {
            // Ω_rel = ω_p - ½ ω_f
            const util::Vec3d wHalf = vortF.empty() ? util::Vec3d{0,0,0}
                : util::Vec3d{ 0.5*vortF[P.cell].x, 0.5*vortF[P.cell].y, 0.5*vortF[P.cell].z };
            const util::Vec3d Wr{ P.omega.x - wHalf.x,
                                  P.omega.y - wHalf.y,
                                  P.omega.z - wHalf.z };
            // F_M = (π/8) d_p³ ρ_f (Ω_rel × u_rel)
            const double K = (M_PI / 8.0) * std::pow(spec_.diameter, 3.0) * rhoF_;
            Fm = { K * (Wr.y * uRel.z - Wr.z * uRel.y),
                   K * (Wr.z * uRel.x - Wr.x * uRel.z),
                   K * (Wr.x * uRel.y - Wr.y * uRel.x) };
        }

        // Newton's second law, explicit Euler (small dt regime adequate for
        // micron-scale particles; RK4 hook reserved for next pass).
        const double invM = 1.0 / std::max(P.mass, 1e-30);
        P.v.x += dt * (Fd.x + Fg.x + Fs.x + Fm.x) * invM;
        P.v.y += dt * (Fd.y + Fg.y + Fs.y + Fm.y) * invM;
        P.v.z += dt * (Fd.z + Fg.z + Fs.z + Fm.z) * invM;

        // --- Particle rotation EoM ---
        // I_p dω_p/dt = -C_ω (π/64) ρ_f d_p^5 |Ω_rel| Ω_rel    with
        //   Ω_rel = ω_p - ½ ω_f and  C_ω from Dennis et al. correlation:
        //   C_ω = 12.9/√Re_ω + 128.4/Re_ω   (Re_ω = ρ_f |Ω_rel| d_p² / μ_f)
        if (spec_.trackRotation && P.inertia > 0.0) {
            const util::Vec3d wHalf = vortF.empty() ? util::Vec3d{0,0,0}
                : util::Vec3d{ 0.5*vortF[P.cell].x, 0.5*vortF[P.cell].y, 0.5*vortF[P.cell].z };
            const util::Vec3d Wr{ P.omega.x - wHalf.x,
                                  P.omega.y - wHalf.y,
                                  P.omega.z - wHalf.z };
            const double WrMag = std::sqrt(dot3(Wr, Wr));
            if (WrMag > 1e-30) {
                const double Reo = rhoF_ * WrMag * spec_.diameter * spec_.diameter
                                 / std::max(muF_, 1e-30);
                const double Cw = (Reo > 1e-6)
                    ? (12.9 / std::sqrt(Reo) + 128.4 / Reo)
                    : 1e6;
                const double K = Cw * (M_PI / 64.0) * rhoF_
                               * std::pow(spec_.diameter, 5.0) * WrMag;
                const double invI = 1.0 / P.inertia;
                P.omega.x -= dt * K * Wr.x * invI;
                P.omega.y -= dt * K * Wr.y * invI;
                P.omega.z -= dt * K * Wr.z * invI;
            }
        }

        P.x.x += dt * P.v.x;
        P.x.y += dt * P.v.y;
        P.x.z += dt * P.v.z;

        // Two-way coupling: reaction on the fluid (Newton's 3rd law).
        // Source per unit volume (N/m³). Accumulated per cell.
        const double invV = 1.0 / std::max(C.volume[P.cell], 1e-30);
        srcX_[P.cell] -= (Fd.x + Fs.x + Fm.x) * invV;
        srcY_[P.cell] -= (Fd.y + Fs.y + Fm.y) * invV;
        srcZ_[P.cell] -= (Fd.z + Fs.z + Fm.z) * invV;

        P.cell = locate(P.x, P.cell);
        if (P.cell == static_cast<meshing::CellId>(-1)) P.active = false;
    }
}

void LagrangianTracker::momentum_source(solver::VectorField& S) const {
    if (S.x.size() != srcX_.size()) S.resize(srcX_.size());
    std::copy(srcX_.begin(), srcX_.end(), S.x.begin());
    std::copy(srcY_.begin(), srcY_.end(), S.y.begin());
    std::copy(srcZ_.begin(), srcZ_.end(), S.z.begin());
}

std::size_t LagrangianTracker::live() const noexcept {
    std::size_t n = 0;
    for (const auto& p : parts_) if (p.active) ++n;
    return n;
}

}  // namespace simall::particles
