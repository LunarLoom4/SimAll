// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/DemCollision.cpp
// Implementation of Hertz-Mindlin DEM contact. Uniform spatial hash grid
// for particle-particle pair enumeration, mesh-face traversal for
// particle-wall contact.
// =============================================================================
#include "particles/DemCollision.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles {

namespace {
inline double sgn(double x) { return (x > 0) - (x < 0); }
}

void DemCollision::initialize(const meshing::Mesh& mesh, const DemMaterial& mat) {
    mesh_ = &mesh; mat_  = mat;
}

void DemCollision::contact_pair(double dt, ParticleState& a, ParticleState& b,
                                double dA, double dB, double mA, double mB) {
    const util::Vec3d r = b.x - a.x;
    const double dist = r.norm();
    const double rsum = 0.5*(dA + dB);
    if (dist >= rsum || dist < 1e-18) return;

    const util::Vec3d n = r * (1.0 / dist);
    const double delta = rsum - dist;                // overlap

    // Effective properties (Hertz).
    const double Eeff = mat_.youngsModulus
                       / (2.0 * (1.0 - mat_.poissonRatio*mat_.poissonRatio));
    const double Reff = (dA*dB) / (2.0*(dA + dB));
    const double mEff = (mA*mB) / (mA + mB);

    const double kn = (4.0/3.0) * Eeff * std::sqrt(Reff);
    const double Fn_spring = kn * std::pow(delta, 1.5);

    // Critical damping from coefficient of restitution e:
    //   ln(e) / sqrt(π² + ln²(e))
    const double le = std::log(std::max(mat_.restitution, 1e-6));
    const double beta = -le / std::sqrt(M_PI*M_PI + le*le);
    const double knLin = 2.0 * Eeff * std::sqrt(Reff * delta);   // linearised
    const double cn = 2.0 * beta * std::sqrt(mEff * knLin);

    const util::Vec3d vrel = b.v - a.v;
    const double vrn = vrel.dot(n);
    const util::Vec3d vt{vrel.x - vrn*n.x, vrel.y - vrn*n.y, vrel.z - vrn*n.z};

    const double Fn = Fn_spring + cn * vrn;          // attractive damping ok

    // Tangential: Mindlin no-slip stiffness ≈ 8 G* sqrt(R* δ); use Coulomb cap.
    const double Geff = Eeff / (2.0*(1.0 + mat_.poissonRatio));
    const double kt = 8.0 * Geff * std::sqrt(Reff * delta);
    const double vtMag = vt.norm();
    const double Ft_cap = mat_.friction * std::abs(Fn);
    const double Ft_lin = kt * vtMag * dt;
    const double Ft = std::min(Ft_lin, Ft_cap);
    util::Vec3d Ftvec{0,0,0};
    if (vtMag > 1e-15) {
        const double s = -Ft / vtMag;
        Ftvec = util::Vec3d{vt.x*s, vt.y*s, vt.z*s};
    }
    const util::Vec3d Fnvec{-Fn * n.x, -Fn * n.y, -Fn * n.z};
    const util::Vec3d Ftot = Fnvec + Ftvec;

    const double aFac = dt / mA;
    const double bFac = dt / mB;
    a.v.x -= Ftot.x * aFac; a.v.y -= Ftot.y * aFac; a.v.z -= Ftot.z * aFac;
    b.v.x += Ftot.x * bFac; b.v.y += Ftot.y * bFac; b.v.z += Ftot.z * bFac;
    ++lastContacts_;
}

void DemCollision::contact_wall(double dt, ParticleState& p, double dP, double mP,
                                const util::Vec3d& nWall, double pen) {
    if (pen <= 0.0) return;
    const double Eeff = mat_.youngsModulus / (1.0 - mat_.poissonRatio*mat_.poissonRatio);
    const double Reff = 0.5*dP;
    const double kn = (4.0/3.0) * Eeff * std::sqrt(Reff);
    const double Fn_spring = kn * std::pow(pen, 1.5);
    const double le = std::log(std::max(mat_.restitution, 1e-6));
    const double beta = -le / std::sqrt(M_PI*M_PI + le*le);
    const double knLin = 2.0 * Eeff * std::sqrt(Reff * pen);
    const double cn = 2.0 * beta * std::sqrt(mP * knLin);
    const double vn = p.v.dot(nWall);
    const double Fn = Fn_spring + cn * std::max(0.0, vn);

    // Tangential.
    const util::Vec3d vt{p.v.x - vn*nWall.x, p.v.y - vn*nWall.y, p.v.z - vn*nWall.z};
    const double vtMag = vt.norm();
    const double Geff = Eeff / (2.0*(1.0 + mat_.poissonRatio));
    const double kt = 8.0 * Geff * std::sqrt(Reff * pen);
    const double Ft_cap = mat_.friction * std::abs(Fn);
    const double Ft = std::min(kt * vtMag * dt, Ft_cap);
    util::Vec3d Ftvec{0,0,0};
    if (vtMag > 1e-15) {
        const double s = -Ft / vtMag;
        Ftvec = util::Vec3d{vt.x*s, vt.y*s, vt.z*s};
    }
    const util::Vec3d Ftot{-Fn*nWall.x + Ftvec.x,
                           -Fn*nWall.y + Ftvec.y,
                           -Fn*nWall.z + Ftvec.z};
    const double f = dt / mP;
    p.v.x += Ftot.x * f; p.v.y += Ftot.y * f; p.v.z += Ftot.z * f;
    // Push out of penetration to avoid sticking on next sub-step.
    p.x.x += nWall.x * pen; p.x.y += nWall.y * pen; p.x.z += nWall.z * pen;
    ++lastContacts_;
}

void DemCollision::apply(double dt, std::vector<ParticleState>& parts,
                         const ParticleSpec& spec) {
    lastContacts_ = 0;
    if (parts.empty() || !mesh_) return;

    // --- Particle-particle pair search via uniform spatial hash. ---
    const double cell = std::max(1e-9, spec.diameter * 1.05);
    std::unordered_map<HashKey, std::vector<std::size_t>, HashHasher> grid;
    grid.reserve(parts.size() * 2);
    for (std::size_t i = 0; i < parts.size(); ++i) {
        if (!parts[i].active) continue;
        HashKey k{
            static_cast<std::int32_t>(std::floor(parts[i].x.x / cell)),
            static_cast<std::int32_t>(std::floor(parts[i].x.y / cell)),
            static_cast<std::int32_t>(std::floor(parts[i].x.z / cell))};
        grid[k].push_back(i);
    }

    auto pairWith = [&](std::size_t i, std::size_t j) {
        if (j <= i) return;
        contact_pair(dt, parts[i], parts[j],
                     spec.diameter, spec.diameter,
                     parts[i].mass, parts[j].mass);
    };

    for (const auto& [k, bucket] : grid) {
        for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj)
        for (int dk = -1; dk <= 1; ++dk) {
            HashKey n{k.i+di, k.j+dj, k.k+dk};
            auto it = grid.find(n);
            if (it == grid.end()) continue;
            for (std::size_t i : bucket) for (std::size_t j : it->second) {
                pairWith(i, j);
            }
        }
    }

    // --- Particle-wall contact: scan host cell's boundary faces. ---
    const auto& C = mesh_->cells();
    const auto& F = mesh_->faces();
    const double radius = 0.5 * spec.diameter;
    for (auto& p : parts) {
        if (!p.active) continue;
        if (p.cell == static_cast<meshing::CellId>(-1)) continue;
        const int s = C.faceOffsets[p.cell], e = C.faceOffsets[p.cell+1];
        for (int kf = s; kf < e; ++kf) {
            const meshing::FaceId fid = C.faceIndices[kf];
            if (F.neighbor[fid] != meshing::kBoundaryCell) continue;
            util::Vec3d n{F.areaX[fid], F.areaY[fid], F.areaZ[fid]};
            const double a = n.norm(); if (a < 1e-30) continue;
            n.x /= a; n.y /= a; n.z /= a;
            // Outward normal from cell centroid sign convention.
            const double sx = F.centroidX[fid] - C.centroidX[p.cell];
            const double sy = F.centroidY[fid] - C.centroidY[p.cell];
            const double sz = F.centroidZ[fid] - C.centroidZ[p.cell];
            if (n.x*sx + n.y*sy + n.z*sz < 0) { n.x = -n.x; n.y = -n.y; n.z = -n.z; }
            const double dn = (p.x.x - F.centroidX[fid])*n.x
                            + (p.x.y - F.centroidY[fid])*n.y
                            + (p.x.z - F.centroidZ[fid])*n.z;
            const double pen = radius - dn;
            if (pen > 0) {
                // Reverse normal to point into the fluid (away from wall).
                util::Vec3d nIn{-n.x, -n.y, -n.z};
                contact_wall(dt, p, spec.diameter, p.mass, nIn, pen);
            }
        }
    }
    (void)sgn;
}

}  // namespace simall::particles
