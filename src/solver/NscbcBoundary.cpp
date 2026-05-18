// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/NscbcBoundary.cpp
//
// Implementation of LODI (Local One-Dimensional Inviscid) characteristic
// boundary updates. For each NSCBC zone we sweep all owned boundary faces
// and accumulate the conservative-variable time derivative
//   ∂U/∂t = -(d1, d2, d3, d4, d5) · ∂U/∂L
// where d_i are the standard LODI vectors (Poinsot & Lele 1992 eq. 17).
// The derivative ∂φ/∂x_n is approximated by (φ_face − φ_cell) / (n·(x_f-x_c)).
// =============================================================================
#include "solver/NscbcBoundary.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver {

namespace {
const NscbcZoneSpec* lookup(const std::vector<NscbcZoneSpec>& v, meshing::ZoneId z) {
    for (const auto& s : v) if (s.zone == z) return &s;
    return nullptr;
}
}

void NscbcBoundary::initialize(const meshing::Mesh& m, std::vector<NscbcZoneSpec> z,
                               double gamma, double Rgas) {
    mesh_ = &m; zones_ = std::move(z); gamma_ = gamma; Rgas_ = Rgas;
}

std::size_t NscbcBoundary::apply(FieldRegistry& F) {
    if (!mesh_) return 0;
    auto* rho  = F.find_scalar("rho");
    auto* p    = F.find_scalar("p");
    auto* T    = F.find_scalar("T");
    auto* U    = F.find_vector("U");
    auto* drho = F.find_scalar("drho_nscbc");
    auto* drhU = F.find_vector("drhoU_nscbc");
    auto* drhE = F.find_scalar("drhoE_nscbc");
    if (!rho || !p || !T || !U) return 0;
    const std::size_t nC = mesh_->cells().size();
    if (!drho) { drho = &F.scalar("drho_nscbc", nC); }
    if (!drhU) { drhU = &F.vector("drhoU_nscbc", nC); }
    if (!drhE) { drhE = &F.scalar("drhoE_nscbc", nC); }

    const auto& C = mesh_->cells();
    const auto& Ff = mesh_->faces();
    std::size_t hits = 0;

    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.neighbor[f] != meshing::kBoundaryCell) continue;
        const auto* spec = lookup(zones_, Ff.boundaryZone[f]);
        if (!spec) continue;
        const meshing::CellId c = Ff.owner[f];

        util::Vec3d n{Ff.areaX[f], Ff.areaY[f], Ff.areaZ[f]};
        const double a = n.norm(); if (a < 1e-30) continue;
        n.x /= a; n.y /= a; n.z /= a;
        const double sx = Ff.centroidX[f] - C.centroidX[c];
        const double sy = Ff.centroidY[f] - C.centroidY[c];
        const double sz = Ff.centroidZ[f] - C.centroidZ[c];
        if (n.x*sx + n.y*sy + n.z*sz < 0) { n.x = -n.x; n.y = -n.y; n.z = -n.z; }
        const double dn = std::max(1e-30, n.x*sx + n.y*sy + n.z*sz);

        const double rhoC = (*rho)[c];
        const double pC   = (*p)[c];
        const double TC   = (*T)[c];
        const double uxC  = U->x[c], uyC = U->y[c], uzC = U->z[c];
        const double cSnd = std::sqrt(std::max(1e-12, gamma_ * Rgas_ * TC));
        const double un   = uxC*n.x + uyC*n.y + uzC*n.z;
        const double M    = std::abs(un) / cSnd;

        // Tangent basis.
        util::Vec3d t1 = (std::abs(n.x) < 0.9) ? util::Vec3d{1,0,0} : util::Vec3d{0,1,0};
        const double tdot = t1.x*n.x + t1.y*n.y + t1.z*n.z;
        t1.x -= tdot*n.x; t1.y -= tdot*n.y; t1.z -= tdot*n.z;
        const double tn = t1.norm(); t1.x /= tn; t1.y /= tn; t1.z /= tn;
        const util::Vec3d t2 = n.cross(t1);

        const double ut1 = uxC*t1.x + uyC*t1.y + uzC*t1.z;
        const double ut2 = uxC*t2.x + uyC*t2.y + uzC*t2.z;

        // Face values (target).
        double rhoF = rhoC, pF = pC, unF = un, ut1F = ut1, ut2F = ut2, TF = TC;
        switch (spec->type) {
            case NscbcType::SubsonicOutflow:
                pF = spec->pInf;
                break;
            case NscbcType::SubsonicInflow: {
                TF   = spec->TInf;
                const util::Vec3d ui = spec->uInf;
                unF  = ui.x*n.x + ui.y*n.y + ui.z*n.z;
                ut1F = ui.x*t1.x + ui.y*t1.y + ui.z*t1.z;
                ut2F = ui.x*t2.x + ui.y*t2.y + ui.z*t2.z;
                rhoF = pC / (Rgas_ * TF);
                break;
            }
            case NscbcType::NonReflectingWall:
                unF = 0.0;
                break;
        }

        // Gradients (one-sided towards face).
        const double dp_dn   = (pF   - pC)   / dn;
        const double drho_dn = (rhoF - rhoC) / dn;
        const double dun_dn  = (unF  - un)   / dn;
        const double dut1_dn = (ut1F - ut1)  / dn;
        const double dut2_dn = (ut2F - ut2)  / dn;

        // Characteristic amplitudes.
        const double l1c = (un - cSnd);
        const double l5c = (un + cSnd);
        double L1 = l1c * (dp_dn - rhoC*cSnd*dun_dn);
        double L2 = un  * (cSnd*cSnd * drho_dn - dp_dn);
        double L3 = un  * dut1_dn;
        double L4 = un  * dut2_dn;
        double L5 = l5c * (dp_dn + rhoC*cSnd*dun_dn);

        // LODI relaxation for incoming waves.
        const double K = spec->sigma * (1.0 - M*M) * cSnd / std::max(1e-12, spec->length);
        if (spec->type == NscbcType::SubsonicOutflow) {
            if (un > 0)        L1 = K * (pC - spec->pInf);   // incoming
            else                L5 = K * (pC - spec->pInf);
        } else if (spec->type == NscbcType::SubsonicInflow) {
            // Inflow: 4 incoming waves to specify (un, ut1, ut2, T).
            if (un > 0) {
                L5 = K * rhoC*cSnd * (un - unF);
                L3 = K * (ut1 - ut1F);
                L4 = K * (ut2 - ut2F);
                L2 = K * (TC - TF) * rhoC * cSnd*cSnd / TC;
            } else {
                L1 = K * rhoC*cSnd * (unF - un);
            }
        }

        // Conservative-variable derivatives (d1..d5).
        const double d1 = (L2 + 0.5*(L5 + L1)) / (cSnd*cSnd);
        const double d2 = 0.5*(L5 + L1);
        const double d3 = 0.5*(L5 - L1) / (rhoC*cSnd);
        const double d4 = L3;
        const double d5 = L4;

        (*drho)[c]  -= d1;
        const double dmx = d1*uxC + rhoC*(d3*n.x + d4*t1.x + d5*t2.x);
        const double dmy = d1*uyC + rhoC*(d3*n.y + d4*t1.y + d5*t2.y);
        const double dmz = d1*uzC + rhoC*(d3*n.z + d4*t1.z + d5*t2.z);
        drhU->x[c] -= dmx;
        drhU->y[c] -= dmy;
        drhU->z[c] -= dmz;
        const double H = (gamma_*pC/((gamma_-1.0)*rhoC))
                       + 0.5*(uxC*uxC + uyC*uyC + uzC*uzC);
        const double dE = 0.5*(uxC*uxC + uyC*uyC + uzC*uzC)*d1
                        + d2/(gamma_-1.0)
                        + rhoC*(uxC*(d3*n.x+d4*t1.x+d5*t2.x)
                              + uyC*(d3*n.y+d4*t1.y+d5*t2.y)
                              + uzC*(d3*n.z+d4*t1.z+d5*t2.z));
        (*drhE)[c] -= dE;
        (void)H;
        ++hits;
    }
    return hits;
}

}  // namespace simall::solver
