// =============================================================================
// SimAll Beta - Particles Subsystem
// File   : src/particles/ChargedParticleField.cpp
// =============================================================================
#include "particles/ChargedParticleField.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::particles {

void ChargedParticleField::initialize(const meshing::Mesh& mesh,
                                      ChargedParticleProps props) {
    mesh_ = &mesh;
    p_    = props;
    SIMALL_LOG_INFO("Particles",
        "Charged field init: q=", p_.charge_per_parcel,
        " useReg=", p_.useFieldRegistry,
        " mag=", p_.includeMagnetic);
}

double ChargedParticleField::apply(double dt, LagrangianTracker& tracker,
                                   solver::FieldRegistry& F) {
    if (!mesh_) return 0.0;
    const auto* Ex = p_.useFieldRegistry ? F.find_scalar("E_x") : nullptr;
    const auto* Ey = p_.useFieldRegistry ? F.find_scalar("E_y") : nullptr;
    const auto* Ez = p_.useFieldRegistry ? F.find_scalar("E_z") : nullptr;
    auto& parts = tracker.mutable_particles();
    double dvMax = 0.0;
    const std::size_t nC = mesh_->cells().size();
    for (auto& pt : parts) {
        if (!pt.active || pt.mass <= 0.0) continue;
        util::Vec3d E = p_.E_uniform;
        if (p_.useFieldRegistry && Ex && Ey && Ez && pt.cell >= 0
            && static_cast<std::size_t>(pt.cell) < nC) {
            E = { (*Ex)[pt.cell], (*Ey)[pt.cell], (*Ez)[pt.cell] };
        }
        const double qm  = p_.charge_per_parcel / std::max(pt.mass, 1e-30);
        util::Vec3d a{ qm * E.x, qm * E.y, qm * E.z };
        if (p_.includeMagnetic) {
            const auto& B = p_.B_uniform;
            // a += (q/m) (v × B)
            a.x += qm * (pt.v.y * B.z - pt.v.z * B.y);
            a.y += qm * (pt.v.z * B.x - pt.v.x * B.z);
            a.z += qm * (pt.v.x * B.y - pt.v.y * B.x);
        }
        const double dvx = a.x * dt, dvy = a.y * dt, dvz = a.z * dt;
        pt.v.x += dvx; pt.v.y += dvy; pt.v.z += dvz;
        const double dvMag = std::sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
        dvMax = std::max(dvMax, dvMag);
    }
    return dvMax;
}

}  // namespace simall::particles
