// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/LeeMassTransfer.cpp
// =============================================================================
#include "multiphase/LeeMassTransfer.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase {

void LeeMassTransfer::initialize(const meshing::Mesh& mesh, const LeeProps& props) {
    mesh_ = &mesh; p_ = props;
}

double LeeMassTransfer::apply(solver::FieldRegistry& F) {
    if (!mesh_) return 0.0;
    const auto& C = mesh_->cells();
    const std::size_t nC = C.size();

    const auto* T  = F.find_scalar("T");
    auto*       av = F.find_scalar("alpha_v");
    if (!T || !av) return 0.0;
    auto& Sav   = F.scalar("S_alpha_v", nC);
    auto& Sm    = F.scalar("S_mass",    nC);
    auto& Sen   = F.scalar("S_energy",  nC);

    double total = 0.0;
    const double Tsat = p_.T_sat;
    for (std::size_t c = 0; c < nC; ++c) {
        const double alphaV = std::clamp((*av)[c], 0.0, 1.0);
        const double alphaL = 1.0 - alphaV;
        const double dT = (*T)[c] - Tsat;
        double mdot = 0.0;
        if (dT > 0.0) {
            mdot = +p_.c_evap * alphaL * p_.rho_liquid * (dT / Tsat);
        } else if (dT < 0.0) {
            mdot = -p_.c_cond * alphaV * p_.rho_vapor  * (-dT / Tsat);
        }
        Sav[c] = mdot / std::max(1e-30, p_.rho_vapor);
        Sm [c] = mdot;
        Sen[c] = -mdot * p_.L_vap;
        total += mdot * C.volume[c];
    }
    return total;
}

}  // namespace simall::multiphase
