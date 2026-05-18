// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/TfcModel.cpp
// =============================================================================
#include "combustion/TfcModel.hpp"
#include "core/Logger.hpp"
#include "solver/LeastSquaresGradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion {

TfcModel::TfcModel(meshing::Mesh& mesh,
                   solver::FieldRegistry& fields,
                   solver::ILinearSolver& linear)
    : mesh_(mesh), F_(fields), lin_(linear) {}

void TfcModel::configure(TfcProps props,
                         const std::vector<solver::BoundarySpec>& bcs) {
    p_   = props;
    bcs_ = bcs;
    const std::size_t nC = mesh_.cells().size();
    F_.scalar("c", nC);
    F_.scalar("k", nC);
    F_.scalar("epsilon", nC);
    F_.scalar("S_progress",   nC);
    F_.scalar("S_combustion", nC);
    F_.scalar("rho_mix",      nC);
    F_.scalar("T_combust",    nC);
    F_.vector("grad_c",       nC);

    cEq_ = std::make_unique<solver::ScalarTransport>(mesh_, F_, lin_);
    cEq_->set_field("c");
    cEq_->set_density(p_.rho_u);            // unburnt reference density
    cEq_->set_diffusivity(p_.D_t);
    cEq_->set_source(std::string{"S_progress"});
    cEq_->set_urf(0.7);
    for (const auto& b : bcs_) {
        switch (b.type) {
            case solver::BCType::VelocityInlet:
            case solver::BCType::PressureInlet:
            case solver::BCType::MassFlowInlet:
                cEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
                break;
            default:
                cEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann,   0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Combustion",
        "TFC (Zimont) configured: A=", p_.A, " S_L=", p_.S_L,
        " (", nC, " cells)");
}

void TfcModel::step() {
    const std::size_t nC = mesh_.cells().size();
    const auto& c       = *F_.find_scalar("c");
    const auto& k       = *F_.find_scalar("k");
    const auto& epsilon = *F_.find_scalar("epsilon");
    auto& Sw            = *F_.find_scalar("S_progress");
    auto& Sq            = *F_.find_scalar("S_combustion");
    auto& rho           = *F_.find_scalar("rho_mix");
    auto& Tcomb         = *F_.find_scalar("T_combust");
    auto& gC            = *F_.find_vector("grad_c");

    // ∇c.
    solver::LeastSquaresGradient gradOp(mesh_);
    gradOp.evaluate(c, gC);

    const double T_u = 300.0;
    const double T_b = T_u + p_.H_combust / 1004.5;   // c_p ≈ air

    for (std::size_t i = 0; i < nC; ++i) {
        const double kk  = std::max(k[i], 1e-12);
        const double eps = std::max(epsilon[i], 1e-12);
        const double up  = std::sqrt(2.0 * kk / 3.0);
        const double Lt  = std::pow(kk, 1.5) / eps;
        const double Ut  = p_.A * std::pow(up, 0.75) * std::pow(p_.S_L, 0.5)
                          * std::pow(p_.alpha_u, -0.25) * std::pow(Lt, 0.25);
        const double magGc = std::sqrt(gC.x[i]*gC.x[i] + gC.y[i]*gC.y[i] + gC.z[i]*gC.z[i]);
        const double omega = p_.rho_u * Ut * magGc;
        Sw [i] = omega;
        Sq [i] = omega * p_.H_combust;
        // Density and temperature interpolated by Bray-Moss-Libby c blend.
        rho  [i] = p_.rho_u / (1.0 + c[i] * (p_.rho_u / std::max(p_.rho_b, 1e-30) - 1.0));
        Tcomb[i] = T_u + c[i] * (T_b - T_u);
    }
    cEq_->solve_iteration();
    // Clip c ∈ [0,1] after solve to guarantee BML closure validity.
    auto& cw = *F_.find_scalar("c");
    for (std::size_t i = 0; i < nC; ++i) cw[i] = std::clamp(cw[i], 0.0, 1.0);
}

}  // namespace simall::combustion
