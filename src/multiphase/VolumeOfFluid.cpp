// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/VolumeOfFluid.cpp
// =============================================================================
#include "multiphase/VolumeOfFluid.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase {

VolumeOfFluid::VolumeOfFluid(meshing::Mesh& m, solver::FieldRegistry& f,
                             solver::ILinearSolver& l)
    : mesh_(m), F_(f), lin_(l) {
    const std::size_t nC = mesh_.cells().size();
    F_.scalar("alpha", nC);
    F_.scalar("rho",   nC);
    F_.scalar("mu",    nC);
    F_.vector("f_surfTension", nC);
    alphaEq_ = std::make_unique<solver::ScalarTransport>(mesh_, F_, lin_);
    alphaEq_->set_field("alpha");
    alphaEq_->set_density(1.0);    // ρ already baked into alpha equation
    alphaEq_->set_diffusivity(0.0);
}

void VolumeOfFluid::configure(VOFFluid p, VOFFluid s,
                              const std::vector<solver::BoundarySpec>& bcs) {
    f1_ = p; f2_ = s; bcs_ = bcs;
    for (const auto& b : bcs_) {
        // Inlets carry their own α (from b.scalarValue), walls/outlets are
        // zero-gradient.
        if (b.type == solver::BCType::VelocityInlet ||
            b.type == solver::BCType::PressureInlet ||
            b.type == solver::BCType::MassFlowInlet) {
            alphaEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet,
                              b.scalarValue, 0.0});
        } else {
            alphaEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
}

void VolumeOfFluid::step() {
    alphaEq_->solve_iteration();
    auto& a   = *F_.find_scalar("alpha");
    auto& rho = *F_.find_scalar("rho");
    auto& mu  = *F_.find_scalar("mu");
    const std::size_t nC = a.size();
    for (std::size_t c = 0; c < nC; ++c) {
        a[c]   = std::clamp(a[c], 0.0, 1.0);
        rho[c] = a[c] * f1_.density   + (1 - a[c]) * f2_.density;
        mu[c]  = a[c] * f1_.viscosity + (1 - a[c]) * f2_.viscosity;
    }

    // ----- Continuum Surface Force (Brackbill, Kothe & Zemach 1992) -------
    // f_st = σ κ ∇α,   κ = -∇·n̂,   n̂ = ∇α/|∇α|
    auto& fst = *F_.find_vector("f_surfTension");
    std::fill(fst.x.begin(), fst.x.end(), 0.0);
    std::fill(fst.y.begin(), fst.y.end(), 0.0);
    std::fill(fst.z.begin(), fst.z.end(), 0.0);
    if (sigma_ <= 0.0) return;

    solver::LeastSquaresGradient G(mesh_);
    solver::VectorField gA;
    G.evaluate(a, gA);

    util::aligned_vector<double> nx(nC), ny(nC), nz(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        const double mag = std::sqrt(gA.x[c]*gA.x[c] + gA.y[c]*gA.y[c] + gA.z[c]*gA.z[c]);
        const double inv = (mag > 1e-12) ? 1.0 / mag : 0.0;
        nx[c] = gA.x[c] * inv;
        ny[c] = gA.y[c] * inv;
        nz[c] = gA.z[c] * inv;
    }

    solver::VectorField gNx, gNy, gNz;
    G.evaluate(nx, gNx); G.evaluate(ny, gNy); G.evaluate(nz, gNz);
    for (std::size_t c = 0; c < nC; ++c) {
        const double kappa = -(gNx.x[c] + gNy.y[c] + gNz.z[c]);
        fst.x[c] = sigma_ * kappa * gA.x[c];
        fst.y[c] = sigma_ * kappa * gA.y[c];
        fst.z[c] = sigma_ * kappa * gA.z[c];
    }
}

}  // namespace simall::multiphase
