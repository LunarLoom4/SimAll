// =============================================================================
// SimAll Beta - Heat Transfer Subsystem
// File   : src/heat_transfer/ConjugateHeatTransfer.cpp
// =============================================================================
#include "heat_transfer/ConjugateHeatTransfer.hpp"
#include "solver/LinearSolvers.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace simall::heat {

SolidRegion::SolidRegion(meshing::Mesh m, SolidMaterial mat, std::string name)
    : mesh_(std::move(m)), mat_(mat), name_(std::move(name)) {}

void SolidRegion::initialize(double T0) {
    const std::size_t nC = mesh_.cells().size();
    F_.scalar("T", nC);
    auto& T = *F_.find_scalar("T");
    std::fill(T.begin(), T.end(), T0);

    // U is needed by the energy equation (advection term); for solids it
    // is identically zero — register a zero vector field of the right size.
    auto& U = F_.vector("U", nC);
    std::fill(U.x.begin(), U.x.end(), 0.0);
    std::fill(U.y.begin(), U.y.end(), 0.0);
    std::fill(U.z.begin(), U.z.end(), 0.0);

    solver::LinearSolverConfig lc;
    lc.kind          = solver::LinearSolverKind::CG;
    lc.preconditioner= solver::PreconditionerKind::ILU;
    lc.tolerance     = 1e-8;
    lc.maxIterations = 200;
    lin_ = solver::make_linear_solver(lc);

    solver::EnergyOptions eo;
    eo.rho = mat_.rho; eo.cp = mat_.cp; eo.k = mat_.k; eo.PrT = 1.0; eo.urf = 1.0;
    // Default: all boundary zones adiabatic; coupler overrides interface zones.
    Teq_ = std::make_unique<solver::EnergyEquation>(mesh_, F_, *lin_, bcs_, eo);
    SIMALL_LOG_INFO("CHT", "Solid region '", name_, "' initialized: ", nC, " cells");
}

double SolidRegion::step(double /*dt*/) {
    return Teq_ ? Teq_->iterate() : 0.0;
}

void SolidRegion::set_interface_flux(meshing::ZoneId z, double q) {
    if (!Teq_) return;
    Teq_->set_zone_bc({z, solver::ScalarBC::Kind::Neumann, q, 0.0});
}

double SolidRegion::interface_temperature(meshing::ZoneId z) const {
    return Teq_ ? Teq_->mean_zone_temperature(z) : 0.0;
}

// ---------------------------------------------------------------------------
//  Partitioned Dirichlet-Neumann coupling driver.
// ---------------------------------------------------------------------------
void ConjugateHeatTransfer::step(double dt,
                                 solver::EnergyEquation& fluidEnergy,
                                 meshing::Mesh& /*fluidMesh*/,
                                 solver::FieldRegistry& /*fluidFields*/,
                                 double tol, int maxOuter) {
    // Per-interface convergence history (last solid temperature).
    std::vector<double> Ts_prev(ifaces_.size(), 0.0);
    for (std::size_t i = 0; i < ifaces_.size(); ++i) {
        for (auto& s : solids_)
            Ts_prev[i] = s->interface_temperature(ifaces_[i].solidZone);
    }

    for (int outer = 0; outer < maxOuter; ++outer) {
        // 1. Fluid side: Dirichlet on each interface = current solid mean T.
        for (std::size_t i = 0; i < ifaces_.size(); ++i) {
            double Ts = 0.0;
            for (auto& s : solids_) {
                const double v = s->interface_temperature(ifaces_[i].solidZone);
                if (v != 0.0) Ts = v;
            }
            fluidEnergy.set_zone_bc({ifaces_[i].fluidZone,
                solver::ScalarBC::Kind::Dirichlet, Ts, 0.0});
        }
        fluidEnergy.iterate();

        // 2. Solid side: Neumann flux = heat flux that just flowed out of the
        //    fluid (Newton-3rd: inflow into solid).
        for (std::size_t i = 0; i < ifaces_.size(); ++i) {
            const double q = fluidEnergy.mean_zone_heat_flux(ifaces_[i].fluidZone);
            for (auto& s : solids_)
                s->set_interface_flux(ifaces_[i].solidZone, q);
        }
        for (auto& s : solids_) s->step(dt);

        // 3. Convergence test on interface temperature drift.
        double drift = 0.0;
        for (std::size_t i = 0; i < ifaces_.size(); ++i) {
            double Ts = 0.0;
            for (auto& s : solids_) {
                const double v = s->interface_temperature(ifaces_[i].solidZone);
                if (v != 0.0) Ts = v;
            }
            drift = std::max(drift, std::abs(Ts - Ts_prev[i]));
            Ts_prev[i] = Ts;
        }
        SIMALL_LOG_INFO("CHT", "outer ", outer, "  ΔT_iface = ", drift, " K");
        if (drift < tol) {
            SIMALL_LOG_INFO("CHT", "converged in ", outer + 1, " outer iterations");
            return;
        }
    }
    SIMALL_LOG_WARN("CHT", "did NOT converge in ", maxOuter, " outer iterations");
}

}  // namespace simall::heat
