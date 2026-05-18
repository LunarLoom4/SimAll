// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/EnergyEquation.cpp
// =============================================================================
#include "solver/EnergyEquation.hpp"
#include "core/Logger.hpp"

#include <algorithm>

namespace simall::solver {

EnergyEquation::EnergyEquation(meshing::Mesh& m, FieldRegistry& f,
                               ILinearSolver& lin,
                               const std::vector<BoundarySpec>& bcs,
                               EnergyOptions opt)
    : mesh_(m), F_(f), opt_(opt) {
    const std::size_t nC = mesh_.cells().size();
    F_.scalar("T", nC);

    Teq_ = std::make_unique<ScalarTransport>(mesh_, F_, lin);
    Teq_->set_field("T");
    Teq_->set_density(opt_.rho * opt_.cp);     // ρ c_p multiplies storage term
    Teq_->set_diffusivity(opt_.k);             // k_eff updated each iteration
    Teq_->set_urf(opt_.urf);

    for (const auto& b : bcs) {
        switch (b.type) {
            case BCType::PressureInlet:
            case BCType::VelocityInlet:
            case BCType::MassFlowInlet:
                Teq_->add_bc({b.zone, ScalarBC::Kind::Dirichlet, b.scalarValue, 0.0});
                break;
            case BCType::NoSlipWall:
            case BCType::Wall:
            case BCType::MovingWall:
                // Wall: Dirichlet on T if scalarValue>0, otherwise adiabatic
                if (b.scalarValue > 0.0)
                    Teq_->add_bc({b.zone, ScalarBC::Kind::Dirichlet, b.scalarValue, 0.0});
                else
                    Teq_->add_bc({b.zone, ScalarBC::Kind::Neumann, 0.0, 0.0});
                break;
            default:
                Teq_->add_bc({b.zone, ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
}

double EnergyEquation::iterate() {
    // Refresh effective conductivity from current μ_t (if present).
    const std::size_t nC = mesh_.cells().size();
    auto* mut = F_.find_scalar("mut");
    double mutAvg = 0.0;
    if (mut) {
        for (std::size_t c = 0; c < nC; ++c) mutAvg += (*mut)[c];
        mutAvg /= std::max<std::size_t>(nC, 1);
    }
    const double kEff = opt_.k + opt_.cp * mutAvg / std::max(opt_.PrT, 1e-6);
    Teq_->set_diffusivity(kEff);
    return Teq_->solve_iteration();
}

double EnergyEquation::mean_zone_temperature(meshing::ZoneId z) const {
    const auto& F = mesh_.faces();
    const auto* T = const_cast<FieldRegistry&>(F_).find_scalar("T");
    if (!T) return 0.0;
    double sum = 0.0, area = 0.0;
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (F.boundaryZone[f] != z) continue;
        const meshing::CellId o = F.owner[f];
        const double a = std::sqrt(F.areaX[f]*F.areaX[f] + F.areaY[f]*F.areaY[f]
                                 + F.areaZ[f]*F.areaZ[f]);
        sum  += (*T)[o] * a;
        area += a;
    }
    return area > 0.0 ? sum / area : 0.0;
}

double EnergyEquation::mean_zone_heat_flux(meshing::ZoneId z) const {
    const auto& F = mesh_.faces();
    const auto& C = mesh_.cells();
    const auto* T = const_cast<FieldRegistry&>(F_).find_scalar("T");
    if (!T) return 0.0;
    auto* mut = const_cast<FieldRegistry&>(F_).find_scalar("mut");
    double sum = 0.0, area = 0.0;
    const double Tface = mean_zone_temperature(z);
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (F.boundaryZone[f] != z) continue;
        const meshing::CellId o = F.owner[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        const double Afmag = std::sqrt(Ax*Ax + Ay*Ay + Az*Az);
        const double dx = F.centroidX[f] - C.centroidX[o];
        const double dy = F.centroidY[f] - C.centroidY[o];
        const double dz = F.centroidZ[f] - C.centroidZ[o];
        const double dn = std::abs(dx*Ax + dy*Ay + dz*Az) / std::max(Afmag, 1e-30);
        const double kEff = opt_.k + (mut ? opt_.cp * (*mut)[o] / std::max(opt_.PrT, 1e-6) : 0.0);
        // q out of fluid = k · (T_owner - T_face)/dn  (heat flows down ∇T into wall when wall colder)
        const double q = kEff * ((*T)[o] - Tface) / std::max(dn, 1e-30);
        sum  += q * Afmag;
        area += Afmag;
    }
    return area > 0.0 ? sum / area : 0.0;
}

}  // namespace simall::solver
