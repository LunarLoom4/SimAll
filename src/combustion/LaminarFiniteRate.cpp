// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/LaminarFiniteRate.cpp
// =============================================================================
#include "combustion/LaminarFiniteRate.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion {

namespace {
constexpr double Rgas = 8.31446261815324;   // [J/(mol·K)]
}

LaminarFiniteRate::LaminarFiniteRate(meshing::Mesh& m, solver::FieldRegistry& f,
                                     solver::ILinearSolver& l)
    : mesh_(m), F_(f), lin_(l) {}

void LaminarFiniteRate::initialize() {
    const std::size_t nC = mesh_.cells().size();
    F_.scalar("S_combustion", nC);
    Yeqs_.clear(); srcY_.clear();
    Yeqs_.reserve(species_.size());
    srcY_.reserve(species_.size());

    for (const auto& sp : species_) {
        const std::string Yname = "Y_" + sp.name;
        const std::string sName = "__src_" + sp.name;
        F_.scalar(Yname, nC);
        F_.scalar(sName, nC);
        auto eq = std::make_unique<solver::ScalarTransport>(mesh_, F_, lin_);
        eq->set_field(Yname);
        eq->set_density(rho_);
        eq->set_diffusivity(rho_ * sp.diffusivity);
        eq->set_source(*F_.find_scalar(sName));
        eq->set_urf(0.8);
        for (const auto& b : bcs_) {
            switch (b.type) {
                case solver::BCType::VelocityInlet:
                case solver::BCType::PressureInlet:
                case solver::BCType::MassFlowInlet:
                    eq->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet,
                                b.scalarValue, 0.0});
                    break;
                default:
                    eq->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            }
        }
        Yeqs_.push_back(std::move(eq));
        srcY_.emplace_back(nC, 0.0);
    }
    SIMALL_LOG_INFO("Combustion", "Initialized ", species_.size(),
        " species, ", reactions_.size(), " reactions");
}

void LaminarFiniteRate::step() {
    const std::size_t nC = mesh_.cells().size();
    const std::size_t nS = species_.size();
    const auto* T = F_.find_scalar("T");

    // Reset per-species sources and combined enthalpy release.
    auto& Scomb = *F_.find_scalar("S_combustion");
    std::fill(Scomb.begin(), Scomb.end(), 0.0);
    for (auto& v : srcY_) std::fill(v.begin(), v.end(), 0.0);

    // Fetch Y arrays.
    std::vector<util::aligned_vector<double>*> Y(nS, nullptr);
    for (std::size_t k = 0; k < nS; ++k)
        Y[k] = F_.find_scalar("Y_" + species_[k].name);

    // Evaluate reaction rates and accumulate ω̇_k.
    for (std::size_t c = 0; c < nC; ++c) {
        const double Tc = T ? std::max((*T)[c], 200.0) : 300.0;
        const double invRT = 1.0 / (Rgas * Tc);
        // Local molar concentrations [X_j] = ρ Y_j / W_j.
        std::vector<double> conc(nS);
        for (std::size_t k = 0; k < nS; ++k) {
            const double Yk = Y[k] ? std::max((*Y[k])[c], 0.0) : 0.0;
            conc[k] = rho_ * Yk / std::max(species_[k].molarMass, 1e-30);
        }
        for (const auto& r : reactions_) {
            double q = r.A * std::pow(Tc, r.beta) * std::exp(-r.Ea * invRT);
            for (std::size_t k = 0; k < nS; ++k) {
                const double ord = (k < r.order.size()) ? r.order[k] : r.nuR[k];
                if (ord > 0.0) q *= std::pow(conc[k], ord);
            }
            for (std::size_t k = 0; k < nS; ++k) {
                const double nuNet = (k < r.nuP.size() ? r.nuP[k] : 0.0)
                                   - (k < r.nuR.size() ? r.nuR[k] : 0.0);
                if (nuNet == 0.0) continue;
                const double wdot_k = nuNet * q * species_[k].molarMass;  // [kg/m³/s]
                srcY_[k][c] += wdot_k;
                Scomb[c]    -= species_[k].formationEnthalpy * wdot_k;
            }
        }
    }

    // Copy sources into the FieldRegistry scratch storage that ScalarTransport
    // already references via set_source (which captured the field pointer).
    for (std::size_t k = 0; k < nS; ++k) {
        auto& dst = *F_.find_scalar("__src_" + species_[k].name);
        std::copy(srcY_[k].begin(), srcY_[k].end(), dst.begin());
    }

    // Solve each species transport equation (one outer pass).
    for (auto& eq : Yeqs_) eq->solve_iteration();

    // Clamp & renormalise Y_k so Σ Y_k ≤ 1 (preserves chemistry-physics
    // consistency in early iterations).
    for (std::size_t c = 0; c < nC; ++c) {
        double sum = 0.0;
        for (std::size_t k = 0; k < nS; ++k) {
            (*Y[k])[c] = std::max((*Y[k])[c], 0.0);
            sum += (*Y[k])[c];
        }
        if (sum > 1.0) {
            const double s = 1.0 / sum;
            for (std::size_t k = 0; k < nS; ++k) (*Y[k])[c] *= s;
        }
    }
}

}  // namespace simall::combustion
