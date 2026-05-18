// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/TransportedPdf.cpp
// =============================================================================
#include "combustion/TransportedPdf.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace simall::combustion {

namespace {
constexpr double Ru = 8.314462618;     // J/(mol·K)
}

void TransportedPdf::initialize(const meshing::Mesh& mesh,
                                const std::vector<std::string>& speciesNames,
                                TransportedPdfParams params) {
    mesh_  = &mesh;
    p_     = params;
    names_ = speciesNames;
    rng_.seed(p_.rngSeed);
    const std::size_t nC = mesh.cells().size();
    cellParticles_.assign(nC, {});
    for (auto& cell : cellParticles_) {
        cell.resize(p_.particlesPerCell);
        for (auto& pt : cell) {
            pt.Y.assign(names_.size(), 0.0);
            if (!pt.Y.empty()) pt.Y[0] = 1.0;     // start in pure first species
        }
    }
    SIMALL_LOG_INFO("Combustion",
        "TransportedPdf initialised: ", nC, " cells × ",
        p_.particlesPerCell, " particles × ", names_.size(), " species");
}

void TransportedPdf::apply_iem(double dt, std::vector<PdfParticle>& particles,
                               const std::vector<double>& Ybar) {
    // IEM (Interaction by Exchange with the Mean / Linear Mean Square Estimation):
    //   dψ_i/dt = -(C_phi / (2 τ_mix)) (ψ_i - ⟨ψ⟩)
    // Here τ_mix is folded into the time step; use C_phi/2 as decay rate.
    const double alpha = 1.0 - std::exp(-p_.C_phi * 0.5 * dt);
    for (auto& pt : particles) {
        for (std::size_t k = 0; k < pt.Y.size(); ++k) {
            pt.Y[k] += alpha * (Ybar[k] - pt.Y[k]);
        }
    }
}

void TransportedPdf::apply_modified_curl(double dt,
                                         std::vector<PdfParticle>& particles) {
    // Modified Curl (Janicka et al. 1979): pick Np * (3/2) * C_phi * dt random
    // pairs, partially mix each pair toward their midpoint with a uniform
    // random weight a ∈ [0,1].
    const std::size_t Np = particles.size();
    if (Np < 2) return;
    const double rate   = 1.5 * p_.C_phi * dt;
    const auto   nPairs = static_cast<std::size_t>(std::max(0.0, rate * Np));
    std::uniform_int_distribution<std::size_t> uidx(0, Np - 1);
    std::uniform_real_distribution<double>     ua(0.0, 1.0);
    for (std::size_t m = 0; m < nPairs; ++m) {
        std::size_t i = uidx(rng_), j = uidx(rng_);
        if (i == j) continue;
        const double a = ua(rng_);
        for (std::size_t k = 0; k < particles[i].Y.size(); ++k) {
            const double mid = 0.5 * (particles[i].Y[k] + particles[j].Y[k]);
            particles[i].Y[k] += a * (mid - particles[i].Y[k]);
            particles[j].Y[k] += a * (mid - particles[j].Y[k]);
        }
    }
}

void TransportedPdf::apply_chemistry(double dt,
                                     std::vector<PdfParticle>& particles,
                                     const std::vector<ChemkinReaction>& reactions) {
    if (reactions.empty()) return;
    // Build species index map.
    std::unordered_map<std::string, std::size_t> sidx;
    for (std::size_t i = 0; i < names_.size(); ++i) sidx[names_[i]] = i;

    for (auto& pt : particles) {
        std::vector<double> dY(pt.Y.size(), 0.0);
        for (const auto& r : reactions) {
            // Arrhenius rate at T_ref (linearised): k = A T^β exp(-Ea/(R T))
            const double T = std::max(p_.T_reference, 300.0);
            const double k = r.fwd.A * std::pow(T, r.fwd.beta)
                           * std::exp(-r.fwd.Ea / (Ru * T));
            // Forward product of [X_j]^{order} with concentrations from Y / Mw.
            double q = k;
            for (const auto& [name, nu] : r.reactants) {
                auto it = sidx.find(name);
                if (it == sidx.end()) { q = 0; break; }
                q *= std::pow(std::max(pt.Y[it->second], 0.0), nu);
            }
            for (const auto& [name, nu] : r.reactants) {
                auto it = sidx.find(name);
                if (it != sidx.end()) dY[it->second] -= nu * q * dt;
            }
            for (const auto& [name, nu] : r.products) {
                auto it = sidx.find(name);
                if (it != sidx.end()) dY[it->second] += nu * q * dt;
            }
        }
        for (std::size_t k = 0; k < pt.Y.size(); ++k) {
            pt.Y[k] = std::clamp(pt.Y[k] + dY[k], 0.0, 1.0);
        }
        // Renormalise Σ Y_k = 1.
        double s = 0.0;
        for (double y : pt.Y) s += y;
        if (s > 1e-30) for (double& y : pt.Y) y /= s;
    }
}

void TransportedPdf::step(double dt, solver::FieldRegistry& F,
                          const std::vector<ChemkinReaction>& reactions) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    if (cellParticles_.size() != nC) return;

    // Pre-allocate field channels for each species.
    std::vector<solver::ScalarField*> ybar(names_.size(), nullptr);
    for (std::size_t k = 0; k < names_.size(); ++k)
        ybar[k] = &F.scalar("Y_" + names_[k], nC);

    std::normal_distribution<double> gauss(0.0, 1.0);
    for (std::size_t c = 0; c < nC; ++c) {
        auto& parts = cellParticles_[c];
        // 1. Ensemble mean before mixing.
        std::vector<double> mean(names_.size(), 0.0);
        for (const auto& pt : parts)
            for (std::size_t k = 0; k < pt.Y.size(); ++k)
                mean[k] += pt.Y[k];
        const double invN = 1.0 / static_cast<double>(parts.size());
        for (auto& v : mean) v *= invN;

        // 2. Random walk on composition (Wiener increment in Y-space) — small
        // amplitude proportional to √dt × σ_Y; σ_Y from cell-mean variance.
        double sigma = 0.0;
        for (const auto& pt : parts)
            for (std::size_t k = 0; k < pt.Y.size(); ++k)
                sigma += (pt.Y[k] - mean[k]) * (pt.Y[k] - mean[k]);
        sigma = std::sqrt(sigma / std::max<double>(parts.size() * names_.size(), 1));
        const double walkAmpl = sigma * std::sqrt(dt);
        for (auto& pt : parts) {
            for (auto& y : pt.Y) {
                y = std::clamp(y + walkAmpl * gauss(rng_), 0.0, 1.0);
            }
            pt.age += dt;
        }

        // 3. Mixing step.
        if (p_.mixing == PdfMixingModel::IEM) apply_iem(dt, parts, mean);
        else                                   apply_modified_curl(dt, parts);

        // 4. Chemistry step.
        apply_chemistry(dt, parts, reactions);

        // 5. Re-average and write back to fields.
        std::fill(mean.begin(), mean.end(), 0.0);
        for (const auto& pt : parts)
            for (std::size_t k = 0; k < pt.Y.size(); ++k)
                mean[k] += pt.Y[k];
        for (auto& v : mean) v *= invN;
        for (std::size_t k = 0; k < names_.size(); ++k) (*ybar[k])[c] = mean[k];
    }
}

}  // namespace simall::combustion
