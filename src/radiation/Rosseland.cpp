// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/Rosseland.cpp
// =============================================================================
#include "radiation/Rosseland.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::radiation {

namespace {
constexpr double kSigmaSB = 5.670374419e-8;
}

bool Rosseland::initialize(const meshing::Mesh& mesh,
                           solver::FieldRegistry& F,
                           const RosselandProps& props) {
    mesh_ = &mesh; F_ = &F; p_ = props;
    const std::size_t nC = mesh.cells().size();
    F.scalar("T",     nC);
    F.scalar("S_rad", nC);
    SIMALL_LOG_INFO("Radiation",
        "Rosseland init: κ=", p_.absorption,
        " σ_s=", p_.scattering, " n=", p_.refractiveN);
    return true;
}

double Rosseland::apply() {
    if (!mesh_ || !F_) return 0.0;
    const auto& C  = mesh_->cells();
    const auto& Ff = mesh_->faces();
    const std::size_t nC = C.size();
    const auto* T  = F_->find_scalar("T");
    auto*       Sr = F_->find_scalar("S_rad");
    if (!T || !Sr) return 0.0;

    const double beta = std::max(p_.absorption + p_.scattering, 1e-30);
    const double n2   = p_.refractiveN * p_.refractiveN;

    // Per-cell radiative conductivity k_R = 16 σ n² T³ / (3 β).
    std::vector<double> kR(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        const double Tc = std::clamp((*T)[c], p_.T_min, p_.T_max);
        kR[c] = 16.0 * kSigmaSB * n2 * Tc * Tc * Tc / (3.0 * beta);
    }

    // ∇·(k_R ∇T) by orthogonal face-flux divergence on interior faces.
    std::vector<double> divFlux(nC, 0.0);
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const auto o = Ff.owner[f];
        const auto n = Ff.neighbor[f];
        if (n == meshing::kBoundaryCell) continue;
        const double dx = C.centroidX[n] - C.centroidX[o];
        const double dy = C.centroidY[n] - C.centroidY[o];
        const double dz = C.centroidZ[n] - C.centroidZ[o];
        const double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
        const double aMag = std::sqrt(Ff.areaX[f]*Ff.areaX[f]
                                    + Ff.areaY[f]*Ff.areaY[f]
                                    + Ff.areaZ[f]*Ff.areaZ[f]);
        if (d < 1e-30 || aMag < 1e-30) continue;
        const double kFace = 0.5 * (kR[o] + kR[n]);
        const double gradT = ((*T)[n] - (*T)[o]) / d;
        const double flux  = kFace * gradT * aMag;       // W
        divFlux[o] += flux;                              // +k∇T·n out of o
        divFlux[n] -= flux;
    }

    // Convert flux per cell to volumetric source (W/m³) and write/add.
    double maxAbs = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        const double V = std::max(C.volume[c], 1e-30);
        const double s = divFlux[c] / V;
        if (p_.accumulate) (*Sr)[c] += s; else (*Sr)[c] = s;
        maxAbs = std::max(maxAbs, std::abs(s));
    }
    return maxAbs;
}

}  // namespace simall::radiation
