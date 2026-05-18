// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/SlfmFlamelet.cpp
// =============================================================================
#include "combustion/SlfmFlamelet.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion {

namespace {

std::size_t lower_index(const std::vector<double>& x, double v) {
    if (v <= x.front()) return 0;
    if (v >= x.back())  return x.size() - 2;
    auto it = std::upper_bound(x.begin(), x.end(), v);
    return static_cast<std::size_t>(std::distance(x.begin(), it)) - 1;
}

}  // namespace

void SlfmFlamelet::build_analytic(std::size_t NZ, std::size_t NlogChi,
                                  double Z_st, double T_un, double T_ad,
                                  double chi_q,
                                  double log10chi_min, double log10chi_max) {
    SlfmTable t;
    t.axis.Z.resize(NZ);
    t.axis.log10Chi.resize(NlogChi);
    for (std::size_t i = 0; i < NZ; ++i) t.axis.Z[i] = double(i) / double(NZ - 1);
    for (std::size_t j = 0; j < NlogChi; ++j)
        t.axis.log10Chi[j] = log10chi_min
                           + (log10chi_max - log10chi_min) * double(j) / double(NlogChi - 1);

    const double Hc = 50.0e6;     // J/kg fuel LHV (CH4)
    std::vector<double> T(NZ * NlogChi), rho(NZ * NlogChi), Q(NZ * NlogChi);

    for (std::size_t j = 0; j < NlogChi; ++j) {
        const double chi  = std::pow(10.0, t.axis.log10Chi[j]);
        const double dext = std::exp(-chi / std::max(chi_q, 1e-30));   // 1 → unstrained, →0 at extinction
        for (std::size_t i = 0; i < NZ; ++i) {
            const double Z   = t.axis.Z[i];
            const double Teq = (Z <= Z_st)
                ? T_un + (T_ad - T_un) * Z / std::max(Z_st, 1e-12)
                : T_un + (T_ad - T_un) * std::max(0.0,
                                                  1.0 - (Z - Z_st) / std::max(1.0 - Z_st, 1e-12));
            const double Tcell = T_un + dext * (Teq - T_un);
            const double rhocell = 353.0 / std::max(Tcell, 1e-3);  // ideal-air ρ ≈ p/(R T)
            const double Q_dot = dext * Hc * 0.5 * Z * (1.0 - Z);
            T  [i * NlogChi + j] = Tcell;
            rho[i * NlogChi + j] = rhocell;
            Q  [i * NlogChi + j] = Q_dot;
        }
    }
    t.var["T"]            = std::move(T);
    t.var["rho"]          = std::move(rho);
    t.var["S_combustion"] = std::move(Q);
    table_ = std::move(t);
    SIMALL_LOG_INFO("Combustion",
        "SLFM analytic table built: ", NZ, " × ", NlogChi);
}

double SlfmFlamelet::lookup(const std::string& var, double Z, double chi) const {
    auto it = table_.var.find(var);
    if (it == table_.var.end()) return 0.0;
    const auto& Zaxis = table_.axis.Z;
    const auto& Caxis = table_.axis.log10Chi;
    if (Zaxis.empty() || Caxis.empty()) return 0.0;
    const std::size_t Nc = Caxis.size();
    const double logChi = std::log10(std::max(chi, 1e-30));
    const std::size_t i0 = lower_index(Zaxis, Z);
    const std::size_t j0 = lower_index(Caxis, logChi);
    const double tz = std::clamp((Z - Zaxis[i0]) / std::max(Zaxis[i0+1]-Zaxis[i0], 1e-30), 0.0, 1.0);
    const double tc = std::clamp((logChi - Caxis[j0]) / std::max(Caxis[j0+1]-Caxis[j0], 1e-30), 0.0, 1.0);
    const auto& V = it->second;
    const double v00 = V[i0       * Nc + j0    ];
    const double v10 = V[(i0 + 1) * Nc + j0    ];
    const double v01 = V[i0       * Nc + j0 + 1];
    const double v11 = V[(i0 + 1) * Nc + j0 + 1];
    return (1 - tz) * (1 - tc) * v00 + tz * (1 - tc) * v10
         + (1 - tz) *      tc  * v01 + tz *      tc  * v11;
}

void SlfmFlamelet::initialize(const meshing::Mesh& mesh, solver::FieldRegistry& F) {
    mesh_ = &mesh;
    const std::size_t nC = mesh.cells().size();
    F.scalar("Z", nC);
    F.scalar("chi", nC);
    F.scalar("T", nC);
    F.scalar("rho", nC);
    F.scalar("S_combustion", nC);
    if (table_.var.empty()) build_analytic();
    SIMALL_LOG_INFO("Combustion", "SLFM initialised on ", nC, " cells");
}

void SlfmFlamelet::apply(solver::FieldRegistry& F) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& Z   = *F.find_scalar("Z");
    const auto& chi = *F.find_scalar("chi");
    auto& T   = *F.find_scalar("T");
    auto& rho = *F.find_scalar("rho");
    auto& Sq  = *F.find_scalar("S_combustion");
    for (std::size_t i = 0; i < nC; ++i) {
        T  [i] = lookup("T",            Z[i], chi[i]);
        rho[i] = lookup("rho",          Z[i], chi[i]);
        Sq [i] = lookup("S_combustion", Z[i], chi[i]);
    }
}

}  // namespace simall::combustion
