// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/FgmFlamelet.cpp
// =============================================================================
#include "combustion/FgmFlamelet.hpp"
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

void FgmFlamelet::build_analytic(std::size_t NZ, std::size_t Nc,
                                 double Z_st, double T_un,
                                 double T_ad, double rho_un, double rho_ad) {
    FgmTable t;
    t.axis.Z.resize(NZ);
    t.axis.c.resize(Nc);
    for (std::size_t i = 0; i < NZ; ++i) t.axis.Z[i] = double(i) / double(NZ - 1);
    for (std::size_t i = 0; i < Nc; ++i) t.axis.c[i] = double(i) / double(Nc - 1);

    std::vector<double> T(NZ * Nc), rho(NZ * Nc), wdot(NZ * Nc), Q(NZ * Nc);
    const double Hc = 50.0e6;   // J/kg fuel LHV (methane)
    for (std::size_t iZ = 0; iZ < NZ; ++iZ) {
        const double Z = t.axis.Z[iZ];
        // Burke-Schumann linear: T_eq(Z) ramp triangular about Z_st
        const double T_eq = (Z <= Z_st)
            ? T_un + (T_ad - T_un) * Z / std::max(Z_st, 1e-12)
            : T_un + (T_ad - T_un) * std::max(0.0, 1.0 - (Z - Z_st) / std::max(1.0 - Z_st, 1e-12));
        for (std::size_t ic = 0; ic < Nc; ++ic) {
            const double c = t.axis.c[ic];
            const double Tcell  = T_un + c * (T_eq - T_un);
            const double rhocell = rho_un + c * (rho_ad - rho_un);
            // Damköhler-style reaction rate peaks at c≈0.7, scaled by Z·(1-Z)
            const double w = std::max(0.0, std::sin(M_PI * c)) * 4.0 * Z * (1.0 - Z) * 100.0;
            T  [iZ * Nc + ic] = Tcell;
            rho[iZ * Nc + ic] = std::max(rhocell, 1e-3);
            wdot[iZ * Nc + ic] = w;
            Q  [iZ * Nc + ic] = w * Hc;
        }
    }
    t.var["T"]            = std::move(T);
    t.var["rho"]          = std::move(rho);
    t.var["S_progress"]   = std::move(wdot);
    t.var["S_combustion"] = std::move(Q);
    table_ = std::move(t);
    SIMALL_LOG_INFO("Combustion",
        "FGM analytic table built: ", NZ, " × ", Nc,
        " (vars=", table_.var.size(), ")");
}

double FgmFlamelet::lookup(const std::string& var, double Z, double c) const {
    auto it = table_.var.find(var);
    if (it == table_.var.end() || table_.axis.Z.empty() || table_.axis.c.empty())
        return 0.0;
    const std::size_t Nc = table_.axis.c.size();
    const std::size_t i0 = lower_index(table_.axis.Z, Z);
    const std::size_t j0 = lower_index(table_.axis.c, c);
    const double      z0 = table_.axis.Z[i0],   z1 = table_.axis.Z[i0 + 1];
    const double      c0 = table_.axis.c[j0],   c1 = table_.axis.c[j0 + 1];
    const double      tz = std::clamp((Z - z0) / std::max(z1 - z0, 1e-30), 0.0, 1.0);
    const double      tc = std::clamp((c - c0) / std::max(c1 - c0, 1e-30), 0.0, 1.0);
    const auto&       T  = it->second;
    const double      v00 = T[i0       * Nc + j0    ];
    const double      v10 = T[(i0 + 1) * Nc + j0    ];
    const double      v01 = T[i0       * Nc + j0 + 1];
    const double      v11 = T[(i0 + 1) * Nc + j0 + 1];
    return (1 - tz) * (1 - tc) * v00 + tz * (1 - tc) * v10
         + (1 - tz) *      tc  * v01 + tz *      tc  * v11;
}

void FgmFlamelet::initialize(const meshing::Mesh& mesh,
                             solver::FieldRegistry& F) {
    mesh_ = &mesh;
    const std::size_t nC = mesh.cells().size();
    F.scalar("Z", nC);
    F.scalar("c", nC);
    F.scalar("T", nC);
    F.scalar("rho", nC);
    F.scalar("S_progress",   nC);
    F.scalar("S_combustion", nC);
    if (table_.var.empty()) build_analytic();
    SIMALL_LOG_INFO("Combustion", "FGM initialised on ", nC, " cells");
}

void FgmFlamelet::apply(solver::FieldRegistry& F) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& Z = *F.find_scalar("Z");
    const auto& c = *F.find_scalar("c");
    auto& T   = *F.find_scalar("T");
    auto& rho = *F.find_scalar("rho");
    auto& Sw  = *F.find_scalar("S_progress");
    auto& Sq  = *F.find_scalar("S_combustion");
    for (std::size_t i = 0; i < nC; ++i) {
        T  [i] = lookup("T",            Z[i], c[i]);
        rho[i] = lookup("rho",          Z[i], c[i]);
        Sw [i] = lookup("S_progress",   Z[i], c[i]);
        Sq [i] = lookup("S_combustion", Z[i], c[i]);
    }
}

}  // namespace simall::combustion
