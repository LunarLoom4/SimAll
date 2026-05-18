// =============================================================================
// SimAll Beta - Combustion Subsystem
// File   : src/combustion/SootMomMethod.cpp
// =============================================================================
#include "combustion/SootMomMethod.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::combustion {

namespace {
constexpr double Ru   = 8.314462618;
constexpr double Avo  = 6.02214076e23;
constexpr double PI   = 3.14159265358979323846;
constexpr double kB   = 1.380649e-23;
constexpr double C_mw = 12.011e-3;
constexpr double v_C  = 1.66e-29;          // m³  primary carbon atom volume
}

void SootMomMethod::initialize(const meshing::Mesh& mesh, MomicProps props) {
    mesh_ = &mesh;
    p_    = props;
    const std::size_t nC = mesh.cells().size();
    SIMALL_LOG_INFO("Combustion",
        "SootMomMethod (MOMIC K=", MOMIC_K, ") initialised on ", nC, " cells");
}

double SootMomMethod::fractional_moment(const std::array<double, MOMIC_K>& M, double r) {
    // Log-linear interpolation between adjacent integer moments:
    //   log10 M(r) = log10 M_i + (r - i) (log10 M_{i+1} - log10 M_i)
    // Bracket r between [floor(r), floor(r)+1] and clamp to [0, K-1].
    if (r <= 0.0) return std::max(M[0], 1e-300);
    if (r >= static_cast<double>(MOMIC_K - 1)) return std::max(M[MOMIC_K - 1], 1e-300);
    const int i = static_cast<int>(std::floor(r));
    const double t = r - i;
    const double l0 = std::log10(std::max(M[i],     1e-300));
    const double l1 = std::log10(std::max(M[i + 1], 1e-300));
    return std::pow(10.0, l0 + t * (l1 - l0));
}

double SootMomMethod::apply(solver::FieldRegistry& F) {
    if (!mesh_) return 0.0;
    const std::size_t nC = mesh_->cells().size();
    const auto* T   = F.find_scalar("T");
    if (!T) return 0.0;
    const auto* rhoF = F.find_scalar("rho_mix");
    const auto* Yac  = F.find_scalar("Y_C2H2");
    const auto* Yo2  = F.find_scalar("Y_O2");
    const auto* Yoh  = F.find_scalar("Y_OH");

    // Field accessors for the K moments and their sources.
    std::array<solver::ScalarField*, MOMIC_K> Mfield{}, Sfield{};
    for (int r = 0; r < MOMIC_K; ++r) {
        Mfield[r] = &F.scalar("M" + std::to_string(r) + "_soot", nC);
        Sfield[r] = &F.scalar("S_M" + std::to_string(r) + "_soot", nC);
    }

    const auto& Cc = mesh_->cells();
    double totalM0 = 0.0;
    for (std::size_t i = 0; i < nC; ++i) {
        const double Tc  = std::max((*T)[i], 300.0);
        const double rho = rhoF ? std::max((*rhoF)[i], 1e-6) : 1.0;
        const double yA  = Yac  ? std::max((*Yac)[i], 0.0) : 0.0;
        const double yO2 = Yo2  ? std::max((*Yo2)[i], 0.0) : 0.0;
        const double yOH = Yoh  ? std::max((*Yoh)[i], 0.0) : 0.0;

        const double cA  = rho * yA  / p_.Mw_C2H2;       // mol/m³
        const double cO2 = rho * yO2 / 31.998e-3;
        const double cOH = rho * yOH / 17.008e-3;

        std::array<double, MOMIC_K> M{};
        for (int r = 0; r < MOMIC_K; ++r) M[r] = std::max((*Mfield[r])[i], 0.0);

        // 1) Inception (acetylene → C16 PAH primary): adds R_inc particles of
        //    fixed volume v_0 = 16·v_C; contributes to all moments via v^r.
        const double R_inc = p_.C_inc * cA * std::exp(-p_.T_inc / Tc);  // 1/(m³·s)
        const double v0    = 16.0 * v_C;

        // 2) Surface growth (HACA): mass-flux per surface area · A_s
        const double R_sg  = p_.C_sg * std::exp(-p_.T_sg / Tc) * cA;     // m/s
        // Effective particle volume from M_1 / M_0.
        const double v_avg = (M[0] > 1e-30) ? M[1] / M[0] : v0;
        const double d_avg = std::cbrt(6.0 * v_avg / PI);
        // Total surface area density ≈ π · M_2 (after Frenklach 2002 closure
        // since d ~ v^{1/3} ⇒ s ~ v^{2/3} ⇒ Σ N s = (π) · M_{2/3} scaled).
        const double M_23  = fractional_moment(M, 2.0/3.0);
        const double A_s   = PI * std::cbrt(36.0 * PI) * M_23;

        // 3) Oxidation (O2 + OH) — both reduce volume.
        const double R_ox_O2 = p_.C_ox_O2 * cO2 * std::exp(-p_.T_ox_O2 / Tc); // mol/(m²·s) per A_s
        const double R_ox_OH = p_.C_ox_OH * cOH * std::sqrt(Tc);
        const double dvol_dt = (R_sg - (R_ox_O2 + R_ox_OH) * C_mw / p_.rho_soot) * A_s;

        // 4) Coagulation (free-molecular Smoluchowski): only impacts M_0 (sink)
        //    and higher moments through M_{1/6} and M_{2/3} weights.
        const double K_fm = std::sqrt(6.0 * kB * Tc / p_.rho_soot)
                          * std::pow(3.0 / (4.0 * PI), 1.0/6.0);
        const double M_16 = fractional_moment(M, 1.0/6.0);
        const double R_coag_M0 = -0.5 * K_fm * M_16 * M_23;             // 1/(m³·s)

        // Source assembly per moment.
        // Inception: dM_r/dt += R_inc · v0^r
        // Surface growth+oxidation acts on M_r by r·dv_dt scaling at v_avg.
        // Coagulation acts on M_0 primarily.
        Sfield[0]->at(i) += R_inc + R_coag_M0;
        for (int r = 1; r < MOMIC_K; ++r) {
            Sfield[r]->at(i) += R_inc * std::pow(v0, r)
                              + r * std::pow(std::max(v_avg, 1e-30), r - 1) * dvol_dt * M[0];
        }
        totalM0 += (R_inc + R_coag_M0) * Cc.volume[i];
    }
    return totalM0;
}

}  // namespace simall::combustion
