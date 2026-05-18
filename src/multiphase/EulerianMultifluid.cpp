// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/EulerianMultifluid.cpp
// =============================================================================
#include "multiphase/EulerianMultifluid.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase {

EulerianMultifluid::EulerianMultifluid(meshing::Mesh& mesh,
                                       solver::FieldRegistry& fields,
                                       solver::ILinearSolver& linear)
    : mesh_(mesh), F_(fields), lin_(linear) {}

void EulerianMultifluid::configure(std::vector<EulerianPhase> phases,
                                   const std::vector<solver::BoundarySpec>& bcs,
                                   double gx, double gy, double gz) {
    phases_ = std::move(phases);
    bcs_    = bcs;
    g_[0]=gx; g_[1]=gy; g_[2]=gz;

    const std::size_t nC = mesh_.cells().size();
    alphaEq_.clear();
    alphaEq_.reserve(phases_.size());

    // Register per-phase fields (alpha, U, T, h) and α-transport equation.
    double aDefault = 1.0 / std::max<std::size_t>(phases_.size(), 1);
    for (std::size_t k = 0; k < phases_.size(); ++k) {
        const auto& nm = phases_[k].name;
        auto& a  = F_.scalar("alpha_" + nm, nC);
        std::fill(a.begin(), a.end(), aDefault);
        F_.vector("U_"     + nm, nC);
        F_.scalar("T_"     + nm, nC);
        F_.scalar("h_"     + nm, nC);
        F_.scalar("rho_"   + nm, nC);
        F_.scalar("__em_src_alpha_" + nm, nC);   // mass-transfer source

        auto eq = std::make_unique<solver::ScalarTransport>(mesh_, F_, lin_);
        eq->set_field("alpha_" + nm);
        eq->set_density(phases_[k].rho);
        eq->set_source(std::string{"__em_src_alpha_" + nm});
        eq->set_urf(0.7);
        for (const auto& b : bcs_) {
            switch (b.type) {
                case solver::BCType::NoSlipWall:
                case solver::BCType::Wall:
                    eq->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
                    break;
                case solver::BCType::VelocityInlet:
                case solver::BCType::PressureInlet:
                case solver::BCType::MassFlowInlet:
                    eq->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, aDefault, 0.0});
                    break;
                default:
                    eq->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            }
        }
        alphaEq_.push_back(std::move(eq));
    }

    // Build drag-pair table for ordered (i<j) phase pairs.
    drag_.clear();
    for (std::size_t i = 0; i < phases_.size(); ++i)
        for (std::size_t j = i+1; j < phases_.size(); ++j)
            drag_.push_back({i, j, 0.0});

    // Mixture-drag coupling fields (one explicit momentum source per phase).
    for (std::size_t k = 0; k < phases_.size(); ++k) {
        F_.vector("S_mom_" + phases_[k].name, nC);
        F_.scalar("S_eng_" + phases_[k].name, nC);
    }
    F_.scalar("rho_mix", nC);

    SIMALL_LOG_INFO("Multiphase",
        "EulerianMultifluid configured: ", phases_.size(),
        " phases, ", drag_.size(), " interphase pairs (", nC, " cells)");
}

void EulerianMultifluid::update_mixture_density() {
    const std::size_t nC = mesh_.cells().size();
    auto& rhoM = *F_.find_scalar("rho_mix");
    std::fill(rhoM.begin(), rhoM.end(), 0.0);
    for (std::size_t k = 0; k < phases_.size(); ++k) {
        const auto& a   = *F_.find_scalar("alpha_" + phases_[k].name);
        auto&       rhk = *F_.find_scalar("rho_"   + phases_[k].name);
        for (std::size_t c = 0; c < nC; ++c) {
            rhk[c]   = phases_[k].rho;
            rhoM[c] += a[c] * phases_[k].rho;
        }
    }
}

void EulerianMultifluid::update_drag_coefficients() {
    const std::size_t nC = mesh_.cells().size();
    const auto& rhoM = *F_.find_scalar("rho_mix");
    for (auto& dp : drag_) {
        const std::size_t i = dp.i, j = dp.j;
        const auto& Ui = *F_.find_vector("U_" + phases_[i].name);
        const auto& Uj = *F_.find_vector("U_" + phases_[j].name);
        // dispersed phase = the one marked dispersed (or smaller diameter)
        const std::size_t disp = phases_[i].dispersed ? i : j;
        const std::size_t cont = (disp == i) ? j : i;
        (void)cont;
        const double dp_dia = std::max(phases_[disp].diameter, 1e-12);
        const double mu_c   = std::max(phases_[cont].mu, 1e-12);

        double Kavg = 0.0;
        for (std::size_t c = 0; c < nC; ++c) {
            const double dux = Ui.x[c]-Uj.x[c];
            const double duy = Ui.y[c]-Uj.y[c];
            const double duz = Ui.z[c]-Uj.z[c];
            const double Ur  = std::sqrt(dux*dux + duy*duy + duz*duz);
            const double Re  = rhoM[c] * Ur * dp_dia / mu_c;
            const double Cd  = (Re < 1000.0)
                ? (24.0 / std::max(Re, 1e-6)) * (1.0 + 0.15 * std::pow(std::max(Re, 1e-6), 0.687))
                : 0.44;
            // K_ij = (3/4) ρ_m C_d |U_r| / d_p
            Kavg += 0.75 * rhoM[c] * Cd * Ur / dp_dia;
        }
        dp.K = (nC > 0) ? Kavg / static_cast<double>(nC) : 0.0;
    }
}

void EulerianMultifluid::apply_interphase_momentum() {
    const std::size_t nC = mesh_.cells().size();
    // Reset per-phase momentum sources.
    for (std::size_t k = 0; k < phases_.size(); ++k) {
        auto& S = *F_.find_vector("S_mom_" + phases_[k].name);
        std::fill(S.x.begin(), S.x.end(), 0.0);
        std::fill(S.y.begin(), S.y.end(), 0.0);
        std::fill(S.z.begin(), S.z.end(), 0.0);
        // Gravity contribution α_k ρ_k g.
        const auto& a = *F_.find_scalar("alpha_" + phases_[k].name);
        const double rho = phases_[k].rho;
        for (std::size_t c = 0; c < nC; ++c) {
            S.x[c] += a[c] * rho * g_[0];
            S.y[c] += a[c] * rho * g_[1];
            S.z[c] += a[c] * rho * g_[2];
        }
    }
    // Drag: K_ij (U_j − U_i) added to phase i, opposite to phase j.
    for (const auto& dp : drag_) {
        const auto& Ui = *F_.find_vector("U_" + phases_[dp.i].name);
        const auto& Uj = *F_.find_vector("U_" + phases_[dp.j].name);
        auto&       Si = *F_.find_vector("S_mom_" + phases_[dp.i].name);
        auto&       Sj = *F_.find_vector("S_mom_" + phases_[dp.j].name);
        for (std::size_t c = 0; c < nC; ++c) {
            const double Mx = dp.K * (Uj.x[c] - Ui.x[c]);
            const double My = dp.K * (Uj.y[c] - Ui.y[c]);
            const double Mz = dp.K * (Uj.z[c] - Ui.z[c]);
            Si.x[c] += Mx; Si.y[c] += My; Si.z[c] += Mz;
            Sj.x[c] -= Mx; Sj.y[c] -= My; Sj.z[c] -= Mz;
        }
    }
}

void EulerianMultifluid::apply_interphase_energy() {
    const std::size_t nC = mesh_.cells().size();
    // Interphase heat transfer Q_ij = h_ij A_ij (T_j − T_i)
    // h_ij from Ranz-Marshall: Nu = 2 + 0.6 Re^{1/2} Pr^{1/3}, h = Nu k_c / d_p.
    for (std::size_t k = 0; k < phases_.size(); ++k) {
        auto& S = *F_.find_scalar("S_eng_" + phases_[k].name);
        std::fill(S.begin(), S.end(), 0.0);
    }
    const auto& rhoM = *F_.find_scalar("rho_mix");
    for (const auto& dp : drag_) {
        const std::size_t i = dp.i, j = dp.j;
        const auto& Ti = *F_.find_scalar("T_" + phases_[i].name);
        const auto& Tj = *F_.find_scalar("T_" + phases_[j].name);
        const auto& ai = *F_.find_scalar("alpha_" + phases_[i].name);
        const auto& Ui = *F_.find_vector("U_" + phases_[i].name);
        const auto& Uj = *F_.find_vector("U_" + phases_[j].name);
        auto&       Si = *F_.find_scalar("S_eng_" + phases_[i].name);
        auto&       Sj = *F_.find_scalar("S_eng_" + phases_[j].name);

        const std::size_t disp = phases_[i].dispersed ? i : j;
        const std::size_t cont = (disp == i) ? j : i;
        const double d_p = std::max(phases_[disp].diameter, 1e-12);
        const double k_c = std::max(phases_[cont].k_thermal, 1e-12);
        const double mu_c= std::max(phases_[cont].mu, 1e-12);
        const double Pr  = mu_c * phases_[cont].cp / k_c;

        // Interfacial area density Ai = 6 α_disp / d_p (spherical dispersed phase)
        for (std::size_t c = 0; c < nC; ++c) {
            const double dux = Ui.x[c]-Uj.x[c];
            const double duy = Ui.y[c]-Uj.y[c];
            const double duz = Ui.z[c]-Uj.z[c];
            const double Ur  = std::sqrt(dux*dux + duy*duy + duz*duz);
            const double Re  = rhoM[c] * Ur * d_p / mu_c;
            const double Nu  = 2.0 + 0.6 * std::sqrt(std::max(Re,0.0)) * std::cbrt(Pr);
            const double h   = Nu * k_c / d_p;
            const double a_disp = (disp == i) ? ai[c]
                                              : (*F_.find_scalar("alpha_" + phases_[disp].name))[c];
            const double Ai  = 6.0 * a_disp / d_p;
            const double Q   = h * Ai * (Tj[c] - Ti[c]);
            Si[c] += Q;
            Sj[c] -= Q;
        }
    }
}

void EulerianMultifluid::step(double /*dt*/) {
    update_mixture_density();
    update_drag_coefficients();
    apply_interphase_momentum();
    apply_interphase_energy();
    // α-equations advanced one outer iteration each.  The ScalarTransport
    // diffusivity is set to 0 (pure advection) — interface compression is
    // delegated to the VOF / IsoAdvector pathway when relevant.
    for (auto& eq : alphaEq_) {
        eq->set_diffusivity(0.0);
        eq->solve_iteration();
    }
    // Volume-fraction sum normalisation (clip + renormalise so Σα_k = 1).
    const std::size_t nC = mesh_.cells().size();
    for (std::size_t c = 0; c < nC; ++c) {
        double s = 0.0;
        for (auto& ph : phases_) {
            auto& a = *F_.find_scalar("alpha_" + ph.name);
            a[c] = std::clamp(a[c], 0.0, 1.0);
            s += a[c];
        }
        const double inv = (s > 1e-30) ? 1.0/s : 0.0;
        for (auto& ph : phases_) {
            auto& a = *F_.find_scalar("alpha_" + ph.name);
            a[c] *= inv;
        }
    }
}

}  // namespace simall::multiphase
