// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/LRR_Reynolds_Stress.cpp
// =============================================================================
#include "turbulence/LRR_Reynolds_Stress.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence {

namespace {
constexpr double Cmu   = 0.09;
constexpr double sigK  = 1.0;
constexpr double sigE  = 1.3;
constexpr double C1    = 1.8;
constexpr double C2    = 0.6;
constexpr double C1eps = 1.44;
constexpr double C2eps = 1.92;
constexpr double kFloor = 1.0e-12;
constexpr double eFloor = 1.0e-12;

// Mapping of independent R_ij components.
// idx 0:R11(xx), 1:R22(yy), 2:R33(zz), 3:R12(xy), 4:R13(xz), 5:R23(yz)
constexpr int IJ[6][2] = {{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}};
const char* const RNAMES[6] = {"R11","R22","R33","R12","R13","R23"};
const char* const SRC_R[6]  = {"__rsm_src_R11","__rsm_src_R22","__rsm_src_R33",
                               "__rsm_src_R12","__rsm_src_R13","__rsm_src_R23"};
constexpr const char* SRC_E = "__rsm_src_eps";
}  // namespace

void LRR_Reynolds_Stress_Full::initialize(meshing::Mesh& m, solver::FieldRegistry& f) {
    mesh_ = &m;
    const std::size_t nC = m.cells().size();

    solver::LinearSolverConfig c;
    c.kind = solver::LinearSolverKind::GMRES;
    c.preconditioner = solver::PreconditionerKind::ILU;
    c.tolerance = 1e-6; c.maxIterations = 100; c.restart = 20;
    lin_ = solver::make_linear_solver(c);

    for (int i = 0; i < 6; ++i) {
        f.scalar(RNAMES[i], nC);
        f.scalar(SRC_R[i],  nC);
    }
    f.scalar("epsilon", nC);
    f.scalar(SRC_E,     nC);
    f.scalar("mut",     nC);

    // Initial isotropic state: R11 = R22 = R33 = ⅔ k0, off-diagonal = 0.
    const double k0 = 1.0e-4, e0 = 1.0e-4;
    for (int i = 0; i < 3; ++i) {
        auto& R = *f.find_scalar(RNAMES[i]);
        std::fill(R.begin(), R.end(), (2.0/3.0)*k0);
    }
    for (int i = 3; i < 6; ++i) {
        auto& R = *f.find_scalar(RNAMES[i]);
        std::fill(R.begin(), R.end(), 0.0);
    }
    auto& E = *f.find_scalar("epsilon");
    std::fill(E.begin(), E.end(), e0);

    mut_.assign(nC, 0.0);

    for (int i = 0; i < 6; ++i) {
        Req_[i] = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
        Req_[i]->set_field(RNAMES[i]);
        Req_[i]->set_density(rho_);
        Req_[i]->set_source(std::string{SRC_R[i]});
        Req_[i]->set_urf(0.5);
        for (const auto& b : bcs_) {
            switch (b.type) {
                case solver::BCType::NoSlipWall:
                case solver::BCType::Wall:
                    Req_[i]->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
                    break;
                case solver::BCType::VelocityInlet:
                case solver::BCType::PressureInlet:
                case solver::BCType::MassFlowInlet:
                    if (i < 3)  // diagonal
                        Req_[i]->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-3*(2.0/3.0), 0.0});
                    else
                        Req_[i]->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 0.0, 0.0});
                    break;
                default:
                    Req_[i]->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
            }
        }
    }
    eEq_ = std::make_unique<solver::ScalarTransport>(m, f, *lin_);
    eEq_->set_field("epsilon");
    eEq_->set_density(rho_);
    eEq_->set_source(std::string{SRC_E});
    eEq_->set_urf(0.5);
    for (const auto& b : bcs_) {
        switch (b.type) {
            case solver::BCType::NoSlipWall:
            case solver::BCType::Wall:
                eEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
                break;
            case solver::BCType::VelocityInlet:
            case solver::BCType::PressureInlet:
            case solver::BCType::MassFlowInlet:
                eEq_->add_bc({b.zone, solver::ScalarBC::Kind::Dirichlet, 1.0e-3, 0.0});
                break;
            default:
                eEq_->add_bc({b.zone, solver::ScalarBC::Kind::Neumann, 0.0, 0.0});
        }
    }
    SIMALL_LOG_INFO("Turbulence", "LRR Reynolds-Stress (7 eqs) initialised (cells=", nC, ")");
}

void LRR_Reynolds_Stress_Full::solve(double /*dt*/, solver::FieldRegistry& f) {
    if (!mesh_) return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U  = *f.find_vector("U");
    auto& mut = *f.find_scalar("mut");
    auto& E   = *f.find_scalar("epsilon");
    auto& Se  = *f.find_scalar(SRC_E);

    // Pointers to all R components and source fields.
    std::array<util::aligned_vector<double>*, 6> R{};
    std::array<util::aligned_vector<double>*, 6> Sr{};
    for (int i = 0; i < 6; ++i) {
        R[i]  = f.find_scalar(RNAMES[i]);
        Sr[i] = f.find_scalar(SRC_R[i]);
    }

    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx); G.evaluate(U.y, gUy); G.evaluate(U.z, gUz);

    for (std::size_t c = 0; c < nC; ++c) {
        // Build full R tensor (symmetric).
        double Rt[3][3];
        Rt[0][0] = (*R[0])[c]; Rt[1][1] = (*R[1])[c]; Rt[2][2] = (*R[2])[c];
        Rt[0][1] = Rt[1][0] = (*R[3])[c];
        Rt[0][2] = Rt[2][0] = (*R[4])[c];
        Rt[1][2] = Rt[2][1] = (*R[5])[c];

        const double k_c = 0.5 * (Rt[0][0] + Rt[1][1] + Rt[2][2]);
        const double kf  = std::max(k_c, kFloor);
        const double ef  = std::max(E[c], eFloor);

        // ∂U_i/∂x_j matrix.
        double dU[3][3] = {
            {gUx.x[c], gUx.y[c], gUx.z[c]},
            {gUy.x[c], gUy.y[c], gUy.z[c]},
            {gUz.x[c], gUz.y[c], gUz.z[c]}
        };

        // P_ij = − R_ik dU_j/dx_k − R_jk dU_i/dx_k
        double Pij[3][3];
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0.0;
                for (int kk = 0; kk < 3; ++kk)
                    s -= Rt[i][kk]*dU[j][kk] + Rt[j][kk]*dU[i][kk];
                Pij[i][j] = s;
            }
        const double Pk = 0.5 * (Pij[0][0] + Pij[1][1] + Pij[2][2]);

        // Φ_ij = − C1 ρ ε/k (R_ij − ⅔ δ k)  − C2 (P_ij − ⅔ δ P)
        const double e_over_k = ef / kf;

        // Source for each R_ij = production + pressure-strain − ε_ij(=⅔δ_ijε).
        for (int idx = 0; idx < 6; ++idx) {
            const int i = IJ[idx][0], j = IJ[idx][1];
            const double dij = (i == j) ? 1.0 : 0.0;
            const double Phi1 = - C1 * rho_ * e_over_k * (Rt[i][j] - (2.0/3.0)*dij*kf);
            const double Phi2 = - C2 * (Pij[i][j] - (2.0/3.0)*dij*Pk);
            const double dis  = (2.0/3.0) * dij * rho_ * ef;
            (*Sr[idx])[c] = rho_*Pij[i][j] + Phi1 + Phi2 - dis;
        }

        // ε source: (ε/k)(C1ε P_k − C2ε ρ ε)
        Se[c] = e_over_k * (C1eps * rho_ * Pk - C2eps * rho_ * ef);

        // μ_t (for momentum) from trace.
        const double mut_c = rho_ * Cmu * kf*kf / ef;
        mut[c]  = mut_c;
        mut_[c] = mut_c;
    }

    double mutAvg = 0.0;
    for (std::size_t c = 0; c < nC; ++c) mutAvg += mut[c];
    mutAvg /= std::max<std::size_t>(nC, 1);

    // R-equations use generalized gradient diffusion D = (ν + ν_t/σ_k) ∇R.
    for (int idx = 0; idx < 6; ++idx) {
        Req_[idx]->set_diffusivity(mu_ + mutAvg / sigK);
        Req_[idx]->solve_iteration();
    }
    eEq_->set_diffusivity(mu_ + mutAvg / sigE);
    eEq_->solve_iteration();

    // Realisability clipping: diagonals ≥ 0; |R_ij| ≤ √(R_ii R_jj).
    for (std::size_t c = 0; c < nC; ++c) {
        for (int i = 0; i < 3; ++i)
            (*R[i])[c] = std::max((*R[i])[c], 0.0);
        // Off-diagonals: idx 3=xy(0,1), 4=xz(0,2), 5=yz(1,2)
        for (int idx = 3; idx < 6; ++idx) {
            const int i = IJ[idx][0], j = IJ[idx][1];
            const double lim = std::sqrt(std::max((*R[i])[c]*(*R[j])[c], 0.0));
            (*R[idx])[c] = std::clamp((*R[idx])[c], -lim, lim);
        }
        E[c] = std::max(E[c], eFloor);
    }
}

}  // namespace simall::turbulence
