// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/SimpTopology.cpp
// =============================================================================
#include "optimization/SimpTopology.hpp"

#include "core/Logger.hpp"
#include "solver/CSRMatrix.hpp"
#include "solver/LinearSolvers.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace simall::optimization
{

bool SimpTopology::initialize(const meshing::Mesh& mesh,
                              solver::FieldRegistry& F,
                              const SimpConfig& cfg)
{
    cfg_ = cfg;
    mesh_ = &mesh;
    const std::size_t nC = mesh.cells().size();
    auto& rho = F.scalar("rho", nC);
    auto& rhoT = F.scalar("rhoTilde", nC);
    auto& kSimp = F.scalar("K_simp", nC);
    auto& dC = F.scalar("dC_drho", nC);
    auto& dCt = F.scalar("dC_drhoTilde", nC);
    std::fill(rho.begin(), rho.end(), cfg.volumeFraction);
    std::fill(rhoT.begin(), rhoT.end(), cfg.volumeFraction);
    std::fill(kSimp.begin(), kSimp.end(), cfg.K_min);
    std::fill(dC.begin(), dC.end(), 0.0);
    std::fill(dCt.begin(), dCt.end(), 0.0);

    if (cfg.filterRadius > 0.0) {
        solver::LinearSolverConfig lc;
        lc.tolerance = 1e-9;
        lc.maxIterations = 400;
        lc.restart = 30;
        lc.kind = solver::LinearSolverKind::CG;
        lc.preconditioner = solver::PreconditionerKind::Jacobi;
        filterSolver_ = solver::make_cg(lc);
    }
    return true;
}

void SimpTopology::apply_helmholtz_filter(const meshing::Mesh& mesh, solver::FieldRegistry& F)
{
    const auto& C = mesh.cells();
    const auto& Ff = mesh.faces();
    const std::size_t nC = C.size();
    auto* rho = F.find_scalar("rho");
    auto* rhoT = F.find_scalar("rhoTilde");
    if (!rho || !rhoT)
        return;

    if (cfg_.filterRadius <= 0.0) {
        for (std::size_t c = 0; c < nC; ++c)
            (*rhoT)[c] = (*rho)[c];
        return;
    }
    // Assemble FV discretisation of  −r² ∇²ρ̃ + ρ̃ = ρ  with zero-flux BCs.
    const double r2 = cfg_.filterRadius * cfg_.filterRadius;
    // Build CSR row layout: diagonal + one entry per interior face neighbour.
    std::vector<std::vector<std::pair<int, double>>> rows(nC);
    for (std::size_t c = 0; c < nC; ++c)
        rows[c].emplace_back(static_cast<int>(c), C.volume[c]);

    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const auto o = Ff.owner[f];
        const auto n = Ff.neighbor[f];
        if (n == meshing::kBoundaryCell)
            continue;
        const double a = std::sqrt(Ff.areaX[f] * Ff.areaX[f] + Ff.areaY[f] * Ff.areaY[f]
                                   + Ff.areaZ[f] * Ff.areaZ[f]);
        const double dx = C.centroidX[n] - C.centroidX[o];
        const double dy = C.centroidY[n] - C.centroidY[o];
        const double dz = C.centroidZ[n] - C.centroidZ[o];
        const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (d <= 0 || a <= 0)
            continue;
        const double coeff = r2 * a / d;
        // owner row
        bool found = false;
        for (auto& [j, v] : rows[o])
            if (j == static_cast<int>(o)) {
                v += coeff;
                found = true;
                break;
            }
        if (!found)
            rows[o].emplace_back(static_cast<int>(o), coeff);
        rows[o].emplace_back(static_cast<int>(n), -coeff);
        // neighbour row
        found = false;
        for (auto& [j, v] : rows[n])
            if (j == static_cast<int>(n)) {
                v += coeff;
                found = true;
                break;
            }
        if (!found)
            rows[n].emplace_back(static_cast<int>(n), coeff);
        rows[n].emplace_back(static_cast<int>(o), -coeff);
    }
    // Coalesce duplicates per row & sort by column.
    solver::CSRMatrix A;
    A.rowPtr.resize(nC + 1);
    int nnz = 0;
    for (std::size_t i = 0; i < nC; ++i) {
        auto& r = rows[i];
        std::sort(r.begin(), r.end(), [](auto& a, auto& b) { return a.first < b.first; });
        // Merge duplicates.
        std::vector<std::pair<int, double>> m;
        m.reserve(r.size());
        for (auto& e : r) {
            if (!m.empty() && m.back().first == e.first)
                m.back().second += e.second;
            else
                m.push_back(e);
        }
        r.swap(m);
        nnz += static_cast<int>(r.size());
    }
    A.colIdx.resize(nnz);
    A.values.resize(nnz);
    int k = 0;
    for (std::size_t i = 0; i < nC; ++i) {
        A.rowPtr[i] = k;
        for (auto& e : rows[i]) {
            A.colIdx[k] = e.first;
            A.values[k] = e.second;
            ++k;
        }
    }
    A.rowPtr[nC] = k;

    util::aligned_vector<double> b(nC), x(nC);
    for (std::size_t i = 0; i < nC; ++i) {
        b[i] = (*rho)[i] * C.volume[i];
        x[i] = (*rhoT)[i];
    }
    filterSolver_->solve(A, b, x);
    for (std::size_t i = 0; i < nC; ++i)
        (*rhoT)[i] = std::clamp(x[i], 0.0, 1.0);
}

double SimpTopology::step(const meshing::Mesh& mesh, solver::FieldRegistry& F)
{
    const auto& C = mesh.cells();
    const std::size_t nC = C.size();
    apply_helmholtz_filter(mesh, F);

    auto* rho = F.find_scalar("rho");
    auto* rhoT = F.find_scalar("rhoTilde");
    auto* kSimp = F.find_scalar("K_simp");
    auto* dC = F.find_scalar("dC_drho");
    auto* dCt = F.find_scalar("dC_drhoTilde");
    if (!rho || !rhoT || !kSimp || !dC || !dCt)
        return 0.0;

    const double p = cfg_.penalty;
    const double dK = cfg_.K_max - cfg_.K_min;
    for (std::size_t c = 0; c < nC; ++c)
        (*kSimp)[c] = cfg_.K_min + dK * std::pow((*rhoT)[c], p);

    // dC/drho_tilde = dC/dK * dK/drhoT  ; the caller stored ∂C/∂ρ in dC_drho
    // assuming it already accounted for chain rule via the property. Here we
    // additionally apply the Helmholtz-filter adjoint (identity-substitute
    // is sufficient when filterRadius=0 — for non-zero radius we approximate
    // by symmetry of the Helmholtz operator so dC/drho = (I − r²∇²)^-1 dC/drhoT.
    // To stay implementation-faithful we reuse the same solve as the forward
    // filter applied to the sensitivity field.
    if (cfg_.filterRadius > 0.0) {
        // Swap rho/rhoT temporarily: filter dC into dCt.
        for (std::size_t c = 0; c < nC; ++c)
            (*rho)[c] = (*dC)[c];
        apply_helmholtz_filter(mesh, F);
        for (std::size_t c = 0; c < nC; ++c) {
            (*dCt)[c] = (*rhoT)[c];
        }
        // Restore: rho/rhoT will be reinitialised by OC update below.
    } else {
        for (std::size_t c = 0; c < nC; ++c)
            (*dCt)[c] = (*dC)[c];
    }

    // OC update with bisection on Lagrange multiplier λ for volume constraint.
    double Vtot = 0.0;
    for (std::size_t c = 0; c < nC; ++c)
        Vtot += C.volume[c];
    const double Vtarget = cfg_.volumeFraction * Vtot;

    util::aligned_vector<double> rhoNew(nC);
    double l1 = 0.0, l2 = 1e9;
    for (int it = 0; it < 64 && (l2 - l1) > 1e-9 * (l1 + l2 + 1e-30); ++it) {
        const double lmid = 0.5 * (l1 + l2);
        double V = 0.0;
        for (std::size_t c = 0; c < nC; ++c) {
            const double s = -(*dCt)[c] / std::max(1e-12, lmid);
            const double bf = std::pow(std::max(0.0, s), cfg_.damping);
            const double old = (*rho)[c];
            double v = std::clamp(
                old * bf, std::max(cfg_.rhoMin, old - cfg_.move), std::min(1.0, old + cfg_.move));
            rhoNew[c] = v;
            V += v * C.volume[c];
        }
        if (V > Vtarget)
            l1 = lmid;
        else
            l2 = lmid;
    }

    double change = 0.0;
    for (std::size_t c = 0; c < nC; ++c) {
        change += std::abs(rhoNew[c] - (*rho)[c]);
        (*rho)[c] = rhoNew[c];
    }
    return change / static_cast<double>(nC);
}

} // namespace simall::optimization
