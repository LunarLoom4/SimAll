// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/P1Radiation.cpp
// =============================================================================
#include "radiation/P1Radiation.hpp"
#include "solver/CSRMatrix.hpp"
#include "solver/LinearSolvers.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_set>

namespace { constexpr double kPi = 3.14159265358979323846; }

namespace simall::radiation {

namespace {
constexpr double kSigmaSB = 5.670374419e-8;   // Stefan-Boltzmann [W/(m²·K⁴)]
}

bool P1Radiation::initialize(const meshing::Mesh& mesh,
                             solver::FieldRegistry& F,
                             const P1Props& props)
{
    mesh_ = &mesh; F_ = &F; p_ = props;
    const std::size_t nC = mesh.cells().size();
    auto& G    = F.scalar("G",     nC);
    auto& Srad = F.scalar("S_rad", nC);
    auto& T    = F.scalar("T",     nC);
    (void)T;
    std::fill(G.begin(),    G.end(),    0.0);
    std::fill(Srad.begin(), Srad.end(), 0.0);

    solver::LinearSolverConfig cfg;
    cfg.kind = solver::LinearSolverKind::CG;
    cfg.preconditioner = solver::PreconditionerKind::ILU;
    cfg.tolerance = 1e-8; cfg.maxIterations = 500;
    lin_ = solver::make_cg(cfg);
    return true;
}

double P1Radiation::solve() {
    const auto& M  = *mesh_;
    const auto& C  = M.cells();
    const auto& Ff = M.faces();
    const std::size_t nC = C.size();

    const auto* T = F_->find_scalar("T");
    auto*       G = F_->find_scalar("G");
    auto*       S = F_->find_scalar("S_rad");
    if (!T || !G || !S) return -1.0;

    const double a  = p_.absorption;
    const double ss = p_.scattering;
    const double A  = p_.asymmetry;
    const double Gam = 1.0 / (3.0 * (a + ss - A * ss + 1e-30));

    std::unordered_set<meshing::ZoneId> wallZones;
    for (const auto& w : walls_) wallZones.insert(w.zone);

    // Assemble CSR for  −Γ ∇²G + a G = 4 π a I_b   in finite-volume form.
    std::vector<std::vector<std::pair<int,double>>> rows(nC);
    util::aligned_vector<double> b(nC, 0.0);

    for (std::size_t c = 0; c < nC; ++c)
        rows[c].emplace_back(static_cast<int>(c), a * C.volume[c]);

    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const auto o = Ff.owner[f];
        const auto n = Ff.neighbor[f];
        const double areaMag = std::sqrt(Ff.areaX[f]*Ff.areaX[f]
                                       + Ff.areaY[f]*Ff.areaY[f]
                                       + Ff.areaZ[f]*Ff.areaZ[f]);
        if (n != meshing::kBoundaryCell) {
            const double dx = C.centroidX[n]-C.centroidX[o];
            const double dy = C.centroidY[n]-C.centroidY[o];
            const double dz = C.centroidZ[n]-C.centroidZ[o];
            const double d  = std::sqrt(dx*dx+dy*dy+dz*dz);
            if (d <= 0 || areaMag <= 0) continue;
            const double coef = Gam * areaMag / d;
            rows[o].emplace_back(static_cast<int>(o),  coef);
            rows[o].emplace_back(static_cast<int>(n), -coef);
            rows[n].emplace_back(static_cast<int>(n),  coef);
            rows[n].emplace_back(static_cast<int>(o), -coef);
        } else {
            const auto z = Ff.boundaryZone[f];
            if (!wallZones.count(z)) continue;
            // Find matching wall spec.
            const P1WallBC* wp = nullptr;
            for (const auto& w : walls_) if (w.zone == z) { wp = &w; break; }
            if (!wp) continue;
            const double eps = std::clamp(wp->emissivity, 1e-3, 1.0);
            const double coef = areaMag * eps / (2.0 * (2.0 - eps));
            const double Eb = 4.0 * kSigmaSB * std::pow(wp->temperature, 4);
            rows[o].emplace_back(static_cast<int>(o), coef);
            b[o] += coef * Eb;
        }
    }
    // Source term from emission.
    for (std::size_t c = 0; c < nC; ++c) {
        const double Ib = kSigmaSB * std::pow((*T)[c], 4);
        b[c] += 4.0 * kPi * a * Ib * C.volume[c];
    }
    // Coalesce.
    solver::CSRMatrix Aop;
    Aop.rowPtr.resize(nC + 1);
    int nnz = 0;
    for (std::size_t i = 0; i < nC; ++i) {
        auto& r = rows[i];
        std::sort(r.begin(), r.end(), [](auto&a,auto&b){return a.first<b.first;});
        std::vector<std::pair<int,double>> m; m.reserve(r.size());
        for (auto& e : r) {
            if (!m.empty() && m.back().first == e.first) m.back().second += e.second;
            else m.push_back(e);
        }
        r.swap(m); nnz += static_cast<int>(r.size());
    }
    Aop.colIdx.resize(nnz); Aop.values.resize(nnz);
    int k = 0;
    for (std::size_t i = 0; i < nC; ++i) {
        Aop.rowPtr[i] = k;
        for (auto& e : rows[i]) { Aop.colIdx[k] = e.first; Aop.values[k] = e.second; ++k; }
    }
    Aop.rowPtr[nC] = k;

    util::aligned_vector<double> x(nC);
    for (std::size_t i = 0; i < nC; ++i) x[i] = (*G)[i];
    const int iters = lin_->solve(Aop, b, x);
    for (std::size_t i = 0; i < nC; ++i) (*G)[i] = std::max(0.0, x[i]);

    // Radiative source.
    for (std::size_t c = 0; c < nC; ++c) {
        const double Ib4 = 4.0 * kSigmaSB * std::pow((*T)[c], 4);
        (*S)[c] = a * ((*G)[c] - Ib4);
    }
    (void)iters;
    return lin_->last_residual();
}

}  // namespace simall::radiation
