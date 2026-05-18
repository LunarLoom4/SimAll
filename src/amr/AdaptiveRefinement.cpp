// =============================================================================
// SimAll Beta - AMR Subsystem
// File   : src/amr/AdaptiveRefinement.cpp
// =============================================================================
#include "amr/AdaptiveRefinement.hpp"
#include "solver/Gradient.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace simall::amr {

void AdaptiveRefinement::compute_indicator(const meshing::Mesh& mesh,
                                           const solver::FieldRegistry& fields,
                                           AmrOptions opt) {
    const std::size_t nC = mesh.cells().size();
    ind_.assign(nC, 0.0);

    solver::LeastSquaresGradient G(mesh);
    auto& f = const_cast<solver::FieldRegistry&>(fields);

    auto h = [&](std::size_t c) { return std::cbrt(std::max(mesh.cells().volume[c], 1e-30)); };

    switch (opt.indicator) {
    case Indicator::GradientMagnitude: {
        const auto* phi = f.find_scalar(opt.fieldName);
        if (!phi) throw std::runtime_error("AMR: scalar field not found: " + opt.fieldName);
        solver::VectorField g;
        G.evaluate(*phi, g);
        for (std::size_t c = 0; c < nC; ++c) {
            const double gm = std::sqrt(g.x[c]*g.x[c] + g.y[c]*g.y[c] + g.z[c]*g.z[c]);
            ind_[c] = gm * h(c);
        }
        break;
    }
    case Indicator::SecondDerivative: {
        const auto* phi = f.find_scalar(opt.fieldName);
        if (!phi) throw std::runtime_error("AMR: scalar field not found: " + opt.fieldName);
        solver::VectorField g, gxx, gyy, gzz;
        G.evaluate(*phi, g);
        G.evaluate(g.x, gxx);
        G.evaluate(g.y, gyy);
        G.evaluate(g.z, gzz);
        for (std::size_t c = 0; c < nC; ++c) {
            const double lap = gxx.x[c] + gyy.y[c] + gzz.z[c];
            const double hh = h(c);
            ind_[c] = std::abs(lap) * hh * hh;
        }
        break;
    }
    case Indicator::Jump: {
        const auto* phi = f.find_scalar(opt.fieldName);
        if (!phi) throw std::runtime_error("AMR: scalar field not found: " + opt.fieldName);
        const auto& F = mesh.faces();
        const auto& C = mesh.cells();
        for (std::size_t c = 0; c < nC; ++c) {
            double mx = 0.0;
            const int s = C.faceOffsets[c];
            const int e = C.faceOffsets[c+1];
            for (int k = s; k < e; ++k) {
                const meshing::FaceId fid = C.faceIndices[k];
                const meshing::CellId nb  = (F.owner[fid] == c) ? F.neighbor[fid]
                                                                : F.owner[fid];
                if (nb == meshing::kBoundaryCell) continue;
                mx = std::max(mx, std::abs((*phi)[nb] - (*phi)[c]));
            }
            ind_[c] = mx;
        }
        break;
    }
    case Indicator::Vorticity: {
        const auto* U = f.find_vector("U");
        if (!U) throw std::runtime_error("AMR: vector field U not found");
        solver::VectorField gUx, gUy, gUz;
        G.evaluate(U->x, gUx);
        G.evaluate(U->y, gUy);
        G.evaluate(U->z, gUz);
        for (std::size_t c = 0; c < nC; ++c) {
            const double wx = gUz.y[c] - gUy.z[c];
            const double wy = gUx.z[c] - gUz.x[c];
            const double wz = gUy.x[c] - gUx.y[c];
            ind_[c] = std::sqrt(wx*wx + wy*wy + wz*wz) * h(c);
        }
        break;
    }}
    SIMALL_LOG_INFO("AMR", "indicator computed (", nC, " cells)");
}

void AdaptiveRefinement::mark(const std::vector<std::int32_t>& level,
                              std::vector<RefineFlag>& flags,
                              AmrOptions opt) const {
    const std::size_t nC = ind_.size();
    flags.assign(nC, RefineFlag::Keep);
    if (nC == 0) return;

    // Sort indices by indicator (ascending) and use partition thresholds.
    std::vector<std::size_t> order(nC);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(),
        [&](std::size_t a, std::size_t b){ return ind_[a] < ind_[b]; });

    const std::size_t nCoarse = static_cast<std::size_t>(opt.coarsenFraction * nC);
    const std::size_t nRefine = static_cast<std::size_t>(opt.refineFraction  * nC);

    for (std::size_t k = 0; k < nCoarse && k < nC; ++k) {
        const std::size_t c = order[k];
        if (level.empty() || level[c] > opt.minLevel)
            flags[c] = RefineFlag::Coarsen;
    }
    for (std::size_t k = 0; k < nRefine && k < nC; ++k) {
        const std::size_t c = order[nC - 1 - k];
        if (level.empty() || level[c] < opt.maxLevel)
            flags[c] = RefineFlag::Refine;
    }
}

}  // namespace simall::amr
