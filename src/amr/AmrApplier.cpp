// =============================================================================
// SimAll Beta - AMR Subsystem
// File   : src/amr/AmrApplier.cpp
// =============================================================================
#include "amr/AmrApplier.hpp"

#include <algorithm>
#include <numeric>

namespace simall::amr {

AmrCycleResult AmrApplier::apply_cycle(const std::vector<double>& ind,
                                         std::vector<std::uint8_t>& level,
                                         const std::vector<std::size_t>& rowPtr,
                                         const std::vector<std::size_t>& idx,
                                         std::function<void(const std::vector<AmrAction>&)> apply) {
    AmrCycleResult r;
    const std::size_t n = ind.size();
    if (n == 0) return r;
    if (level.size() != n) level.assign(n, 0);

    // Rank cells by indicator.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b){ return ind[a] > ind[b]; });

    const std::size_t nRef = std::size_t(double(n) * opt_.refineFraction);
    const std::size_t nCrs = std::size_t(double(n) * opt_.coarsenFraction);

    std::vector<AmrAction> action(n, AmrAction::Keep);
    for (std::size_t k = 0; k < nRef && k < n; ++k) {
        const std::size_t c = order[k];
        if (level[c] < opt_.maxLevel) action[c] = AmrAction::Refine;
    }
    for (std::size_t k = 0; k < nCrs && k < n; ++k) {
        const std::size_t c = order[n - 1 - k];
        if (level[c] > 0) action[c] = AmrAction::Coarsen;
    }

    // 2:1 enforcement: if a cell is to coarsen but a neighbour is +2 levels
    // higher, keep it; if a cell at level L has a neighbour at L-1 marked
    // for coarsen but our cell isn't marked for refine, leave neighbour Keep.
    if (opt_.enforceTwoToOne && !rowPtr.empty() && rowPtr.size() == n + 1) {
        bool changed = true;
        int  guard = 0;
        while (changed && guard++ < 4) {
            changed = false;
            for (std::size_t c = 0; c < n; ++c) {
                if (action[c] == AmrAction::Coarsen) {
                    for (std::size_t k = rowPtr[c]; k < rowPtr[c + 1]; ++k) {
                        const std::size_t nb = idx[k];
                        if (level[nb] > level[c]) { action[c] = AmrAction::Keep; changed = true; break; }
                    }
                }
                if (action[c] == AmrAction::Refine) {
                    for (std::size_t k = rowPtr[c]; k < rowPtr[c + 1]; ++k) {
                        const std::size_t nb = idx[k];
                        if (action[nb] == AmrAction::Coarsen
                            && std::uint8_t(level[c] + 1) > level[nb]) {
                            action[nb] = AmrAction::Keep;
                            changed = true;
                        }
                    }
                }
            }
        }
    }

    // Update level vector in-place so the caller's bookkeeping reflects the
    // mark decisions even if `apply` is a no-op (useful in tests).
    for (std::size_t c = 0; c < n; ++c) {
        if      (action[c] == AmrAction::Refine)  { ++level[c]; ++r.nRefined; }
        else if (action[c] == AmrAction::Coarsen) { --level[c]; ++r.nCoarsened; }
        else                                       { ++r.nKept; }
    }
    if (apply) apply(action);

    double sum = 0.0, mx = 0.0;
    for (auto v : ind) { sum += v; mx = std::max(mx, v); }
    r.maxIndicator  = mx;
    r.meanIndicator = sum / double(n);
    return r;
}

}  // namespace simall::amr
