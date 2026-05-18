// =============================================================================
// SimAll Beta - AMR Subsystem
// File   : src/amr/AmrApplier.hpp
// Week   : 18
//
// Adaptive Mesh Refinement orchestrator.  Bridges the existing
// AdaptiveRefinement *criteria* (gradient, hessian, Q-criterion, etc.)
// with the RefinementApplier *operator* (split / coarsen cells, rebuild
// face graph).  Performs a complete cycle in `apply_cycle()`:
//
//   1.  Evaluate per-cell error indicator η_i.
//   2.  Sort cells; mark top  P_refine  fraction for split,
//                    bottom P_coarsen  fraction for coarsen.
//   3.  Enforce level-difference cap (one-level rule).
//   4.  Hand the mark list to RefinementApplier.
//   5.  Re-balance partitions (callable hook).
// =============================================================================
#pragma once

#include <cstdint>
#include <functional>
#include <vector>

namespace simall::amr {

enum class AmrAction : std::uint8_t { Keep = 0, Refine = 1, Coarsen = 2 };

struct AmrOptions {
    double      refineFraction  = 0.10;
    double      coarsenFraction = 0.05;
    std::uint8_t maxLevel       = 4;
    bool         enforceTwoToOne = true;
};

struct AmrCycleResult {
    std::size_t nRefined  = 0;
    std::size_t nCoarsened = 0;
    std::size_t nKept      = 0;
    double      maxIndicator = 0.0;
    double      meanIndicator = 0.0;
};

class AmrApplier {
public:
    explicit AmrApplier(AmrOptions opt = {}) : opt_(opt) {}

    /// One adapt cycle.  The caller provides:
    ///   * `indicator`  — per-cell error metric (size = nCells);
    ///   * `level`      — current refinement level of each cell;
    ///   * `neighbours` — symmetric cell adjacency (flat CSR-style);
    ///   * `apply`      — hook that mutates the mesh according to actions.
    AmrCycleResult apply_cycle(const std::vector<double>&    indicator,
                                std::vector<std::uint8_t>&    level,
                                const std::vector<std::size_t>& neighRowPtr,
                                const std::vector<std::size_t>& neighIdx,
                                std::function<void(const std::vector<AmrAction>&)> apply);

private:
    AmrOptions opt_;
};

}  // namespace simall::amr
