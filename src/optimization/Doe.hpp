// =============================================================================
// SimAll Beta - Optimization Subsystem
// File   : src/optimization/Doe.hpp
// Week   : 18
//
// Design of Experiments samplers.  All return an N × dim matrix in
// row-major order whose rows are space-filling design points scaled into
// the user-supplied box [xLow, xUp].
//
//   * Full Factorial (k levels per dimension)
//   * Latin Hypercube Sampling (LHS) with optimised maximin distance
//   * Sobol low-discrepancy sequence (up to 6 dimensions natively;
//     higher-D extension by direction-number table)
// =============================================================================
#pragma once

#include <cstdint>
#include <vector>

namespace simall::optimization {

[[nodiscard]] std::vector<std::vector<double>> doe_full_factorial(
        const std::vector<double>& xLow,
        const std::vector<double>& xUp,
        std::size_t                levels);

[[nodiscard]] std::vector<std::vector<double>> doe_latin_hypercube(
        const std::vector<double>& xLow,
        const std::vector<double>& xUp,
        std::size_t                nSamples,
        std::uint64_t              seed = 0xA5A5A5A5);

[[nodiscard]] std::vector<std::vector<double>> doe_sobol(
        const std::vector<double>& xLow,
        const std::vector<double>& xUp,
        std::size_t                nSamples);

}  // namespace simall::optimization
