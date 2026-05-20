// =============================================================================
// SimAll Beta - IBM Subsystem
// File   : src/ibm/DirectForcingIbm.hpp
// Week   : 18
//
// Direct-Forcing Immersed Boundary Method (Mohd-Yusof 1997, Uhlmann 2005).
// For each Lagrangian marker point on the immersed surface we:
//
//   1.  Interpolate the *predicted* Eulerian velocity u* to the marker
//       using a regularised δ-kernel (Roma 3-point kernel by default).
//   2.  Compute the forcing  f = (U_B − u*) / Δt   so that the Eulerian
//       velocity matches the body velocity U_B after one explicit step.
//   3.  Spread f back to the Eulerian grid via the same kernel.
//
// The implementation is grid-agnostic: it works against any cell-centred
// Eulerian field for which the caller supplies the cell centroid and the
// uniform spacing h.  Boundary-fitted accuracy is recovered by setting
// the kernel width to ~3·h.
// =============================================================================
#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace simall::ibm
{

using CellIdx = std::uint32_t;
using MarkerIdx = std::uint32_t;

struct MarkerPoint
{
    std::array<double, 3> position;
    std::array<double, 3> bodyVelocity;
    double area = 1.0;
};

struct EulerianCell
{
    std::array<double, 3> centroid;
    double volume = 1.0;
};

struct DirectForcingOptions
{
    double spacing = 1.0;     // characteristic h
    double kernelWidth = 3.0; // in units of h
    double dt = 1.0;
};

struct DirectForcingResult
{
    std::vector<std::array<double, 3>> markerForce;   // per-marker f
    std::vector<std::array<double, 3>> eulerianForce; // per-cell, scattered
};

[[nodiscard]] DirectForcingResult compute_direct_forcing(
    const std::vector<MarkerPoint>& markers,
    const std::vector<EulerianCell>& cells,
    const std::vector<std::array<double, 3>>& predictedVelocity,
    DirectForcingOptions opt);

} // namespace simall::ibm
