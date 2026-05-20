// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/PopulationBalance.hpp
// Phase  : 12 — Population Balance Model (PBM) via Quadrature Method of
// Moments (QMOM, McGraw 1997).
//
// Number-density transport for N moments m_k = ∫ L^k n(L) dL,  k = 0..2N-1:
//
//   ∂m_k/∂t + ∇·(u m_k) = B_k - D_k + G_k
//
// where B_k, D_k, G_k are the moment source terms from particle birth
// (coalescence/breakage), death, and growth respectively.
//
// The product-difference (PD) algorithm of Gordon (1968) recovers the
// abscissae L_i and weights w_i from the moments so that the source
// terms can be evaluated by Gauss-Christoffel quadrature.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <array>
#include <functional>
#include <string>
#include <vector>

namespace simall::multiphase
{

struct PbmProps
{
    int numMoments = 4;       // even, ≥ 2; ⇒ N = numMoments/2 nodes
    double growthRateG = 0.0; // dL/dt (constant; user can refine)
    /// Aggregation kernel β(L_i, L_j) [m³/s].
    std::function<double(double, double)> aggregationKernel;
    /// Breakage frequency g(L) [1/s] and fragment distribution PDF β(L, L′).
    std::function<double(double)> breakageFrequency;
    std::function<double(double, double)> daughterDistribution;
};

class PopulationBalance
{
public:
    void initialize(const meshing::Mesh& mesh, solver::FieldRegistry& fields, PbmProps props);

    /// Advance moments by `dt` using a first-order forward Euler integration
    /// of the moment source terms (transport handled separately by a
    /// ScalarTransport call on each "mN" field by the outer solver).
    void integrate_sources(double dt);

    /// Recover (weights, abscissae) at cell `c`. Returns true on success.
    bool quadrature_nodes(std::size_t c,
                          std::vector<double>& weights,
                          std::vector<double>& abscissae) const;

private:
    /// Wheeler / product-difference algorithm. m has size 2N.
    static bool pd_algorithm(const std::vector<double>& m,
                             std::vector<double>& w,
                             std::vector<double>& L);

    std::string moment_name(int k) const;

    const meshing::Mesh* mesh_ = nullptr;
    solver::FieldRegistry* fields_ = nullptr;
    PbmProps props_{};
};

} // namespace simall::multiphase
