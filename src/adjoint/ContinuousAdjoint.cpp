// =============================================================================
// SimAll Beta - Adjoint Subsystem
// File   : src/adjoint/ContinuousAdjoint.cpp
// =============================================================================
#include "adjoint/ContinuousAdjoint.hpp"

#include <algorithm>
#include <cmath>

namespace simall::adjoint
{

std::vector<double> assemble_shape_sensitivity(const std::vector<SurfaceFace>& faces,
                                               const std::vector<AdjointFlowSample>& s)
{
    std::vector<double> out(faces.size(), 0.0);
    const std::size_t n = std::min(faces.size(), s.size());
    for (std::size_t f = 0; f < n; ++f) {
        const auto& fc = faces[f];
        const auto& a = s[f];
        // Wall-shape sensitivity for impermeable wall:
        //   G = -p_ψ + 2 μ (∂u_ψ · n / ∂n)   simplified to  G = -p_ψ + 2μ(u_ψ · n)
        // when the adjoint flow gradient is replaced by the bulk adjoint velocity.
        const double udotn =
            a.uPsi[0] * fc.normal[0] + a.uPsi[1] * fc.normal[1] + a.uPsi[2] * fc.normal[2];
        out[f] = (-a.pPsi + 2.0 * a.mu * udotn) * fc.area;
    }
    return out;
}

std::vector<double> assemble_volume_sensitivity(const std::vector<double>& p,
                                                const std::vector<double>& pPsi)
{
    const std::size_t n = std::min(p.size(), pPsi.size());
    std::vector<double> out(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        out[i] = -pPsi[i] * p[i];
    return out;
}

std::vector<double> accumulate_node_sensitivity(std::size_t nNodes,
                                                const std::vector<SurfaceFace>& faces,
                                                const std::vector<double>& faceSens)
{
    std::vector<double> out(nNodes, 0.0);
    std::vector<double> w(nNodes, 0.0);
    const std::size_t n = std::min(faces.size(), faceSens.size());
    for (std::size_t f = 0; f < n; ++f) {
        const auto& fc = faces[f];
        const double share = (fc.nNodes > 0) ? 1.0 / double(fc.nNodes) : 0.0;
        for (std::uint8_t k = 0; k < fc.nNodes; ++k) {
            const NodeIdx ni = fc.nodes[k];
            if (ni < nNodes) {
                out[ni] += faceSens[f] * share;
                w[ni] += fc.area * share;
            }
        }
    }
    for (std::size_t i = 0; i < nNodes; ++i)
        if (w[i] > 1e-300)
            out[i] /= w[i];
    return out;
}

} // namespace simall::adjoint
