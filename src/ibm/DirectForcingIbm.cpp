// =============================================================================
// SimAll Beta - IBM Subsystem
// File   : src/ibm/DirectForcingIbm.cpp
// =============================================================================
#include "ibm/DirectForcingIbm.hpp"

#include <algorithm>
#include <cmath>

namespace simall::ibm
{

namespace
{

// Roma 3-point regularised δ kernel: support |r/h| ≤ 1.5, ∫δ dx = 1.
double roma_kernel(double r, double h)
{
    const double q = std::abs(r) / std::max(h, 1e-30);
    if (q < 0.5) {
        return (1.0 / (3.0 * h)) * (1.0 + std::sqrt(std::max(1.0 - 3.0 * q * q, 0.0)));
    } else if (q < 1.5) {
        const double t = 5.0 - 3.0 * q;
        const double inside = std::max(1.0 - 3.0 * (1.0 - q) * (1.0 - q), 0.0);
        return (1.0 / (6.0 * h)) * (t - std::sqrt(inside));
    }
    return 0.0;
}

double kernel3d(const std::array<double, 3>& a, const std::array<double, 3>& b, double h)
{
    return roma_kernel(a[0] - b[0], h) * roma_kernel(a[1] - b[1], h) * roma_kernel(a[2] - b[2], h);
}

} // namespace

DirectForcingResult compute_direct_forcing(const std::vector<MarkerPoint>& markers,
                                           const std::vector<EulerianCell>& cells,
                                           const std::vector<std::array<double, 3>>& uStar,
                                           DirectForcingOptions opt)
{
    DirectForcingResult r;
    r.markerForce.assign(markers.size(), {0, 0, 0});
    r.eulerianForce.assign(cells.size(), {0, 0, 0});
    if (markers.empty() || cells.empty() || uStar.size() != cells.size())
        return r;
    const double h = std::max(opt.spacing, 1e-30);
    const double half = opt.kernelWidth * h;
    const double dt = std::max(opt.dt, 1e-30);

    // Step 1+2: interpolate u* at each marker, compute force.
    for (std::size_t m = 0; m < markers.size(); ++m) {
        std::array<double, 3> uInterp = {0, 0, 0};
        for (std::size_t c = 0; c < cells.size(); ++c) {
            if (std::abs(markers[m].position[0] - cells[c].centroid[0]) > half)
                continue;
            if (std::abs(markers[m].position[1] - cells[c].centroid[1]) > half)
                continue;
            if (std::abs(markers[m].position[2] - cells[c].centroid[2]) > half)
                continue;
            const double w = kernel3d(markers[m].position, cells[c].centroid, h) * cells[c].volume;
            uInterp[0] += w * uStar[c][0];
            uInterp[1] += w * uStar[c][1];
            uInterp[2] += w * uStar[c][2];
        }
        for (int k = 0; k < 3; ++k)
            r.markerForce[m][k] = (markers[m].bodyVelocity[k] - uInterp[k]) / dt;
    }

    // Step 3: spread marker forces back to grid.
    for (std::size_t m = 0; m < markers.size(); ++m) {
        for (std::size_t c = 0; c < cells.size(); ++c) {
            if (std::abs(markers[m].position[0] - cells[c].centroid[0]) > half)
                continue;
            if (std::abs(markers[m].position[1] - cells[c].centroid[1]) > half)
                continue;
            if (std::abs(markers[m].position[2] - cells[c].centroid[2]) > half)
                continue;
            const double w = kernel3d(markers[m].position, cells[c].centroid, h) * markers[m].area;
            for (int k = 0; k < 3; ++k)
                r.eulerianForce[c][k] += w * r.markerForce[m][k];
        }
    }
    return r;
}

} // namespace simall::ibm
