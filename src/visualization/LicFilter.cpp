// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/LicFilter.cpp
// =============================================================================
#include "visualization/LicFilter.hpp"

#include <algorithm>
#include <cmath>

namespace simall::visualization
{

namespace
{

/// Trilinear sample, but in 2D.  Returns (vx, vy) in grid-cell units per
/// world-unit step (we keep the API in grid coords).
struct V2
{
    double x;
    double y;
};

inline V2 sample_vec(const RegularGrid2& g, double fx, double fy)
{
    if (fx < 0)
        fx = 0;
    if (fy < 0)
        fy = 0;
    if (fx > g.nx - 1)
        fx = g.nx - 1;
    if (fy > g.ny - 1)
        fy = g.ny - 1;
    const std::uint32_t i0 = static_cast<std::uint32_t>(fx);
    const std::uint32_t j0 = static_cast<std::uint32_t>(fy);
    const std::uint32_t i1 = std::min<std::uint32_t>(i0 + 1, g.nx - 1);
    const std::uint32_t j1 = std::min<std::uint32_t>(j0 + 1, g.ny - 1);
    const double tx = fx - i0, ty = fy - j0;
    auto idx = [&](std::uint32_t i, std::uint32_t j) {
        return static_cast<std::size_t>(j) * g.nx + i;
    };
    const double v00x = g.vx[idx(i0, j0)], v00y = g.vy[idx(i0, j0)];
    const double v10x = g.vx[idx(i1, j0)], v10y = g.vy[idx(i1, j0)];
    const double v01x = g.vx[idx(i0, j1)], v01y = g.vy[idx(i0, j1)];
    const double v11x = g.vx[idx(i1, j1)], v11y = g.vy[idx(i1, j1)];
    const double x0x = v00x * (1 - tx) + v10x * tx;
    const double x1x = v01x * (1 - tx) + v11x * tx;
    const double x0y = v00y * (1 - tx) + v10y * tx;
    const double x1y = v01y * (1 - tx) + v11y * tx;
    return {x0x * (1 - ty) + x1x * ty, x0y * (1 - ty) + x1y * ty};
}

/// Tiny SplitMix64 RNG so the noise field is deterministic across platforms.
inline std::uint64_t splitmix(std::uint64_t& s)
{
    s += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = s;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

} // namespace

Image LicFilter::generate(const RegularGrid2& g, const LicConfig& cfgIn)
{
    Image img;
    if (!g.valid())
        return img;
    LicConfig cfg = cfgIn;
    if (cfg.width == 0)
        cfg.width = g.nx;
    if (cfg.height == 0)
        cfg.height = g.ny;
    img.resize(cfg.width, cfg.height);

    // Generate per-pixel noise texture once.
    std::vector<float> noise(static_cast<std::size_t>(cfg.width) * cfg.height);
    {
        std::uint64_t s = cfg.noiseSeed ? cfg.noiseSeed : 1ULL;
        for (auto& n : noise) {
            n = static_cast<float>(splitmix(s) >> 40) / 16777216.0f; // [0,1)
        }
    }

    auto noise_at = [&](double fx, double fy) -> double {
        if (fx < 0)
            fx = 0;
        if (fy < 0)
            fy = 0;
        if (fx > cfg.width - 1)
            fx = cfg.width - 1;
        if (fy > cfg.height - 1)
            fy = cfg.height - 1;
        const std::uint32_t i = static_cast<std::uint32_t>(fx);
        const std::uint32_t j = static_cast<std::uint32_t>(fy);
        return noise[static_cast<std::size_t>(j) * cfg.width + i];
    };

    // Map output (i,j) → grid (gx,gy)
    const double sx = static_cast<double>(g.nx - 1) / std::max<double>(cfg.width - 1, 1);
    const double sy = static_cast<double>(g.ny - 1) / std::max<double>(cfg.height - 1, 1);

    for (std::uint32_t j = 0; j < cfg.height; ++j) {
        for (std::uint32_t i = 0; i < cfg.width; ++i) {
            double gx = i * sx;
            double gy = j * sy;
            double acc = noise_at(i, j);
            double weight = 1.0;

            // Walk forward + backward in field space.
            for (int dir = 0; dir < 2; ++dir) {
                double fx = gx, fy = gy;
                const double sign = (dir == 0) ? 1.0 : -1.0;
                for (std::uint32_t k = 0; k < cfg.stepsPerSide; ++k) {
                    V2 v = sample_vec(g, fx, fy);
                    const double mag = std::sqrt(v.x * v.x + v.y * v.y);
                    if (mag < 1e-12)
                        break;
                    fx += sign * cfg.stepSize * v.x / mag;
                    fy += sign * cfg.stepSize * v.y / mag;
                    const double ix = fx / sx;
                    const double iy = fy / sy;
                    if (ix < 0 || iy < 0 || ix > cfg.width - 1 || iy > cfg.height - 1)
                        break;
                    acc += noise_at(ix, iy);
                    weight += 1.0;
                }
            }
            const double normalized = acc / weight;
            const std::uint8_t b =
                static_cast<std::uint8_t>(std::clamp(normalized, 0.0, 1.0) * 255.0);
            const std::size_t off = (static_cast<std::size_t>(j) * cfg.width + i) * 4u;
            img.pixels[off + 0] = b;
            img.pixels[off + 1] = b;
            img.pixels[off + 2] = b;
            img.pixels[off + 3] = 255;
        }
    }
    return img;
}

} // namespace simall::visualization
