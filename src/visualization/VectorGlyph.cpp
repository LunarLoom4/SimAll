// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/VectorGlyph.cpp
// =============================================================================
#include "visualization/VectorGlyph.hpp"

namespace simall::visualization
{

GlyphSet VectorGlyph::sample(const VectorSampler& sampler,
                             const std::vector<util::Vec3d>& seeds,
                             const GlyphConfig& cfg)
{
    GlyphSet out;
    if (!sampler)
        return out;
    out.anchor.reserve(seeds.size());
    out.direction.reserve(seeds.size());
    out.magnitude.reserve(seeds.size());
    for (const auto& p : seeds) {
        if (out.anchor.size() >= cfg.maxAnchors)
            break;
        auto v = sampler(p);
        if (!v)
            continue;
        const double mag = v->norm();
        if (mag <= cfg.minMagnitude)
            continue;
        util::Vec3d dir = *v;
        if (cfg.normalizeDirection && mag > 0.0) {
            dir = util::Vec3d{dir.x / mag, dir.y / mag, dir.z / mag};
        }
        out.anchor.push_back(p);
        out.direction.push_back(dir);
        out.magnitude.push_back(mag);
    }
    return out;
}

GlyphSet VectorGlyph::sample_lattice(const VectorSampler& sampler,
                                     const util::BoundingBox& box,
                                     std::uint32_t nx,
                                     std::uint32_t ny,
                                     std::uint32_t nz,
                                     const GlyphConfig& cfg)
{
    std::vector<util::Vec3d> seeds;
    if (!box.valid() || nx == 0 || ny == 0 || nz == 0) {
        return VectorGlyph::sample(sampler, seeds, cfg);
    }
    seeds.reserve(static_cast<std::size_t>(nx) * ny * nz);
    const util::Vec3d e = box.extent();
    const double dx = nx > 1 ? e.x / static_cast<double>(nx - 1) : 0.0;
    const double dy = ny > 1 ? e.y / static_cast<double>(ny - 1) : 0.0;
    const double dz = nz > 1 ? e.z / static_cast<double>(nz - 1) : 0.0;
    for (std::uint32_t k = 0; k < nz; ++k) {
        for (std::uint32_t j = 0; j < ny; ++j) {
            for (std::uint32_t i = 0; i < nx; ++i) {
                seeds.push_back({box.min.x + i * dx, box.min.y + j * dy, box.min.z + k * dz});
            }
        }
    }
    return VectorGlyph::sample(sampler, seeds, cfg);
}

} // namespace simall::visualization
