// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/VolumeRaycast.cpp
// =============================================================================
#include "visualization/VolumeRaycast.hpp"

#include <algorithm>
#include <cmath>

namespace simall::visualization {

namespace {

inline util::Vec3d safe_normalize(const util::Vec3d& v, const util::Vec3d& fb) {
    const double n = v.norm();
    return (n > 1e-30) ? util::Vec3d{v.x/n, v.y/n, v.z/n} : fb;
}

/// Trilinear sample at world position p.  Returns 0 when out of domain.
inline double sample_trilinear(const RegularGrid3& g, const util::Vec3d& p) {
    const double fx = (p.x - g.origin.x) / g.spacing.x;
    const double fy = (p.y - g.origin.y) / g.spacing.y;
    const double fz = (p.z - g.origin.z) / g.spacing.z;
    if (fx < 0 || fy < 0 || fz < 0) return 0.0;
    if (fx > g.nx - 1 || fy > g.ny - 1 || fz > g.nz - 1) return 0.0;
    const std::uint32_t i0 = static_cast<std::uint32_t>(fx);
    const std::uint32_t j0 = static_cast<std::uint32_t>(fy);
    const std::uint32_t k0 = static_cast<std::uint32_t>(fz);
    const std::uint32_t i1 = std::min(i0 + 1, g.nx - 1);
    const std::uint32_t j1 = std::min(j0 + 1, g.ny - 1);
    const std::uint32_t k1 = std::min(k0 + 1, g.nz - 1);
    const double tx = fx - i0, ty = fy - j0, tz = fz - k0;
    auto V = [&](std::uint32_t i, std::uint32_t j, std::uint32_t k) {
        return static_cast<double>(g.values[(static_cast<std::size_t>(k) * g.ny + j) * g.nx + i]);
    };
    const double c00 = V(i0,j0,k0)*(1-tx) + V(i1,j0,k0)*tx;
    const double c10 = V(i0,j1,k0)*(1-tx) + V(i1,j1,k0)*tx;
    const double c01 = V(i0,j0,k1)*(1-tx) + V(i1,j0,k1)*tx;
    const double c11 = V(i0,j1,k1)*(1-tx) + V(i1,j1,k1)*tx;
    const double c0  = c00*(1-ty) + c10*ty;
    const double c1  = c01*(1-ty) + c11*ty;
    return c0*(1-tz) + c1*tz;
}

/// Slab intersection of a ray with an axis-aligned box.  Returns (tNear, tFar);
/// out is invalid (tNear > tFar) when the ray misses.
struct Slab { double tNear; double tFar; };
Slab intersect_aabb(const util::Vec3d& o, const util::Vec3d& d,
                    const util::BoundingBox& b) {
    double tN = -1e30, tF = 1e30;
    for (int ax = 0; ax < 3; ++ax) {
        const double oo = (ax == 0 ? o.x : ax == 1 ? o.y : o.z);
        const double dd = (ax == 0 ? d.x : ax == 1 ? d.y : d.z);
        const double bMin = (ax == 0 ? b.min.x : ax == 1 ? b.min.y : b.min.z);
        const double bMax = (ax == 0 ? b.max.x : ax == 1 ? b.max.y : b.max.z);
        if (std::abs(dd) < 1e-30) {
            if (oo < bMin || oo > bMax) return {1.0, -1.0};
            continue;
        }
        double t1 = (bMin - oo) / dd;
        double t2 = (bMax - oo) / dd;
        if (t1 > t2) std::swap(t1, t2);
        if (t1 > tN) tN = t1;
        if (t2 < tF) tF = t2;
        if (tN > tF) return {1.0, -1.0};
    }
    return {tN, tF};
}

}  // namespace

Image VolumeRaycast::render(const RegularGrid3& g,
                             const Camera& cam,
                             const TransferFunction& tf,
                             const RaycastConfig& cfg)
{
    Image img;
    if (!g.valid()) return img;
    img.resize(cam.viewportWidth, cam.viewportHeight);

    const util::BoundingBox box = g.bbox();
    const double ds = (cfg.sampleStep > 0.0) ? cfg.sampleStep
                    : 0.5 * std::min({g.spacing.x, g.spacing.y, g.spacing.z});

    // Build orthonormal camera basis.
    util::Vec3d fwd   = safe_normalize(cam.focalPoint - cam.position, {0,0,-1});
    util::Vec3d right = safe_normalize(fwd.cross(cam.up), {1,0,0});
    util::Vec3d up    = safe_normalize(right.cross(fwd),  {0,1,0});

    const double aspect = cam.aspect();
    const double halfH  = cam.orthographic
        ? cam.orthoHeight * 0.5
        : std::tan(cam.fovYDegrees * 0.5 * 3.14159265358979323846 / 180.0);
    const double halfW  = halfH * aspect;

    for (std::uint32_t y = 0; y < cam.viewportHeight; ++y) {
        for (std::uint32_t x = 0; x < cam.viewportWidth; ++x) {
            const double sx = (2.0 * (x + 0.5) / cam.viewportWidth  - 1.0);
            const double sy = (1.0 - 2.0 * (y + 0.5) / cam.viewportHeight);
            util::Vec3d origin;
            util::Vec3d dir;
            if (cam.orthographic) {
                origin = cam.position + right * (sx * halfW) + up * (sy * halfH);
                dir    = fwd;
            } else {
                origin = cam.position;
                dir = safe_normalize(fwd + right * (sx * halfW) + up * (sy * halfH),
                                     {0,0,-1});
            }
            Slab s = intersect_aabb(origin, dir, box);

            // Initialise pixel with background colour.
            Color4 px = cfg.background;

            if (s.tNear <= s.tFar) {
                const double tStart = std::max(s.tNear, 0.0);
                Color4 acc{0,0,0,0};
                for (double t = tStart; t <= s.tFar; t += ds) {
                    util::Vec3d p = origin + dir * t;
                    const double sample = sample_trilinear(g, p);
                    Color4 c = tf.sample(sample);
                    const float a = static_cast<float>(
                        std::clamp(c[3] * cfg.opacityScale, 0.0, 1.0));
                    const float oneMinus = 1.0f - acc[3];
                    acc[0] += oneMinus * a * c[0];
                    acc[1] += oneMinus * a * c[1];
                    acc[2] += oneMinus * a * c[2];
                    acc[3] += oneMinus * a;
                    if (cfg.useEarlyTermination &&
                        acc[3] >= static_cast<float>(cfg.earlyOpacity)) break;
                }
                // Composite over background.
                const float oneMinus = 1.0f - acc[3];
                px[0] = acc[0] + oneMinus * cfg.background[0];
                px[1] = acc[1] + oneMinus * cfg.background[1];
                px[2] = acc[2] + oneMinus * cfg.background[2];
                px[3] = acc[3] + oneMinus * cfg.background[3];
            }

            const std::size_t off = (static_cast<std::size_t>(y) * img.width + x) * 4u;
            img.pixels[off + 0] = static_cast<std::uint8_t>(std::clamp(px[0]*255.0f, 0.0f, 255.0f));
            img.pixels[off + 1] = static_cast<std::uint8_t>(std::clamp(px[1]*255.0f, 0.0f, 255.0f));
            img.pixels[off + 2] = static_cast<std::uint8_t>(std::clamp(px[2]*255.0f, 0.0f, 255.0f));
            img.pixels[off + 3] = static_cast<std::uint8_t>(std::clamp(px[3]*255.0f, 0.0f, 255.0f));
        }
    }
    return img;
}

}  // namespace simall::visualization
