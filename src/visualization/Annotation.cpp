// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/Annotation.cpp
// =============================================================================
#include "visualization/Annotation.hpp"

#include <cmath>

namespace simall::visualization
{

namespace
{
inline util::Vec3d safe_normalize(const util::Vec3d& v, const util::Vec3d& fb)
{
    const double n = v.norm();
    return (n > 1e-30) ? util::Vec3d{v.x / n, v.y / n, v.z / n} : fb;
}
} // namespace

std::uint64_t AnnotationLayer::add(const Annotation& a)
{
    std::lock_guard lock(mu_);
    Annotation copy = a;
    copy.id = nextId_++;
    items_.emplace(copy.id, copy);
    return copy.id;
}

void AnnotationLayer::remove(std::uint64_t id)
{
    std::lock_guard lock(mu_);
    items_.erase(id);
}

void AnnotationLayer::clear()
{
    std::lock_guard lock(mu_);
    items_.clear();
}

void AnnotationLayer::set_visible(std::uint64_t id, bool v)
{
    std::lock_guard lock(mu_);
    auto it = items_.find(id);
    if (it != items_.end())
        it->second.visible = v;
}

std::vector<Annotation> AnnotationLayer::snapshot() const
{
    std::lock_guard lock(mu_);
    std::vector<Annotation> out;
    out.reserve(items_.size());
    for (auto& [_, a] : items_)
        out.push_back(a);
    return out;
}

std::size_t AnnotationLayer::size() const noexcept
{
    std::lock_guard lock(mu_);
    return items_.size();
}

std::optional<std::pair<double, double>> AnnotationLayer::project_world(const Camera& cam,
                                                                        const util::Vec3d& world)
{
    util::Vec3d fwd = safe_normalize(cam.focalPoint - cam.position, {0, 0, -1});
    util::Vec3d right = safe_normalize(fwd.cross(cam.up), {1, 0, 0});
    util::Vec3d up = safe_normalize(right.cross(fwd), {0, 1, 0});

    util::Vec3d rel = world - cam.position;
    const double zCam = fwd.dot(rel);
    if (zCam <= 1e-9)
        return std::nullopt; // behind / at camera

    const double xCam = right.dot(rel);
    const double yCam = up.dot(rel);

    const double aspect = cam.aspect();
    const double halfH =
        cam.orthographic ? cam.orthoHeight * 0.5
                         : zCam * std::tan(cam.fovYDegrees * 0.5 * 3.14159265358979323846 / 180.0);
    const double halfW = halfH * aspect;
    if (halfH <= 0.0)
        return std::nullopt;

    const double ndcX = xCam / halfW;
    const double ndcY = yCam / halfH;
    const double px = (ndcX * 0.5 + 0.5) * static_cast<double>(cam.viewportWidth);
    const double py = (1.0 - (ndcY * 0.5 + 0.5)) * static_cast<double>(cam.viewportHeight);
    return std::make_pair(px, py);
}

} // namespace simall::visualization
