// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/CameraController.cpp
// =============================================================================
#include "visualization/CameraController.hpp"

#include <algorithm>
#include <cmath>

namespace simall::visualization {

namespace {
inline util::Vec3d safe_normalize(const util::Vec3d& v, const util::Vec3d& fb) {
    const double n = v.norm();
    return (n > 1e-30) ? util::Vec3d{v.x / n, v.y / n, v.z / n} : fb;
}
}  // namespace

CameraController::CameraController() = default;

void CameraController::set_camera(const Camera& c) {
    std::lock_guard lock(mu_); cam_ = c;
}

Camera CameraController::camera() const {
    std::lock_guard lock(mu_); return cam_;
}

void CameraController::set_mode(InteractionMode m) {
    std::lock_guard lock(mu_); mode_ = m;
}

void CameraController::on_drag(double dx, double dy) {
    std::lock_guard lock(mu_);
    switch (mode_) {
        case InteractionMode::Orbit: apply_orbit(dx, dy); break;
        case InteractionMode::Pan:   apply_pan  (dx, dy); break;
        case InteractionMode::Dolly: apply_dolly(dy);     break;
        case InteractionMode::Walk:                       break;  // ignore
    }
}

void CameraController::on_wheel(double ticks) {
    std::lock_guard lock(mu_);
    // Negate so positive ticks bring the camera closer regardless of mode.
    apply_dolly(-ticks * 50.0);   // treat wheel as ~50 px equivalent
}

void CameraController::on_key(double forward, double right, double up) {
    std::lock_guard lock(mu_);
    if (mode_ != InteractionMode::Walk) return;
    apply_walk(forward, right, up);
}

void CameraController::fit_to(const util::BoundingBox& box, double margin) {
    std::lock_guard lock(mu_);
    if (!box.valid()) return;
    const util::Vec3d c = box.center();
    const util::Vec3d e = box.extent();
    const double diag = std::sqrt(e.x*e.x + e.y*e.y + e.z*e.z) * margin;

    util::Vec3d viewDir = safe_normalize(cam_.focalPoint - cam_.position,
                                          {0, 0, -1});
    const double halfFov = cam_.fovYDegrees * 0.5 * 3.14159265358979323846 / 180.0;
    const double dist = (cam_.orthographic || halfFov <= 0.0)
                            ? diag
                            : 0.5 * diag / std::tan(halfFov);

    cam_.focalPoint = c;
    cam_.position   = c - viewDir * dist;
    cam_.orthoHeight = diag;
    // Loosen near/far around the new framing.
    cam_.nearPlane = std::max(1e-4, 0.01 * diag);
    cam_.farPlane  = std::max(cam_.nearPlane * 10.0, 100.0 * diag);
}

void CameraController::snap_axis(util::Vec3d axis, util::Vec3d up) {
    std::lock_guard lock(mu_);
    axis = safe_normalize(axis, {0, 0, 1});
    up   = safe_normalize(up,   {0, 1, 0});
    const double dist = (cam_.position - cam_.focalPoint).norm();
    cam_.position = cam_.focalPoint - axis * dist;
    cam_.up       = up;
}

void CameraController::save_view(const std::string& name) {
    std::lock_guard lock(mu_);
    saved_[name] = cam_;
}

bool CameraController::restore_view(const std::string& name) {
    std::lock_guard lock(mu_);
    auto it = saved_.find(name);
    if (it == saved_.end()) return false;
    cam_ = it->second;
    return true;
}

std::vector<std::string> CameraController::saved_views() const {
    std::lock_guard lock(mu_);
    std::vector<std::string> out; out.reserve(saved_.size());
    for (auto& [k, _] : saved_) out.push_back(k);
    return out;
}

// ---------------------------------------------------------------------------
//  Apply primitives (caller already holds mu_)
// ---------------------------------------------------------------------------
void CameraController::apply_orbit(double dxPx, double dyPx) {
    const double azim = -dxPx * orbitSens_;
    const double elev = -dyPx * orbitSens_;

    util::Vec3d offset = cam_.position - cam_.focalPoint;
    util::Vec3d up     = safe_normalize(cam_.up, {0, 1, 0});
    util::Vec3d right  = safe_normalize(offset.cross(up), {1, 0, 0});

    // Rodrigues rotation around up by azim.
    auto rotate = [](const util::Vec3d& v, util::Vec3d k, double th) {
        k = safe_normalize(k, {0, 1, 0});
        const double c = std::cos(th), s = std::sin(th);
        const double dot = k.dot(v);
        return v * c
             + k.cross(v) * s
             + k * (dot * (1.0 - c));
    };
    offset = rotate(offset, up,    azim);
    right  = rotate(right,  up,    azim);
    offset = rotate(offset, right, elev);

    cam_.position = cam_.focalPoint + offset;
    cam_.up = safe_normalize(right.cross(offset) * -1.0, up);
}

void CameraController::apply_pan(double dxPx, double dyPx) {
    util::Vec3d offset = cam_.position - cam_.focalPoint;
    const double dist  = offset.norm();
    util::Vec3d viewDir = safe_normalize(offset, {0, 0, 1});
    util::Vec3d up      = safe_normalize(cam_.up, {0, 1, 0});
    util::Vec3d right   = safe_normalize(viewDir.cross(up), {1, 0, 0});
    up = safe_normalize(right.cross(viewDir), up);

    const double k = panSens_ * std::max(1.0, dist);
    util::Vec3d delta = right * (-dxPx * k) + up * (dyPx * k);
    cam_.position   = cam_.position   + delta;
    cam_.focalPoint = cam_.focalPoint + delta;
}

void CameraController::apply_dolly(double dyPxOrTicks) {
    // Positive dyPx (mouse down) → zoom out, negative → in.  Match VTK feel.
    const double factor = std::pow(dollyFactor_, dyPxOrTicks * 0.01);
    util::Vec3d offset = cam_.position - cam_.focalPoint;
    offset = offset * factor;
    if (offset.norm() < 1e-6) offset = util::Vec3d{1e-6, 0, 0};
    cam_.position = cam_.focalPoint + offset;
    if (cam_.orthographic) cam_.orthoHeight *= factor;
}

void CameraController::apply_walk(double forward, double right, double up) {
    util::Vec3d viewDir = safe_normalize(cam_.focalPoint - cam_.position,
                                          {0, 0, -1});
    util::Vec3d upV     = safe_normalize(cam_.up, {0, 1, 0});
    util::Vec3d rt      = safe_normalize(viewDir.cross(upV), {1, 0, 0});
    util::Vec3d delta   = viewDir * forward + rt * right + upV * up;
    cam_.position   = cam_.position   + delta;
    cam_.focalPoint = cam_.focalPoint + delta;
}

}  // namespace simall::visualization
