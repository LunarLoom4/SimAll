// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/CameraController.hpp
// Phase  : 3 / Week 14 (scene-management trio: CameraController)
//
// State machine that maps mouse-style deltas into Camera updates.  No Qt /
// OpenGL dependency — pixel/wheel/key deltas come in, Camera comes out, the
// GUI binds those signals.  Supports orbit (turntable), pan, dolly,
// fit-to-bbox, walk mode, and named-view bookmarking.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <mutex>
#include <string>
#include <unordered_map>

namespace simall::visualization {

enum class InteractionMode : std::uint8_t {
    Orbit,        // LMB drag = azimuth/elevation around focal point
    Pan,          // MMB drag = translate camera + focal point in screen XY
    Dolly,        // RMB drag = move camera along view axis (zoom)
    Walk          // WASD-style; arrow direction in deltas
};

class CameraController {
public:
    CameraController();

    void  set_camera(const Camera& c);
    Camera camera() const;

    void  set_mode(InteractionMode m);
    InteractionMode mode() const noexcept { return mode_; }

    /// Pixel deltas from the GUI.  dx/dy are screen-pixel deltas; the
    /// controller scales them by sensitivity and the current viewport size.
    void on_drag(double dx, double dy);
    void on_wheel(double ticks);              // positive = zoom in
    void on_key  (double forward, double right, double up);   // walk-mode

    /// Frame the camera onto the given world bbox, leaving a margin (1.05x
    /// default).  Keeps current orientation.
    void fit_to(const util::BoundingBox& box, double margin = 1.05);

    /// Snap to a canonical orientation (e.g. +X axis looking towards origin,
    /// up = +Z).
    void snap_axis(util::Vec3d axis, util::Vec3d up);

    /// Named views.
    void                       save_view(const std::string& name);
    bool                       restore_view(const std::string& name);
    std::vector<std::string>   saved_views() const;

    /// Sensitivities (radians/pixel for orbit, world/pixel for pan,
    /// log-zoom per tick for dolly).
    void set_orbit_sensitivity (double radPerPixel) { orbitSens_ = radPerPixel; }
    void set_pan_sensitivity   (double worldPerPixel) { panSens_  = worldPerPixel; }
    void set_dolly_sensitivity (double factorPerTick) { dollyFactor_ = factorPerTick; }

private:
    void apply_orbit(double dxPx, double dyPx);
    void apply_pan  (double dxPx, double dyPx);
    void apply_dolly(double dyPxOrTicks);
    void apply_walk (double forward, double right, double up);

    mutable std::mutex                            mu_;
    Camera                                        cam_;
    InteractionMode                               mode_       = InteractionMode::Orbit;
    double                                        orbitSens_  = 0.01;
    double                                        panSens_    = 0.005;
    double                                        dollyFactor_ = 0.9;
    std::unordered_map<std::string, Camera>       saved_;
};

}  // namespace simall::visualization
