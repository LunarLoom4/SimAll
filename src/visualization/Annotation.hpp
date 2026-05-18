// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/Annotation.hpp
// Phase  : 3 / Week 14 (overlays: text annotations)
//
// Text labels positioned in world or screen space.  Owns nothing graphical
// — the back-end renderer reads the AnnotationLayer at frame time and
// renders text using its own glyph atlas.  Provides a CPU helper to
// project world-anchored annotations to screen coordinates via the camera.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace simall::visualization {

enum class AnchorSpace : std::uint8_t { World, ScreenPixels, ScreenNormalized };

enum class TextAlign : std::uint8_t {
    TopLeft, TopCenter, TopRight,
    MiddleLeft, MiddleCenter, MiddleRight,
    BottomLeft, BottomCenter, BottomRight
};

struct Annotation {
    std::uint64_t  id          = 0;
    std::string    text;
    util::Vec3d    anchor      {0, 0, 0};
    AnchorSpace    space       = AnchorSpace::World;
    TextAlign      align       = TextAlign::MiddleCenter;
    double         fontSizePx  = 14.0;
    Color4         color       {1, 1, 1, 1};
    Color4         background  {0, 0, 0, 0};
    bool           visible     = true;
};

class AnnotationLayer {
public:
    std::uint64_t add(const Annotation& a);
    void          remove(std::uint64_t id);
    void          clear();
    void          set_visible(std::uint64_t id, bool v);

    std::vector<Annotation> snapshot() const;
    std::size_t             size() const noexcept;

    /// Project a world-space anchor through the given Camera.
    /// Returns (xPx, yPx) using origin top-left of the viewport.
    /// nullopt means the anchor is behind the camera.
    static std::optional<std::pair<double, double>>
    project_world(const Camera& cam, const util::Vec3d& world);

private:
    mutable std::mutex                              mu_;
    std::unordered_map<std::uint64_t, Annotation>   items_;
    std::uint64_t                                   nextId_ = 1;
};

}  // namespace simall::visualization
