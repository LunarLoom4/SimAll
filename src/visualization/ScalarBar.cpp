// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ScalarBar.cpp
// =============================================================================
#include "visualization/ScalarBar.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace simall::visualization
{

std::vector<ScalarBar::Tick> ScalarBar::ticks() const
{
    std::vector<Tick> out;
    if (tickCount == 0 || scalarMax <= scalarMin)
        return out;
    out.reserve(tickCount);
    for (std::uint32_t i = 0; i < tickCount; ++i) {
        const double t =
            (tickCount == 1) ? 0.5 : static_cast<double>(i) / static_cast<double>(tickCount - 1);
        const double s = scalarMin + t * (scalarMax - scalarMin);
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(precision) << s;
        if (!units.empty())
            oss << ' ' << units;
        out.push_back({s, oss.str()});
    }
    return out;
}

Image ScalarBar::rasterise_gradient() const
{
    Image img;
    if (lengthPx == 0 || widthPx == 0)
        return img;
    const std::uint32_t W = (orientation == ScalarBarOrientation::Vertical) ? widthPx : lengthPx;
    const std::uint32_t H = (orientation == ScalarBarOrientation::Vertical) ? lengthPx : widthPx;
    img.resize(W, H);

    for (std::uint32_t y = 0; y < H; ++y) {
        for (std::uint32_t x = 0; x < W; ++x) {
            double t = 0.0;
            if (orientation == ScalarBarOrientation::Vertical) {
                t = 1.0
                    - static_cast<double>(y)
                          / static_cast<double>(std::max<std::uint32_t>(H - 1, 1));
            } else {
                t = static_cast<double>(x) / static_cast<double>(std::max<std::uint32_t>(W - 1, 1));
            }
            const double s = scalarMin + t * (scalarMax - scalarMin);
            const Color4 c = tf.sample(s);
            const std::size_t off = (static_cast<std::size_t>(y) * W + x) * 4u;
            img.pixels[off + 0] =
                static_cast<std::uint8_t>(std::clamp(c[0] * 255.0f, 0.0f, 255.0f));
            img.pixels[off + 1] =
                static_cast<std::uint8_t>(std::clamp(c[1] * 255.0f, 0.0f, 255.0f));
            img.pixels[off + 2] =
                static_cast<std::uint8_t>(std::clamp(c[2] * 255.0f, 0.0f, 255.0f));
            img.pixels[off + 3] =
                static_cast<std::uint8_t>(std::clamp(c[3] * 255.0f, 0.0f, 255.0f));
        }
    }
    return img;
}

} // namespace simall::visualization
