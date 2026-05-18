// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ScalarBar.hpp
// Phase  : 3 / Week 14 (overlays: colour-bar legend)
//
// Pure data + a CPU rasteriser that emits an Image of the colour bar.  The
// rasteriser exists so headless verification can confirm the lookup-table /
// tick layout without instantiating Qt.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <string>
#include <vector>

namespace simall::visualization {

enum class ScalarBarOrientation : std::uint8_t { Vertical, Horizontal };

struct ScalarBar {
    std::string          title;
    std::string          units;
    double               scalarMin = 0.0;
    double               scalarMax = 1.0;
    std::uint32_t        tickCount = 5;
    int                  precision = 3;
    ScalarBarOrientation orientation = ScalarBarOrientation::Vertical;
    std::uint32_t        widthPx   = 32;
    std::uint32_t        lengthPx  = 256;
    TransferFunction     tf;

    /// Sample tick positions and the printable text for each.
    struct Tick { double scalar; std::string label; };
    std::vector<Tick> ticks() const;

    /// Rasterise the gradient ramp at the configured resolution.  Returns
    /// an Image (gradient only, no labels — labels are painted by the GUI).
    Image rasterise_gradient() const;
};

}  // namespace simall::visualization
