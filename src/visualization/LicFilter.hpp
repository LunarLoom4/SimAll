// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/LicFilter.hpp
// Phase  : 3 / Week 14 (filters: 2D Line Integral Convolution)
//
// Cabral & Leedom 1993 Line Integral Convolution: per output pixel,
// integrate a streamline forward and backward through a 2D vector field,
// accumulate a per-pixel white-noise texture along that path, normalise.
// Used for dense vector visualisation on planar sections / 2D slices.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <vector>

namespace simall::visualization {

struct RegularGrid2 {
    std::uint32_t      nx = 0, ny = 0;
    double             originX = 0.0, originY = 0.0;
    double             spacingX = 1.0, spacingY = 1.0;
    std::vector<float> vx;       // size nx*ny
    std::vector<float> vy;

    bool valid() const noexcept {
        return nx > 0 && ny > 0 &&
               vx.size() == static_cast<std::size_t>(nx) * ny &&
               vy.size() == vx.size();
    }
};

struct LicConfig {
    std::uint32_t       width    = 0;     // 0 → use grid resolution
    std::uint32_t       height   = 0;
    std::uint32_t       stepsPerSide = 24;
    double              stepSize = 1.0;   // in grid cells
    std::uint64_t       noiseSeed = 0xC0FFEEUL;
};

class LicFilter {
public:
    /// Generate the LIC image.  pixels are grayscale-on-RGBA8 (R=G=B=intensity,
    /// A=255).
    static Image generate(const RegularGrid2& field, const LicConfig& cfg = {});
};

}  // namespace simall::visualization
