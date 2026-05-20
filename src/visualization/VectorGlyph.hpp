// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/VectorGlyph.hpp
// Phase  : 3 / Week 14 (filters: vector glyph emission)
//
// Builds an arrow GlyphSet from a vector sampler at a list of seed points.
// Magnitude is captured separately so the renderer can apply linear or log
// scale modes.  Optional culling drops glyphs below a magnitude threshold
// to declutter low-velocity regions.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

namespace simall::visualization
{

struct GlyphConfig
{
    double minMagnitude = 0.0; // skip glyphs with |v| <= this
    bool normalizeDirection = true;
    double maxAnchors = 1.0e9; // hard upper bound to protect the GUI
};

class VectorGlyph
{
public:
    static GlyphSet sample(const VectorSampler& sampler,
                           const std::vector<util::Vec3d>& seeds,
                           const GlyphConfig& cfg = {});

    /// Emit glyphs on a regular 3D lattice covering `box` with `nx*ny*nz`
    /// sample points.  Useful for quick whole-domain views.
    static GlyphSet sample_lattice(const VectorSampler& sampler,
                                   const util::BoundingBox& box,
                                   std::uint32_t nx,
                                   std::uint32_t ny,
                                   std::uint32_t nz,
                                   const GlyphConfig& cfg = {});
};

} // namespace simall::visualization
