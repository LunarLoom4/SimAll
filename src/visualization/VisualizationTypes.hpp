// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/VisualizationTypes.hpp
// Phase  : 3 / Week 14 (post-processing + viewport data model)
//
// Shared, VTK-free data structures used by every Week-14 module.  Filters,
// the scene graph, the picker, and the recorders all speak in terms of
// these plain-old-data containers so they can be unit-tested without a
// display server and so the actual VTK adapter inside Viewport stays a
// thin translation layer.
// =============================================================================
#pragma once

#include "utilities/MathTypes.hpp"

#include <array>
#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace simall::visualization {

using ActorId  = std::uint64_t;
using NodeId   = std::uint64_t;
using FrameIdx = std::uint64_t;

inline constexpr ActorId kInvalidActorId = 0;
inline constexpr NodeId  kInvalidNodeId  = 0;

/// RGBA in linear floating-point (renderer applies gamma at the end).
using Color4 = std::array<float, 4>;

// ---------------------------------------------------------------------------
//  Geometry views (POD) — every filter/renderer hands one of these around.
// ---------------------------------------------------------------------------

/// Indexed triangle surface. Triangles are stored as flat triples in
/// triIndex (size = 3*numTriangles). pointScalars / triScalars are optional
/// — empty when the consumer is colour-blind to the field.
struct SurfaceMesh {
    std::vector<util::Vec3d> points;
    std::vector<std::int32_t> triIndex;        // 3*nTri
    std::vector<double> pointScalars;          // size == 0 OR points.size()
    std::vector<double> triScalars;            // size == 0 OR triIndex.size()/3

    std::size_t triangle_count() const noexcept { return triIndex.size() / 3; }
    bool empty() const noexcept { return triIndex.empty(); }
};

/// Indexed tetrahedral volume. Tets stored as flat quads (size = 4*nTet).
/// pointScalars[i] is the value at points[i]; required by IsoSurface and
/// SectionCut. cellVectors[i] is an optional flow vector per tet used by
/// StreamlineRk4 in "cell-data" sampling mode.
struct VolumeMesh {
    std::vector<util::Vec3d>  points;
    std::vector<std::int32_t> tetIndex;        // 4*nTet
    std::vector<double>       pointScalars;
    std::vector<util::Vec3d>  cellVectors;     // size == 0 OR nTet

    std::size_t tet_count() const noexcept { return tetIndex.size() / 4; }
    bool empty() const noexcept { return tetIndex.empty(); }
};

/// Polyline soup in CSR layout (lineOffsets has size nLines+1).
struct LineSet {
    std::vector<util::Vec3d>  points;
    std::vector<std::int32_t> lineOffsets;     // CSR; offsets into points
    std::vector<double>       scalars;         // size 0 OR points.size()

    std::size_t line_count() const noexcept {
        return lineOffsets.empty() ? 0u : lineOffsets.size() - 1u;
    }
};

/// Oriented arrow primitives. The renderer draws each one as a
/// scale·direction shaft anchored at `anchor`.
struct GlyphSet {
    std::vector<util::Vec3d> anchor;
    std::vector<util::Vec3d> direction;        // unit if a glyph is meant to be drawn
    std::vector<double>      magnitude;        // |raw_vector| at sample
};

/// RGBA8 framebuffer used by ScreenshotRecorder and the volume / LIC tests.
/// pixels[(y*width + x)*4 + c], rows top-to-bottom.
struct Image {
    std::uint32_t          width  = 0;
    std::uint32_t          height = 0;
    std::vector<std::uint8_t> pixels;          // size = 4*w*h
    void resize(std::uint32_t w, std::uint32_t h) {
        width = w; height = h;
        pixels.assign(static_cast<std::size_t>(4) * w * h, 0);
    }
    bool empty() const noexcept { return width == 0 || height == 0; }
};

// ---------------------------------------------------------------------------
//  Transfer function / colormap
// ---------------------------------------------------------------------------

/// Piecewise-linear scalar→RGBA mapping defined by sorted control points.
struct TransferFunction {
    struct Stop {
        double scalar;
        Color4 color;
    };
    std::vector<Stop> stops;

    /// Sample the gradient at `s`. Returns the first stop colour below the
    /// range and the last stop colour above it (clamped).  Empty TF returns
    /// fully opaque white as a safe default.
    Color4 sample(double s) const {
        if (stops.empty()) return {1, 1, 1, 1};
        if (s <= stops.front().scalar) return stops.front().color;
        if (s >= stops.back().scalar)  return stops.back().color;
        for (std::size_t i = 1; i < stops.size(); ++i) {
            if (s <= stops[i].scalar) {
                const double t =
                    (s - stops[i-1].scalar) /
                    (stops[i].scalar - stops[i-1].scalar);
                Color4 out{};
                for (int c = 0; c < 4; ++c) {
                    out[c] = static_cast<float>(
                        (1.0 - t) * stops[i-1].color[c] + t * stops[i].color[c]);
                }
                return out;
            }
        }
        return stops.back().color;
    }
};

/// Standard "cool→warm" diverging colormap (Kenneth Moreland).  Useful when
/// callers do not bring their own transfer function.
inline TransferFunction make_cool_warm(double low, double high) {
    return TransferFunction{{
        {low,                          {0.230f, 0.299f, 0.754f, 1.0f}},
        {0.5 * (low + high),           {0.865f, 0.865f, 0.865f, 1.0f}},
        {high,                         {0.706f, 0.016f, 0.150f, 1.0f}}
    }};
}

// ---------------------------------------------------------------------------
//  Sampler abstraction (used by StreamlineRk4 / VectorGlyph / LIC)
// ---------------------------------------------------------------------------

/// User-supplied analytic-or-tabulated vector sampler.  Returns nullopt when
/// the point lies outside the field's domain — this terminates the streamline.
using VectorSampler = std::function<std::optional<util::Vec3d>(const util::Vec3d&)>;

/// Scalar counterpart used by IsoSurface (when callers prefer an analytic
/// field) and by Annotation/probe operations.
using ScalarSampler = std::function<std::optional<double>(const util::Vec3d&)>;

/// Camera state shared by CameraController and PickingBridge.  Mirrors
/// Viewport::CameraState but lives in the data-pure layer so picking can be
/// tested without QVTKOpenGLNativeWidget.
struct Camera {
    util::Vec3d position    {0, 0, 10};
    util::Vec3d focalPoint  {0, 0,  0};
    util::Vec3d up          {0, 1,  0};
    double      fovYDegrees  = 30.0;
    double      nearPlane    = 0.1;
    double      farPlane     = 1000.0;
    bool        orthographic = false;
    double      orthoHeight  = 2.0;   // world units; only used when orthographic
    std::uint32_t viewportWidth  = 800;
    std::uint32_t viewportHeight = 600;

    double aspect() const noexcept {
        return viewportHeight == 0 ? 1.0
            : static_cast<double>(viewportWidth) /
              static_cast<double>(viewportHeight);
    }
};

}  // namespace simall::visualization
