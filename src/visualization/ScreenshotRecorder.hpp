// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/ScreenshotRecorder.hpp
// Phase  : 3 / Week 14 (recording: framebuffer → file)
//
// Persists Image objects to disk.  Supports binary PPM (P6) which has no
// external dependencies and round-trips losslessly for verification.  A raw
// RGBA blob writer is included for pipelines that want to hand the bytes
// to ffmpeg / vapoursynth in a streaming fashion.
// =============================================================================
#pragma once

#include "visualization/VisualizationTypes.hpp"

#include <filesystem>
#include <string>

namespace simall::visualization {

class ScreenshotRecorder {
public:
    /// Write an Image to a binary PPM (P6) file.  Alpha is dropped because
    /// PPM is RGB-only; callers wanting alpha should use save_rgba_raw().
    /// Returns false on I/O failure.
    static bool save_ppm(const Image& image, const std::filesystem::path& path);

    /// Write the raw RGBA bytes (width*height*4) preceded by an 8-byte
    /// little-endian header [w:uint32, h:uint32].  Compatible with the
    /// SimAll thumbnail cache.
    static bool save_rgba_raw(const Image& image, const std::filesystem::path& path);

    /// Read back a PPM (P6) file into an Image.  Returns an empty Image on
    /// failure.  Used by tests for round-trip verification.
    static Image load_ppm(const std::filesystem::path& path);
};

}  // namespace simall::visualization
