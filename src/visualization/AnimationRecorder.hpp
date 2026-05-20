// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/AnimationRecorder.hpp
// Phase  : 3 / Week 14 (recording: time-stepped frame queue)
//
// Drives ScreenshotRecorder at fixed intervals.  The renderer pushes
// (simTime, image) pairs into the recorder; the recorder decides whether
// each frame should be persisted given the configured frames-per-second
// schedule and output directory.  No background threads — flush() drains
// synchronously so unit tests can verify file contents deterministically.
// =============================================================================
#pragma once

#include "visualization/ScreenshotRecorder.hpp"
#include "visualization/VisualizationTypes.hpp"

#include <filesystem>
#include <mutex>
#include <string>
#include <vector>

namespace simall::visualization
{

struct AnimationConfig
{
    std::filesystem::path outputDir;
    std::string filenameStem = "frame";
    std::uint32_t fps = 30; // wall-clock framerate
    double startTime = 0.0;
    double endTime = 1e300;
    bool rgbaRaw = false; // false → PPM, true → RGBA blob
    std::uint32_t filenameDigits = 6;
};

struct AnimationStats
{
    FrameIdx framesQueued = 0;
    FrameIdx framesWritten = 0;
    FrameIdx framesDropped = 0; // frames arriving before next due
    double lastSimTime = -1.0;
};

class AnimationRecorder
{
public:
    explicit AnimationRecorder(AnimationConfig cfg);

    /// Submit a frame for the given simulation time.  Returns true when the
    /// frame was kept (not necessarily yet written to disk).
    bool submit(double simTime, Image image);

    /// Drain pending frames to disk in submission order.  Returns the
    /// number of new files created.
    FrameIdx flush();

    AnimationStats stats() const;
    void reset();

private:
    struct Pending
    {
        double t;
        Image img;
        FrameIdx index;
    };

    std::filesystem::path frame_path(FrameIdx idx) const;

    mutable std::mutex mu_;
    AnimationConfig cfg_;
    AnimationStats stats_;
    std::vector<Pending> queue_;
    double nextDueTime_ = 0.0;
};

} // namespace simall::visualization
