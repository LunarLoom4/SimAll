// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/AnimationRecorder.cpp
// =============================================================================
#include "visualization/AnimationRecorder.hpp"

#include <iomanip>
#include <sstream>
#include <system_error>

namespace simall::visualization
{

AnimationRecorder::AnimationRecorder(AnimationConfig cfg) : cfg_(std::move(cfg))
{
    if (!cfg_.outputDir.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(cfg_.outputDir, ec);
    }
    nextDueTime_ = cfg_.startTime;
}

bool AnimationRecorder::submit(double simTime, Image image)
{
    std::lock_guard lock(mu_);
    if (simTime < cfg_.startTime || simTime > cfg_.endTime) {
        ++stats_.framesDropped;
        return false;
    }
    if (simTime + 1e-12 < nextDueTime_) {
        ++stats_.framesDropped;
        stats_.lastSimTime = simTime;
        return false;
    }
    queue_.push_back(Pending{simTime, std::move(image), stats_.framesQueued});
    ++stats_.framesQueued;
    stats_.lastSimTime = simTime;
    // Advance schedule.
    const double dt = (cfg_.fps == 0) ? 0.0 : 1.0 / static_cast<double>(cfg_.fps);
    nextDueTime_ = simTime + dt;
    return true;
}

FrameIdx AnimationRecorder::flush()
{
    std::lock_guard lock(mu_);
    FrameIdx written = 0;
    for (auto& p : queue_) {
        const auto path = frame_path(p.index);
        const bool ok = cfg_.rgbaRaw ? ScreenshotRecorder::save_rgba_raw(p.img, path)
                                     : ScreenshotRecorder::save_ppm(p.img, path);
        if (ok) {
            ++written;
            ++stats_.framesWritten;
        }
    }
    queue_.clear();
    return written;
}

AnimationStats AnimationRecorder::stats() const
{
    std::lock_guard lock(mu_);
    return stats_;
}

void AnimationRecorder::reset()
{
    std::lock_guard lock(mu_);
    stats_ = {};
    queue_.clear();
    nextDueTime_ = cfg_.startTime;
}

std::filesystem::path AnimationRecorder::frame_path(FrameIdx idx) const
{
    std::ostringstream oss;
    oss << cfg_.filenameStem << '_' << std::setfill('0') << std::setw(cfg_.filenameDigits) << idx
        << (cfg_.rgbaRaw ? ".rgba" : ".ppm");
    return cfg_.outputDir / oss.str();
}

} // namespace simall::visualization
