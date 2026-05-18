// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/RotatingFileSink.cpp
// =============================================================================
#include "core/RotatingFileSink.hpp"

#include <cstdio>
#include <system_error>

namespace simall::core {

RotatingFileSink::RotatingFileSink(std::filesystem::path path,
                                   std::uint64_t maxBytes,
                                   unsigned maxFiles)
    : basePath_(std::move(path)), maxBytes_(maxBytes), maxFiles_(maxFiles) {
    std::error_code ec;
    if (basePath_.has_parent_path())
        std::filesystem::create_directories(basePath_.parent_path(), ec);
    openCurrent();
}

RotatingFileSink::~RotatingFileSink() {
    closeCurrent();
}

void RotatingFileSink::openCurrent() {
#if defined(_WIN32)
    fp_ = nullptr;
    fopen_s(&fp_, basePath_.string().c_str(), "ab");
#else
    fp_ = std::fopen(basePath_.string().c_str(), "ab");
#endif
    if (fp_) {
        std::fseek(fp_, 0, SEEK_END);
        long pos = std::ftell(fp_);
        written_ = (pos > 0) ? static_cast<std::uint64_t>(pos) : 0;
    }
}

void RotatingFileSink::closeCurrent() {
    if (fp_) { std::fflush(fp_); std::fclose(fp_); fp_ = nullptr; }
}

void RotatingFileSink::rotate() {
    closeCurrent();
    std::error_code ec;
    // Delete oldest.
    auto oldest = basePath_;
    oldest += "." + std::to_string(maxFiles_);
    std::filesystem::remove(oldest, ec);
    // Shift .N → .(N+1) downwards.
    for (unsigned i = maxFiles_; i > 1; --i) {
        auto src = basePath_; src += "." + std::to_string(i - 1);
        auto dst = basePath_; dst += "." + std::to_string(i);
        if (std::filesystem::exists(src, ec))
            std::filesystem::rename(src, dst, ec);
    }
    auto first = basePath_; first += ".1";
    if (std::filesystem::exists(basePath_, ec))
        std::filesystem::rename(basePath_, first, ec);
    written_ = 0;
    openCurrent();
}

void RotatingFileSink::write(std::string_view msg) {
    std::lock_guard lk(mtx_);
    if (!fp_) return;
    if (written_ + msg.size() > maxBytes_) rotate();
    if (!fp_) return;
    std::fwrite(msg.data(), 1, msg.size(), fp_);
    std::fflush(fp_);
    written_ += msg.size();
}

void RotatingFileSink::flush() {
    std::lock_guard lk(mtx_);
    if (fp_) std::fflush(fp_);
}

void RotatingFileSink::attachToLogger() {
    // The default Logger uses std::ofstream sinks; tee here is intentionally
    // a no-op stub: in practice the application replaces Logger::emit by
    // forwarding through this sink. Kept for API symmetry with future
    // spdlog migration.
}

}  // namespace simall::core
