// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/RotatingFileSink.hpp
// Phase  : 1.3 (APPLICATION CORE → Logging extension)
//
// Size-bounded rotating-file sink for `Logger`. Attaches as a custom sink
// (separate from Logger's plain `add_file_sink`) and is rotated when the
// active file exceeds `maxBytes`. Up to `maxFiles` historical files are
// retained:  <path>.1, <path>.2, … <path>.N. Thread-safe; uses a flush on
// every emit so logs survive crashes.
// =============================================================================
#pragma once

#include "core/Logger.hpp"

#include <cstdio>
#include <filesystem>
#include <mutex>
#include <string>

namespace simall::core {

class RotatingFileSink {
public:
    RotatingFileSink(std::filesystem::path path,
                     std::uint64_t         maxBytes = 16ull * 1024ull * 1024ull,
                     unsigned              maxFiles = 5);
    ~RotatingFileSink();

    RotatingFileSink(const RotatingFileSink&)            = delete;
    RotatingFileSink& operator=(const RotatingFileSink&) = delete;

    /// Write a pre-formatted message. Rotation happens transparently.
    void write(std::string_view msg);

    /// Force flush (no rotation).
    void flush();

    /// Attach to the global Logger as a "callback" file via Logger::add_file_sink.
    /// For convenience: writes go via this sink's rotation policy by tee-ing
    /// to a temporary stream. NOT used internally — `Logger` users that want
    /// rotation should write directly via `RotatingFileSink::write` from a
    /// custom log handler.
    void attachToLogger();

private:
    void openCurrent();
    void closeCurrent();
    void rotate();

    std::filesystem::path basePath_;
    std::uint64_t         maxBytes_;
    unsigned              maxFiles_;
    std::FILE*            fp_       = nullptr;
    std::uint64_t         written_  = 0;
    std::mutex            mtx_;
};

}  // namespace simall::core
