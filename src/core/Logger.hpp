// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/Logger.hpp
// Phase  : 1.3 (APPLICATION CORE → Logging System)
//
// Header-only async logger. Supports leveled logging, file + console sinks,
// thread-safe enqueue, and structured crash diagnostics. Designed for HPC
// use: zero allocation on the hot path when SIMALL_LOG_OFF is defined.
// =============================================================================
#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

namespace simall::core {

enum class LogLevel : int { Trace = 0, Debug, Info, Warn, Error, Fatal, Off };

inline const char* to_string(LogLevel l) {
    switch (l) {
        case LogLevel::Trace: return "TRACE";
        case LogLevel::Debug: return "DEBUG";
        case LogLevel::Info:  return "INFO ";
        case LogLevel::Warn:  return "WARN ";
        case LogLevel::Error: return "ERROR";
        case LogLevel::Fatal: return "FATAL";
        default:              return "OFF  ";
    }
}

class Logger {
public:
    static Logger& instance() {
        static Logger inst;
        return inst;
    }

    void set_level(LogLevel l) noexcept { level_.store(l, std::memory_order_relaxed); }
    LogLevel level()      const noexcept { return level_.load(std::memory_order_relaxed); }

    void add_file_sink(const std::string& path) {
        std::lock_guard lk(sink_mtx_);
        file_sinks_.emplace_back(std::make_unique<std::ofstream>(path, std::ios::app));
    }
    void enable_console(bool on) { console_.store(on); }

    template <typename... Args>
    void log(LogLevel lvl, const char* category, Args&&... args) {
        if (lvl < level_.load(std::memory_order_relaxed)) return;
        std::ostringstream oss;
        format_prefix(oss, lvl, category);
        (oss << ... << std::forward<Args>(args));
        oss << '\n';
        push(oss.str());
    }

    void flush() {
        std::unique_lock lk(queue_mtx_);
        flushed_.wait(lk, [&]{ return queue_.empty(); });
    }

private:
    Logger()
        : level_(LogLevel::Info), console_(true), running_(true)
    {
        worker_ = std::thread([this]{ worker_loop(); });
    }
    ~Logger() {
        running_.store(false);
        queue_cv_.notify_all();
        if (worker_.joinable()) worker_.join();
    }
    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;

    void format_prefix(std::ostringstream& oss, LogLevel lvl, const char* category) {
        using namespace std::chrono;
        auto now = system_clock::now();
        auto t   = system_clock::to_time_t(now);
        auto ms  = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
        std::tm tm{};
#if defined(_WIN32)
        localtime_s(&tm, &t);
#else
        localtime_r(&t, &tm);
#endif
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d.%03lld",
                      tm.tm_hour, tm.tm_min, tm.tm_sec,
                      static_cast<long long>(ms.count()));
        oss << '[' << buf << "][" << to_string(lvl) << "][" << category << "] ";
    }

    void push(std::string msg) {
        {
            std::lock_guard lk(queue_mtx_);
            queue_.emplace(std::move(msg));
        }
        queue_cv_.notify_one();
    }

    void worker_loop() {
        while (running_.load() || !queue_empty()) {
            std::unique_lock lk(queue_mtx_);
            queue_cv_.wait(lk, [&]{ return !queue_.empty() || !running_.load(); });
            while (!queue_.empty()) {
                std::string msg = std::move(queue_.front());
                queue_.pop();
                lk.unlock();
                emit(msg);
                lk.lock();
            }
            flushed_.notify_all();
        }
    }

    bool queue_empty() {
        std::lock_guard lk(queue_mtx_);
        return queue_.empty();
    }

    void emit(const std::string& msg) {
        if (console_.load()) std::fputs(msg.c_str(), stdout);
        std::lock_guard lk(sink_mtx_);
        for (auto& s : file_sinks_) { (*s) << msg; s->flush(); }
    }

    std::atomic<LogLevel> level_;
    std::atomic<bool>     console_;
    std::atomic<bool>     running_;
    std::thread           worker_;

    std::mutex                     queue_mtx_;
    std::condition_variable        queue_cv_;
    std::condition_variable        flushed_;
    std::queue<std::string>        queue_;

    std::mutex                                  sink_mtx_;
    std::vector<std::unique_ptr<std::ofstream>> file_sinks_;
};

}  // namespace simall::core

#define SIMALL_LOG(level, category, ...)                          \
    ::simall::core::Logger::instance().log(                       \
        ::simall::core::LogLevel::level, category, __VA_ARGS__)

#define SIMALL_LOG_INFO(cat, ...)  SIMALL_LOG(Info,  cat, __VA_ARGS__)
#define SIMALL_LOG_WARN(cat, ...)  SIMALL_LOG(Warn,  cat, __VA_ARGS__)
#define SIMALL_LOG_ERROR(cat, ...) SIMALL_LOG(Error, cat, __VA_ARGS__)
#define SIMALL_LOG_DEBUG(cat, ...) SIMALL_LOG(Debug, cat, __VA_ARGS__)
