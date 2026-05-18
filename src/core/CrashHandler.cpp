// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/CrashHandler.cpp
// =============================================================================
#include "core/CrashHandler.hpp"
#include "core/Logger.hpp"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <mutex>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN
#  include <Windows.h>
#  include <DbgHelp.h>
#  pragma comment(lib, "dbghelp.lib")
#else
#  include <csignal>
#  include <unistd.h>
#  if __has_include(<execinfo.h>)
#    include <execinfo.h>
#    define SIMALL_HAVE_BACKTRACE 1
#  endif
#endif

namespace simall::core {

namespace {

std::mutex                gMtx;
CrashHandlerConfig        gCfg;
std::atomic<bool>         gInstalled{false};

std::string timestamp() {
    auto now = std::chrono::system_clock::now();
    auto t   = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%04d%02d%02d-%02d%02d%02d",
                  tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                  tm.tm_hour, tm.tm_min, tm.tm_sec);
    return buf;
}

std::filesystem::path makeDumpPath(std::string_view tag) {
    std::filesystem::path dir = gCfg.dumpDirectory;
    if (dir.empty()) {
        std::error_code ec;
        dir = std::filesystem::current_path(ec);
    }
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    std::string fname = gCfg.appName + "-" + std::string(tag) + "-" + timestamp();
#if defined(_WIN32)
    fname += ".dmp";
#else
    fname += ".log";
#endif
    return dir / fname;
}

#if defined(_WIN32)

LONG WINAPI seh_filter(EXCEPTION_POINTERS* ep) {
    std::lock_guard lk(gMtx);
    auto path = makeDumpPath("crash");
    if (gCfg.onBeforeDump) try { gCfg.onBeforeDump(path.string()); } catch (...) {}

    HANDLE hFile = ::CreateFileW(path.wstring().c_str(), GENERIC_WRITE, 0, nullptr,
                                 CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (hFile != INVALID_HANDLE_VALUE) {
        MINIDUMP_EXCEPTION_INFORMATION info{};
        info.ThreadId          = ::GetCurrentThreadId();
        info.ExceptionPointers = ep;
        info.ClientPointers    = FALSE;
        MINIDUMP_TYPE type = static_cast<MINIDUMP_TYPE>(
            MiniDumpWithDataSegs | MiniDumpWithHandleData |
            MiniDumpWithThreadInfo | MiniDumpWithUnloadedModules);
        ::MiniDumpWriteDump(::GetCurrentProcess(), ::GetCurrentProcessId(),
                            hFile, type, &info, nullptr, nullptr);
        ::CloseHandle(hFile);
    }
    SIMALL_LOG_ERROR("crash", "Unhandled exception. Minidump: ", path.string());
    Logger::instance().flush();
    return EXCEPTION_EXECUTE_HANDLER;
}

#else  // POSIX

void posix_signal(int sig) {
    // Async-signal-safe context: minimal work.
    std::lock_guard lk(gMtx);
    auto path = makeDumpPath("crash");
    if (gCfg.onBeforeDump) try { gCfg.onBeforeDump(path.string()); } catch (...) {}
    FILE* f = std::fopen(path.string().c_str(), "w");
    if (f) {
        std::fprintf(f, "SimAll crash: signal %d\n", sig);
#ifdef SIMALL_HAVE_BACKTRACE
        void* frames[64];
        int n = ::backtrace(frames, 64);
        ::backtrace_symbols_fd(frames, n, fileno(f));
#endif
        std::fclose(f);
    }
    SIMALL_LOG_ERROR("crash", "Caught signal ", sig, " — wrote ", path.string());
    Logger::instance().flush();
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

#endif

}  // namespace

bool CrashHandler::install(CrashHandlerConfig cfg) {
    std::lock_guard lk(gMtx);
    gCfg = std::move(cfg);
    if (gInstalled.exchange(true)) return true;

#if defined(_WIN32)
    ::SetUnhandledExceptionFilter(&seh_filter);
#else
    struct sigaction sa{};
    sa.sa_handler = &posix_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESETHAND;
    for (int s : {SIGSEGV, SIGABRT, SIGFPE, SIGILL, SIGBUS}) {
        ::sigaction(s, &sa, nullptr);
    }
#endif
    return true;
}

std::filesystem::path CrashHandler::writeDumpNow(std::string_view reason) {
    std::lock_guard lk(gMtx);
    auto path = makeDumpPath("manual");
#if defined(_WIN32)
    HANDLE hFile = ::CreateFileW(path.wstring().c_str(), GENERIC_WRITE, 0, nullptr,
                                 CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, nullptr);
    if (hFile != INVALID_HANDLE_VALUE) {
        ::MiniDumpWriteDump(::GetCurrentProcess(), ::GetCurrentProcessId(),
                            hFile, MiniDumpWithThreadInfo, nullptr, nullptr, nullptr);
        ::CloseHandle(hFile);
    }
#else
    FILE* f = std::fopen(path.string().c_str(), "w");
    if (f) {
        std::fprintf(f, "Manual dump: %.*s\n", int(reason.size()), reason.data());
#ifdef SIMALL_HAVE_BACKTRACE
        void* frames[64];
        int n = ::backtrace(frames, 64);
        ::backtrace_symbols_fd(frames, n, fileno(f));
#endif
        std::fclose(f);
    }
#endif
    SIMALL_LOG_INFO("crash", "Manual dump (", reason, ") → ", path.string());
    return path;
}

void CrashHandler::uninstall() {
    std::lock_guard lk(gMtx);
    if (!gInstalled.exchange(false)) return;
#if defined(_WIN32)
    ::SetUnhandledExceptionFilter(nullptr);
#else
    for (int s : {SIGSEGV, SIGABRT, SIGFPE, SIGILL, SIGBUS}) {
        std::signal(s, SIG_DFL);
    }
#endif
}

}  // namespace simall::core
