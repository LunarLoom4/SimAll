// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/CrashHandler.hpp
// Phase  : 1.3 (APPLICATION CORE → Crash Diagnostics)
//
// Cross-platform crash handler. Installs signal handlers (POSIX: SIGSEGV,
// SIGABRT, SIGFPE, SIGILL, SIGBUS) and an unhandled-exception filter on
// Windows that writes a minidump via MiniDumpWriteDump (dbghelp). The
// destination directory and a user "before-crash" callback may be
// configured. The handler is best-effort: it logs the failure, writes a
// minidump or stack trace where supported, and re-raises to abort the
// process.
// =============================================================================
#pragma once

#include <filesystem>
#include <functional>
#include <string>

namespace simall::core
{

struct CrashHandlerConfig
{
    std::filesystem::path dumpDirectory; // empty → cwd
    std::string appName = "SimAll";
    std::function<void(const std::string& reason)> onBeforeDump;
};

class CrashHandler
{
public:
    /// Install once at program start. Idempotent; subsequent calls update
    /// configuration. Returns false if installation failed (rare).
    static bool install(CrashHandlerConfig cfg = {});

    /// Manually trigger a dump (useful for diagnostics on user request).
    /// Returns the path of the dump file on success.
    static std::filesystem::path writeDumpNow(std::string_view reason);

    /// Uninstall (mainly for tests).
    static void uninstall();
};

} // namespace simall::core
