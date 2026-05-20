// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/PythonEngine.cpp
// =============================================================================
#include "scripting/PythonEngine.hpp"

#include <atomic>

namespace simall::scripting
{

namespace
{
std::atomic<IScriptEngine*> g_installed{nullptr};
}

PythonEngine::PythonEngine()
{
    if (python::is_built_in())
        python::start();
}

PythonEngine::~PythonEngine()
{
    if (g_installed.load() == this)
        g_installed.store(nullptr);
    // Note: we do NOT stop the interpreter here; ownership of the interpreter
    // lifetime is in PyBindings.cpp (one process, one interpreter).
}

void PythonEngine::execute(const std::string& src)
{
    auto err = python::exec(src);
    (void) err; // surfaced via the output sink; caller can also drain
}

std::string PythonEngine::repl(const std::string& line)
{
    if (!python::is_built_in())
        return "Python scripting not built in.";
    // First try to eval (expression); fall back to exec (statement) so the
    // REPL prints results for "1+1" but accepts "x = 1".
    if (line.find_first_of("=:;") == std::string::npos) {
        auto r = python::eval(line);
        if (!r.empty())
            return r;
    }
    auto err = python::exec(line);
    if (!err.empty())
        return err;
    return python::drain_captured_output();
}

IScriptEngine* PythonEngine::install()
{
    auto* prev = g_installed.exchange(new PythonEngine());
    // The `engine()` accessor in ScriptEngine.cpp returns a static NullEngine
    // by default; downstream code that wants the real Python engine should
    // call PythonEngine::install() during bootstrap and use the returned
    // pointer (or rely on a dedicated accessor wrapped around g_installed).
    return prev;
}

} // namespace simall::scripting
