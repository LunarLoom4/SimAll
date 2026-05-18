// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/PyBindings.hpp
//
// Public façade for the pybind11-based Python embedding.  Compiles to a
// no-op when SIMALL_HAVE_PYTHON is not defined so the rest of the codebase
// can call `python::start()` / `stop()` / `exec()` / `eval()`
// unconditionally.  When Python is built in, the embedded interpreter is
// initialised on `start()`, the `simall` Python module is registered, and
// the shared `ScriptContext` / `CommandBus` / `UdfHost` singletons are
// exposed.
// =============================================================================
#pragma once

#include <functional>
#include <string>

namespace simall::scripting::python {

[[nodiscard]] bool   is_built_in() noexcept;   // true iff SIMALL_HAVE_PYTHON
[[nodiscard]] bool   is_running()  noexcept;

void                 start();                  // idempotent
void                 stop();                   // idempotent; safe at shutdown

// Execute a code block.  Returns "" on success or the formatted error.
std::string          exec(const std::string& code);

// Evaluate a Python expression and return its repr().  Empty string on error
// (error text is written to the install_output callback if set).
std::string          eval(const std::string& expression);

// Install an output sink — receives stdout / stderr writes from the
// embedded interpreter.  Called from any thread; consumer must serialise.
using OutputSink = std::function<void(const std::string&)>;
void                 install_output(OutputSink sink);

// Read & clear the buffered captured stdout (no-op when no sink set).
std::string          drain_captured_output();

}  // namespace simall::scripting::python
