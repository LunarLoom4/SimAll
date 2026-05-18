// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/PythonEngine.hpp
//
// `IScriptEngine` implementation backed by the embedded Python interpreter.
// When Python is not built in, this engine still installs cleanly but every
// call short-circuits with an explanatory message — preserving the
// `simall::scripting::engine()` contract for the GUI / UDF code.
// =============================================================================
#pragma once

#include "scripting/PyBindings.hpp"
#include "scripting/ScriptEngine.hpp"

namespace simall::scripting {

class PythonEngine : public IScriptEngine {
public:
    PythonEngine();
    ~PythonEngine() override;
    void        execute(const std::string& src) override;
    std::string repl   (const std::string& line) override;

    // Install this engine as the global `engine()` instance.  Returns the
    // previous engine pointer so callers can restore it.
    static IScriptEngine* install();
};

}  // namespace simall::scripting
