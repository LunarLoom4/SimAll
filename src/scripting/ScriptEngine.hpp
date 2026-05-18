// Phase 20 (SCRIPTING). Python embedding via pybind11 is wired in Phase 20;
// here we provide a string-evaluation contract that the Python console
// and UDF system both target.
#pragma once
#include <functional>
#include <string>
namespace simall::scripting {
class IScriptEngine {
public:
    virtual ~IScriptEngine() = default;
    virtual void execute(const std::string& src)       = 0;
    virtual std::string repl(const std::string& line)  = 0;
};
IScriptEngine& engine();   // installed by scripting module on bootstrap
}
