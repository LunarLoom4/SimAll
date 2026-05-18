#include "scripting/ScriptEngine.hpp"
namespace simall::scripting {
namespace {
class NullEngine : public IScriptEngine {
public:
    void execute(const std::string&) override {}
    std::string repl(const std::string&) override { return "scripting backend not yet bound"; }
};
}
IScriptEngine& engine() { static NullEngine e; return e; }
}
