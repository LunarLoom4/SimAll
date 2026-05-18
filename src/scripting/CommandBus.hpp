// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/CommandBus.hpp
//
// Typed command object + thread-safe dispatcher.  Every GUI action (ribbon
// button, dialog OK, dock toggle) emits a `ScriptableCommand`, which the
// dispatcher routes to registered handlers AND to listeners (used by
// `MacroRecorder`).  Scripting code can synthesize commands programmatically
// to drive the application from Python.
// =============================================================================
#pragma once

#include <functional>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::scripting {

struct ScriptableCommand {
    std::string                                       name;          // "Mesh.Surface"
    std::vector<std::pair<std::string, std::string>>  args;          // ordered
    // Stable text form for macro replay: "Mesh.Surface(size=0.01, growth=1.2)"
    [[nodiscard]] std::string to_python_call(std::string_view module = "simall") const;
};

using CommandHandler  = std::function<void(const ScriptableCommand&)>;
using CommandListener = std::function<void(const ScriptableCommand&)>;

class CommandBus {
public:
    static CommandBus& instance();

    // Handlers receive ONE command (first registered handler wins).
    void  register_handler(std::string name, CommandHandler h);
    bool  has_handler(std::string_view name) const;
    void  unregister_handler(std::string_view name);

    // Listeners receive EVERY command (used by macro recorder).
    int   add_listener(CommandListener l);
    void  remove_listener(int id);

    // Dispatch.  Returns true iff a handler ran (listeners still fire).
    bool  dispatch(const ScriptableCommand& cmd);

    void  clear();
    [[nodiscard]] std::vector<std::string> registered_commands() const;

private:
    CommandBus() = default;
    mutable std::mutex                            mu_;
    std::unordered_map<std::string, CommandHandler>  handlers_;
    std::unordered_map<int, CommandListener>         listeners_;
    int                                              nextListenerId_ = 1;
};

}  // namespace simall::scripting
