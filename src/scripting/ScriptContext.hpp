// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/ScriptContext.hpp
//
// Pure-C++ environment shared by every scripting front-end (the GUI
// `PythonConsolePanel`, the UDF host, the macro recorder, and the embedded
// REPL).  The context keeps a string-keyed bag of stringly-typed values
// (e.g. "project.path", "solver.algorithm", "current.zone") so scripting
// callers and headless tests can inspect / mutate the same state without
// linking pybind11.
//
// Concrete pybind11 bindings (PythonEngine) attach to this object, expose
// it as `simall.context` in the `simall` Python module, and route Python
// dict reads/writes through `get`/`set`.
// =============================================================================
#pragma once

#include <chrono>
#include <functional>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::scripting
{

class ScriptContext
{
public:
    static ScriptContext& instance();

    // -- value bag ------------------------------------------------------------
    void set(std::string key, std::string value);
    [[nodiscard]] std::optional<std::string> get(std::string_view key) const;
    [[nodiscard]] std::string get_or(std::string_view key, std::string_view fallback) const;
    bool erase(std::string_view key);
    [[nodiscard]] std::vector<std::string> keys() const;
    void clear();

    // -- working directory ----------------------------------------------------
    void set_working_directory(std::string path);
    [[nodiscard]] std::string working_directory() const;

    // -- change subscription --------------------------------------------------
    // Subscribers receive (key, oldValue, newValue).  Used by the GUI to
    // refresh affected dialogs and by the macro recorder to log writes.
    using Subscriber =
        std::function<void(const std::string&, const std::string&, const std::string&)>;
    int subscribe(Subscriber s);
    void unsubscribe(int id);

    // -- snapshot (deterministic for hashing in tests) ------------------------
    [[nodiscard]] std::string deterministic_snapshot() const;

private:
    ScriptContext() = default;
    mutable std::mutex mu_;
    std::unordered_map<std::string, std::string> bag_;
    std::string cwd_;
    std::unordered_map<int, Subscriber> subs_;
    int nextSubId_ = 1;
};

} // namespace simall::scripting
