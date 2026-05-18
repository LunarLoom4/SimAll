// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/MacroRecorder.hpp
//
// Records `ScriptableCommand`s flowing through `CommandBus` and emits a
// Python script that, when replayed, reproduces the original sequence of
// GUI actions.  Provides start / stop / save / clear; the recorded log is
// kept in memory (a `std::vector<ScriptableCommand>`) so callers can patch
// it before serialisation.
// =============================================================================
#pragma once

#include "scripting/CommandBus.hpp"

#include <string>
#include <vector>

namespace simall::scripting {

class MacroRecorder {
public:
    MacroRecorder();
    ~MacroRecorder();

    void start();
    void stop();
    [[nodiscard]] bool recording() const noexcept { return recording_; }

    void clear();
    [[nodiscard]] size_t size() const noexcept { return commands_.size(); }
    [[nodiscard]] const std::vector<ScriptableCommand>& commands() const noexcept { return commands_; }

    // Manually append a command (useful for tests / synthetic playback).
    void append(ScriptableCommand cmd);

    // -- serialisation --------------------------------------------------------
    [[nodiscard]] std::string to_python_script(std::string_view header = {}) const;
    bool save_to_file(const std::string& path,
                      std::string_view header = {}) const;

private:
    int                                 listenerId_ = 0;
    bool                                recording_  = false;
    std::vector<ScriptableCommand>      commands_;
};

}  // namespace simall::scripting
