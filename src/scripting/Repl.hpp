// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/Repl.hpp
//
// Python-aware line buffer / continuation detector.  The Qt
// `PythonConsolePanel` and the headless `simall` CLI both feed lines into a
// `Repl`; the REPL replies with one of three states for every line:
//
//   * Ready       — buffer holds a complete statement; flush and execute.
//   * NeedsMore   — open bracket / triple-quoted string / trailing `:`
//                   with no dedented follow-up line — emit a `... ` prompt.
//   * SyntaxError — unbalanced brackets that closed too many times.
//
// The detector is pure C++ — no Python interpreter required.  The Python
// engine is what eventually evaluates the accumulated buffer.
// =============================================================================
#pragma once

#include <deque>
#include <string>
#include <vector>

namespace simall::scripting {

enum class ReplState { Ready, NeedsMore, SyntaxError };

struct ReplFeedResult {
    ReplState    state;
    std::string  message;       // empty on Ready / NeedsMore
    std::string  promptHint;    // ">>> " or "... "
};

class Repl {
public:
    Repl() = default;

    // Feed a single physical line (newline terminator NOT included).
    ReplFeedResult feed(std::string line);

    // Return and clear the accumulated buffer (call when state==Ready).
    [[nodiscard]] std::string take_buffer();

    // Reset everything (Ctrl-C).
    void clear();

    // -- history --------------------------------------------------------------
    void push_history(std::string entry);
    [[nodiscard]] const std::vector<std::string>& history() const noexcept { return history_; }
    void clear_history() { history_.clear(); }

    // Up/Down navigation.  Returns the entry at the new cursor or empty
    // when the cursor reaches the bottom edge.
    [[nodiscard]] std::string history_prev();
    [[nodiscard]] std::string history_next();

    [[nodiscard]] bool empty() const noexcept { return buffer_.empty(); }

private:
    [[nodiscard]] ReplState classify() const;

    std::string              buffer_;        // accumulated lines, '\n' separated
    std::vector<std::string> history_;
    int                      historyCursor_ = -1;     // -1 == bottom
};

// -- exposed for unit-testing the analyser ----------------------------------
struct BracketState {
    int parens     = 0;   // ()
    int brackets   = 0;   // []
    int braces     = 0;   // {}
    bool inSingle  = false;
    bool inDouble  = false;
    bool inTriSingle = false;
    bool inTriDouble = false;
    bool trailingColon = false;
    bool trailingBackslash = false;
};
[[nodiscard]] BracketState analyse_python(std::string_view src);

}  // namespace simall::scripting
