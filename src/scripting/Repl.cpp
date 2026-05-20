// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/Repl.cpp
// =============================================================================
#include "scripting/Repl.hpp"

#include <cctype>

namespace simall::scripting
{

BracketState analyse_python(std::string_view src)
{
    BracketState s;
    bool lastWasBackslash = false;
    size_t lastNonWs = std::string::npos;
    for (size_t i = 0; i < src.size(); ++i) {
        const char c = src[i];
        const bool inTri = s.inTriSingle || s.inTriDouble;
        if (!s.inSingle && !s.inDouble && !inTri && c == '#') {
            // comment to end of line
            while (i < src.size() && src[i] != '\n')
                ++i;
            continue;
        }
        // Triple quote detection.
        if (!s.inSingle && !s.inDouble) {
            if (i + 2 < src.size() && src[i] == '"' && src[i + 1] == '"' && src[i + 2] == '"') {
                s.inTriDouble = !s.inTriDouble;
                i += 2;
                continue;
            }
            if (i + 2 < src.size() && src[i] == '\'' && src[i + 1] == '\'' && src[i + 2] == '\'') {
                s.inTriSingle = !s.inTriSingle;
                i += 2;
                continue;
            }
        }
        if (!s.inTriSingle && !s.inTriDouble) {
            if (!s.inDouble && c == '\'' && !lastWasBackslash) {
                s.inSingle = !s.inSingle;
            } else if (!s.inSingle && c == '"' && !lastWasBackslash) {
                s.inDouble = !s.inDouble;
            }
        }
        lastWasBackslash = (c == '\\' && !lastWasBackslash);

        if (!s.inSingle && !s.inDouble && !s.inTriSingle && !s.inTriDouble) {
            switch (c) {
            case '(':
                ++s.parens;
                break;
            case ')':
                --s.parens;
                break;
            case '[':
                ++s.brackets;
                break;
            case ']':
                --s.brackets;
                break;
            case '{':
                ++s.braces;
                break;
            case '}':
                --s.braces;
                break;
            default:
                break;
            }
        }
        if (!std::isspace(static_cast<unsigned char>(c)))
            lastNonWs = i;
    }
    if (lastNonWs != std::string::npos) {
        const char ch = src[lastNonWs];
        s.trailingColon = (ch == ':');
        s.trailingBackslash = (ch == '\\');
    }
    return s;
}

ReplState Repl::classify() const
{
    auto s = analyse_python(buffer_);
    if (s.parens < 0 || s.brackets < 0 || s.braces < 0)
        return ReplState::SyntaxError;
    if (s.parens > 0 || s.brackets > 0 || s.braces > 0)
        return ReplState::NeedsMore;
    if (s.inTriSingle || s.inTriDouble)
        return ReplState::NeedsMore;
    if (s.trailingBackslash)
        return ReplState::NeedsMore;
    if (s.trailingColon) {
        // Block-statement: need at least one more line.  We treat an empty
        // trailing line as the terminator (matches CPython REPL behaviour).
        // If the buffer's last *physical* line is blank, the block ends.
        size_t lastNl = buffer_.find_last_of('\n');
        if (lastNl == std::string::npos)
            return ReplState::NeedsMore;
        std::string_view tail(buffer_.data() + lastNl + 1, buffer_.size() - lastNl - 1);
        for (char c : tail)
            if (!std::isspace(static_cast<unsigned char>(c)))
                return ReplState::NeedsMore;
        return ReplState::Ready;
    }
    return ReplState::Ready;
}

ReplFeedResult Repl::feed(std::string line)
{
    if (!buffer_.empty())
        buffer_.push_back('\n');
    buffer_ += line;
    const ReplState st = classify();
    ReplFeedResult r;
    r.state = st;
    switch (st) {
    case ReplState::Ready:
        r.promptHint = ">>> ";
        break;
    case ReplState::NeedsMore:
        r.promptHint = "... ";
        break;
    case ReplState::SyntaxError:
        r.message = "Unbalanced closing bracket";
        r.promptHint = ">>> ";
        buffer_.clear();
        break;
    }
    return r;
}

std::string Repl::take_buffer()
{
    std::string out;
    out.swap(buffer_);
    return out;
}

void Repl::clear()
{
    buffer_.clear();
    historyCursor_ = -1;
}

void Repl::push_history(std::string entry)
{
    if (entry.empty())
        return;
    if (!history_.empty() && history_.back() == entry) {
        historyCursor_ = -1;
        return;
    }
    history_.push_back(std::move(entry));
    historyCursor_ = -1;
}

std::string Repl::history_prev()
{
    if (history_.empty())
        return {};
    if (historyCursor_ == -1)
        historyCursor_ = int(history_.size()) - 1;
    else if (historyCursor_ > 0)
        --historyCursor_;
    return history_[historyCursor_];
}

std::string Repl::history_next()
{
    if (history_.empty() || historyCursor_ == -1)
        return {};
    if (historyCursor_ + 1 < int(history_.size())) {
        ++historyCursor_;
        return history_[historyCursor_];
    }
    historyCursor_ = -1;
    return {};
}

} // namespace simall::scripting
