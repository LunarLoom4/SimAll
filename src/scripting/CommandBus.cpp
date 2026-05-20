// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/CommandBus.cpp
// =============================================================================
#include "scripting/CommandBus.hpp"

#include <algorithm>
#include <sstream>

namespace simall::scripting
{

namespace
{
bool needs_quoting(std::string_view v)
{
    if (v.empty())
        return true;
    if (v == "true" || v == "false" || v == "None")
        return false;
    bool numeric = true;
    int dotCount = 0;
    for (size_t i = 0; i < v.size(); ++i) {
        const char c = v[i];
        if (c == '-' && i == 0)
            continue;
        if (c == '.') {
            if (++dotCount > 1) {
                numeric = false;
                break;
            }
            continue;
        }
        if (c == 'e' || c == 'E')
            continue;
        if (c == '+')
            continue;
        if (!std::isdigit(static_cast<unsigned char>(c))) {
            numeric = false;
            break;
        }
    }
    return !numeric;
}

std::string python_repr(std::string_view v)
{
    if (!needs_quoting(v))
        return std::string(v);
    std::string out;
    out.reserve(v.size() + 2);
    out.push_back('"');
    for (char c : v) {
        if (c == '\\' || c == '"')
            out.push_back('\\');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}
} // namespace

std::string ScriptableCommand::to_python_call(std::string_view module) const
{
    std::ostringstream os;
    if (!module.empty())
        os << module << '.';
    os << name << '(';
    bool first = true;
    for (auto& [k, v] : args) {
        if (!first)
            os << ", ";
        first = false;
        os << k << '=' << python_repr(v);
    }
    os << ')';
    return os.str();
}

CommandBus& CommandBus::instance()
{
    static CommandBus s;
    return s;
}

void CommandBus::register_handler(std::string name, CommandHandler h)
{
    std::lock_guard lk(mu_);
    handlers_[std::move(name)] = std::move(h);
}

bool CommandBus::has_handler(std::string_view name) const
{
    std::lock_guard lk(mu_);
    return handlers_.find(std::string(name)) != handlers_.end();
}

void CommandBus::unregister_handler(std::string_view name)
{
    std::lock_guard lk(mu_);
    handlers_.erase(std::string(name));
}

int CommandBus::add_listener(CommandListener l)
{
    std::lock_guard lk(mu_);
    const int id = nextListenerId_++;
    listeners_[id] = std::move(l);
    return id;
}

void CommandBus::remove_listener(int id)
{
    std::lock_guard lk(mu_);
    listeners_.erase(id);
}

bool CommandBus::dispatch(const ScriptableCommand& cmd)
{
    CommandHandler handler;
    std::vector<CommandListener> listenerCopy;
    {
        std::lock_guard lk(mu_);
        auto it = handlers_.find(cmd.name);
        if (it != handlers_.end())
            handler = it->second;
        listenerCopy.reserve(listeners_.size());
        for (auto& [_, l] : listeners_)
            listenerCopy.push_back(l);
    }
    for (auto& l : listenerCopy) {
        try {
            l(cmd);
        } catch (...) {
        }
    }
    if (!handler)
        return false;
    try {
        handler(cmd);
    } catch (...) {
    }
    return true;
}

void CommandBus::clear()
{
    std::lock_guard lk(mu_);
    handlers_.clear();
    listeners_.clear();
}

std::vector<std::string> CommandBus::registered_commands() const
{
    std::lock_guard lk(mu_);
    std::vector<std::string> n;
    n.reserve(handlers_.size());
    for (auto& [k, _] : handlers_)
        n.push_back(k);
    std::sort(n.begin(), n.end());
    return n;
}

} // namespace simall::scripting
