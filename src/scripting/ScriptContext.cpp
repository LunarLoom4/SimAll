// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/ScriptContext.cpp
// =============================================================================
#include "scripting/ScriptContext.hpp"

#include <algorithm>
#include <sstream>

namespace simall::scripting {

ScriptContext& ScriptContext::instance() {
    static ScriptContext s;
    return s;
}

void ScriptContext::set(std::string key, std::string value) {
    std::vector<std::pair<int, Subscriber>> notifyList;
    std::string oldValue;
    {
        std::lock_guard lk(mu_);
        auto it = bag_.find(key);
        if (it != bag_.end()) oldValue = it->second;
        bag_[key] = value;
        notifyList.assign(subs_.begin(), subs_.end());
    }
    for (auto& [_, s] : notifyList) {
        try { s(key, oldValue, value); } catch (...) {}
    }
}

std::optional<std::string> ScriptContext::get(std::string_view key) const {
    std::lock_guard lk(mu_);
    auto it = bag_.find(std::string(key));
    if (it == bag_.end()) return std::nullopt;
    return it->second;
}

std::string ScriptContext::get_or(std::string_view key,
                                   std::string_view fallback) const {
    if (auto v = get(key)) return *v;
    return std::string(fallback);
}

bool ScriptContext::erase(std::string_view key) {
    std::lock_guard lk(mu_);
    return bag_.erase(std::string(key)) > 0;
}

std::vector<std::string> ScriptContext::keys() const {
    std::lock_guard lk(mu_);
    std::vector<std::string> out;
    out.reserve(bag_.size());
    for (auto& [k, _] : bag_) out.push_back(k);
    std::sort(out.begin(), out.end());
    return out;
}

void ScriptContext::clear() {
    std::lock_guard lk(mu_);
    bag_.clear();
}

void ScriptContext::set_working_directory(std::string path) {
    std::lock_guard lk(mu_);
    cwd_ = std::move(path);
}

std::string ScriptContext::working_directory() const {
    std::lock_guard lk(mu_);
    return cwd_;
}

int ScriptContext::subscribe(Subscriber s) {
    std::lock_guard lk(mu_);
    const int id = nextSubId_++;
    subs_[id] = std::move(s);
    return id;
}

void ScriptContext::unsubscribe(int id) {
    std::lock_guard lk(mu_);
    subs_.erase(id);
}

std::string ScriptContext::deterministic_snapshot() const {
    std::lock_guard lk(mu_);
    std::vector<std::pair<std::string,std::string>> items(bag_.begin(), bag_.end());
    std::sort(items.begin(), items.end());
    std::ostringstream os;
    for (auto& [k, v] : items) os << k << '=' << v << '\n';
    return os.str();
}

}  // namespace simall::scripting
