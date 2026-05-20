// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/PluginRegistry.cpp
// Phase  : 17
// =============================================================================
#include "core/PluginRegistry.hpp"

#include <algorithm>

namespace simall::core
{

void PluginRegistry::register_plugin(plugins::IPlugin* p, Category cat)
{
    if (!p)
        return;
    std::lock_guard<std::mutex> g(mtx_);
    auto& bucket = buckets_[cat];
    if (std::find(bucket.begin(), bucket.end(), p) == bucket.end())
        bucket.push_back(p);
}

void PluginRegistry::register_auto(plugins::IPlugin* p)
{
    if (!p)
        return;
    if (auto* cp = dynamic_cast<plugins::CategorisedPlugin*>(p)) {
        register_plugin(p, cp->category());
    } else {
        register_plugin(p, Category::Generic);
    }
}

void PluginRegistry::unregister(plugins::IPlugin* p)
{
    if (!p)
        return;
    std::lock_guard<std::mutex> g(mtx_);
    for (auto& kv : buckets_) {
        auto& v = kv.second;
        v.erase(std::remove(v.begin(), v.end(), p), v.end());
    }
}

std::vector<plugins::IPlugin*> PluginRegistry::by_category(Category cat) const
{
    std::lock_guard<std::mutex> g(mtx_);
    auto it = buckets_.find(cat);
    return it == buckets_.end() ? std::vector<plugins::IPlugin*>{} : it->second;
}

std::vector<plugins::IPlugin*> PluginRegistry::all() const
{
    std::lock_guard<std::mutex> g(mtx_);
    std::vector<plugins::IPlugin*> out;
    for (const auto& kv : buckets_)
        out.insert(out.end(), kv.second.begin(), kv.second.end());
    return out;
}

std::size_t PluginRegistry::size() const
{
    std::lock_guard<std::mutex> g(mtx_);
    std::size_t n = 0;
    for (const auto& kv : buckets_)
        n += kv.second.size();
    return n;
}

void PluginRegistry::clear()
{
    std::lock_guard<std::mutex> g(mtx_);
    buckets_.clear();
}

PluginRegistry& PluginRegistry::instance()
{
    static PluginRegistry s;
    return s;
}

} // namespace simall::core
