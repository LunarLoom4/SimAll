// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/CellAdapterRegistry.cpp
// Phase  : 22 Pass 22.5
// =============================================================================
#include "workbench/CellAdapterRegistry.hpp"

#include <algorithm>

namespace simall::workbench
{

bool CellAdapterRegistry::register_factory(std::string id, Factory f)
{
    std::lock_guard<std::mutex> lk(mu_);
    const bool fresh = factories_.find(id) == factories_.end();
    factories_[std::move(id)] = std::move(f);
    return fresh;
}

bool CellAdapterRegistry::unregister(const std::string& id)
{
    std::lock_guard<std::mutex> lk(mu_);
    return factories_.erase(id) != 0;
}

bool CellAdapterRegistry::has(const std::string& id) const
{
    std::lock_guard<std::mutex> lk(mu_);
    return factories_.find(id) != factories_.end();
}

std::unique_ptr<ICellAdapter> CellAdapterRegistry::make(const std::string& id) const
{
    Factory copy;
    {
        std::lock_guard<std::mutex> lk(mu_);
        auto it = factories_.find(id);
        if (it == factories_.end())
            return nullptr;
        copy = it->second; // copy under lock; invoke outside
    }
    if (!copy)
        return nullptr;
    return copy();
}

std::vector<std::string> CellAdapterRegistry::keys() const
{
    std::vector<std::string> out;
    {
        std::lock_guard<std::mutex> lk(mu_);
        out.reserve(factories_.size());
        for (const auto& kv : factories_)
            out.push_back(kv.first);
    }
    std::sort(out.begin(), out.end());
    return out;
}

std::size_t CellAdapterRegistry::size() const
{
    std::lock_guard<std::mutex> lk(mu_);
    return factories_.size();
}

bool CellAdapterRegistry::empty() const
{
    std::lock_guard<std::mutex> lk(mu_);
    return factories_.empty();
}

void CellAdapterRegistry::clear()
{
    std::lock_guard<std::mutex> lk(mu_);
    factories_.clear();
}

} // namespace simall::workbench
