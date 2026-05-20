// =============================================================================
// SimAll Beta — GUI Core (Qt-free)
// File   : src/gui_core/PanelRegistry.cpp
// =============================================================================
#include "gui_core/PanelRegistry.hpp"

namespace simall::gui_core
{

PanelRegistry& PanelRegistry::instance()
{
    static PanelRegistry s;
    return s;
}

void PanelRegistry::register_panel(PanelInfo info, PanelFactory factory)
{
    const std::string id = info.id;
    entries_[id] = Entry{std::move(info), std::move(factory)};
}

bool PanelRegistry::has(std::string_view id) const
{
    return entries_.find(std::string(id)) != entries_.end();
}

const PanelInfo* PanelRegistry::info(std::string_view id) const
{
    auto it = entries_.find(std::string(id));
    return it == entries_.end() ? nullptr : &it->second.info;
}

void* PanelRegistry::create(std::string_view id) const
{
    auto it = entries_.find(std::string(id));
    if (it == entries_.end() || !it->second.factory)
        return nullptr;
    return it->second.factory();
}

std::vector<PanelInfo> PanelRegistry::list() const
{
    std::vector<PanelInfo> out;
    out.reserve(entries_.size());
    for (auto& [_, e] : entries_)
        out.push_back(e.info);
    return out;
}

void PanelRegistry::clear()
{
    entries_.clear();
}

} // namespace simall::gui_core
