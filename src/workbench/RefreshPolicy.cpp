// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/RefreshPolicy.cpp
// Phase  : 22 Pass 22.2
// =============================================================================
#include "workbench/RefreshPolicy.hpp"

namespace simall::workbench
{

std::string_view to_string(RefreshPolicy p) noexcept
{
    switch (p) {
    case RefreshPolicy::ReadyOnly:
        return "ReadyOnly";
    case RefreshPolicy::IncludeStale:
        return "IncludeStale";
    case RefreshPolicy::StopOnFailed:
        return "StopOnFailed";
    case RefreshPolicy::UpdateProject:
        return "UpdateProject";
    }
    return "<custom RefreshPolicy>";
}

} // namespace simall::workbench
