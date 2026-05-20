// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/RefreshPolicy.hpp
// Phase  : 22 Pass 22.2 (state-propagation engine)
//
// RefreshPolicy -- the rule the WorkflowEngine uses when asked to compute
// a "what should I refresh next?" plan.  Pure enum + helpers; no state.
// =============================================================================
#pragma once

#include <cstdint>
#include <string_view>

namespace simall::workbench
{

// ---------------------------------------------------------------------------
// What does the user mean when they click "Refresh"?
//
//   ReadyOnly       : only cells whose required inputs are wired and all
//                     upstream cells are already UpToDate.  This is the
//                     default; matches Workbench's "Refresh" button.
//   IncludeStale    : also include RefreshRequired cells whose upstream is
//                     itself RefreshRequired -- caller must run parents
//                     first.  Useful for "Update Project" semantics where
//                     the engine recurses upstream automatically.
//   StopOnFailed    : never include cells that have any transitive Failed
//                     ancestor.  Combine with either of the above.
// ---------------------------------------------------------------------------
enum class RefreshPolicy : std::uint8_t
{
    ReadyOnly = 0,
    IncludeStale = 1 << 0,
    StopOnFailed = 1 << 1,
    UpdateProject = IncludeStale | StopOnFailed, // convenient combo
};

[[nodiscard]] constexpr bool has_flag(RefreshPolicy p, RefreshPolicy f) noexcept
{
    using U = std::uint8_t;
    return (static_cast<U>(p) & static_cast<U>(f)) == static_cast<U>(f);
}

[[nodiscard]] std::string_view to_string(RefreshPolicy) noexcept;

} // namespace simall::workbench
