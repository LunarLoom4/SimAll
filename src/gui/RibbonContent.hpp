// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/RibbonContent.hpp
//
// Populates a `RibbonBar` with the 13-tab industrial layout from Section 3.8
// of the spec.  Returns a flat `RibbonActions` struct so MainWindow can wire
// signals into its slots without traversing the QToolBars.
//
// Listing each named action in one place lets us inspect / re-bind / hide
// individual buttons from automation tests and plugin code.
// =============================================================================
#pragma once

#include <QAction>
#include <QString>
#include <unordered_map>

namespace simall::gui
{

class RibbonBar;

struct RibbonActions
{
    std::unordered_map<QString, QAction*> byName;
    QAction* operator[](const QString& k) const
    {
        auto it = byName.find(k);
        return it == byName.end() ? nullptr : it->second;
    }
};

// Populate every tab with its default action set; returns the lookup map.
RibbonActions populate_default_ribbon(RibbonBar* ribbon);

} // namespace simall::gui
