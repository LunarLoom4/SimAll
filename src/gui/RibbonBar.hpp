// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/RibbonBar.hpp
// Phase  : 2.1 / Section 3.8 of ultra-detailed spec.
//
// Ribbon implemented with QTabWidget + per-tab QToolBar rows. Heights and
// button dimensions follow Section 3.8 exactly (ribbon=110 px, large=72×72).
// =============================================================================
#pragma once

#include <QTabWidget>
#include <QToolButton>
#include <QAction>

namespace simall::gui {

class RibbonBar : public QTabWidget {
    Q_OBJECT
public:
    explicit RibbonBar(QWidget* parent = nullptr);

    QAction* add_button(const QString& tab, const QString& label, const QString& iconPath = {});

private:
    void add_tab(const QString& name);
};

}  // namespace simall::gui
