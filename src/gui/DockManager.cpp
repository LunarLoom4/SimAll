// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/DockManager.cpp
// =============================================================================
#include "gui/DockManager.hpp"

#include "core/Logger.hpp"

#include <QByteArray>
#include <QDockWidget>
#include <QMainWindow>
#include <QSettings>
#include <QWidget>

namespace simall::gui
{

namespace
{

Qt::DockWidgetArea to_qt(gui_core::DockArea a)
{
    switch (a) {
    case gui_core::DockArea::Left:
        return Qt::LeftDockWidgetArea;
    case gui_core::DockArea::Right:
        return Qt::RightDockWidgetArea;
    case gui_core::DockArea::Top:
        return Qt::TopDockWidgetArea;
    case gui_core::DockArea::Bottom:
        return Qt::BottomDockWidgetArea;
    default:
        return Qt::RightDockWidgetArea;
    }
}

gui_core::DockArea from_qt(Qt::DockWidgetArea a)
{
    switch (a) {
    case Qt::LeftDockWidgetArea:
        return gui_core::DockArea::Left;
    case Qt::RightDockWidgetArea:
        return gui_core::DockArea::Right;
    case Qt::TopDockWidgetArea:
        return gui_core::DockArea::Top;
    case Qt::BottomDockWidgetArea:
        return gui_core::DockArea::Bottom;
    default:
        return gui_core::DockArea::Hidden;
    }
}

} // namespace

DockManager::DockManager(QMainWindow* main, QObject* parent) : QObject(parent), main_(main) {}

QDockWidget* DockManager::add_panel(const QString& panelId,
                                    const QString& title,
                                    QWidget* panel,
                                    gui_core::DockArea area)
{
    auto* dock = new QDockWidget(title, main_);
    dock->setObjectName(panelId); // crucial: QMainWindow::saveState keys on this
    dock->setWidget(panel);
    main_->addDockWidget(to_qt(area), dock);
    docks_[panelId] = dock;
    return dock;
}

void DockManager::show_panel(const QString& panelId, bool visible)
{
    if (auto it = docks_.find(panelId); it != docks_.end())
        it->second->setVisible(visible);
}

bool DockManager::is_visible(const QString& panelId) const
{
    auto it = docks_.find(panelId);
    return it != docks_.end() && it->second->isVisible();
}

gui_core::DockPerspective DockManager::capture(const std::string& name) const
{
    gui_core::DockPerspective p;
    p.name = name;
    p.qtBlob = QString::fromLatin1(main_->saveState().toBase64()).toStdString();
    for (const auto& [id, dock] : docks_) {
        gui_core::DockPlacement pl;
        pl.panelId = id.toStdString();
        pl.visible = dock->isVisible();
        pl.floating = dock->isFloating();
        pl.area = pl.floating ? gui_core::DockArea::Floating : from_qt(main_->dockWidgetArea(dock));
        const QRect g = dock->geometry();
        pl.x = g.x();
        pl.y = g.y();
        pl.w = g.width();
        pl.h = g.height();
        p.placements.push_back(std::move(pl));
    }
    return p;
}

void DockManager::apply_perspective(const gui_core::DockPerspective& p)
{
    // Restore opaque Qt state first (splitter ratios, tabify groups).
    if (!p.qtBlob.empty()) {
        QByteArray blob = QByteArray::fromBase64(QByteArray::fromStdString(p.qtBlob));
        main_->restoreState(blob);
    }
    // Apply explicit per-panel toggles afterwards in case the blob is stale.
    for (const auto& pl : p.placements) {
        auto it = docks_.find(QString::fromStdString(pl.panelId));
        if (it == docks_.end())
            continue;
        QDockWidget* d = it->second;
        d->setVisible(pl.visible);
        if (pl.floating) {
            d->setFloating(true);
            if (pl.w > 0 && pl.h > 0)
                d->setGeometry(pl.x, pl.y, pl.w, pl.h);
        } else if (pl.area != gui_core::DockArea::Hidden
                   && pl.area != gui_core::DockArea::Floating) {
            main_->addDockWidget(to_qt(pl.area), d);
            d->setFloating(false);
        }
    }
}

void DockManager::save_perspective(const std::string& name)
{
    perspectives_[name] = capture(name);
}

bool DockManager::load_perspective(const std::string& name)
{
    auto it = perspectives_.find(name);
    if (it == perspectives_.end())
        return false;
    apply_perspective(it->second);
    return true;
}

std::vector<std::string> DockManager::perspective_names() const
{
    std::vector<std::string> n;
    n.reserve(perspectives_.size());
    for (auto& [k, _] : perspectives_)
        n.push_back(k);
    return n;
}

void DockManager::persist_to_settings(const QString& key) const
{
    QSettings s;
    s.beginGroup(key);
    s.setValue("state", main_->saveState());
    s.setValue("geometry", main_->saveGeometry());
    QStringList names;
    for (auto& [k, p] : perspectives_) {
        names << QString::fromStdString(k);
        s.setValue("perspective/" + QString::fromStdString(k),
                   QString::fromStdString(gui_core::serialize(p)));
    }
    s.setValue("perspective_names", names);
    s.endGroup();
}

void DockManager::restore_from_settings(const QString& key)
{
    QSettings s;
    s.beginGroup(key);
    if (s.contains("geometry"))
        main_->restoreGeometry(s.value("geometry").toByteArray());
    if (s.contains("state"))
        main_->restoreState(s.value("state").toByteArray());
    perspectives_.clear();
    const QStringList names = s.value("perspective_names").toStringList();
    for (const QString& n : names) {
        const QString txt = s.value("perspective/" + n).toString();
        auto p = gui_core::deserialize(txt.toStdString());
        if (p)
            perspectives_[n.toStdString()] = std::move(*p);
    }
    s.endGroup();
}

} // namespace simall::gui
