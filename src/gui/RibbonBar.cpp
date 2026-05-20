#include "gui/RibbonBar.hpp"

#include <QToolBar>

#include <QHash>
#include <QVBoxLayout>
#include <QWidget>

namespace simall::gui
{

namespace
{
constexpr int kRibbonHeight = 110; // Section 3.8
constexpr int kBigBtn = 72;
}

class RibbonTab : public QWidget
{
public:
    explicit RibbonTab(QWidget* parent = nullptr) : QWidget(parent)
    {
        auto* lay = new QHBoxLayout(this);
        lay->setContentsMargins(6, 4, 6, 4);
        lay->setSpacing(4);
        bar_ = new QToolBar(this);
        bar_->setIconSize({32, 32});
        bar_->setToolButtonStyle(Qt::ToolButtonTextUnderIcon);
        lay->addWidget(bar_);
        lay->addStretch();
    }
    QToolBar* bar() { return bar_; }

private:
    QToolBar* bar_;
};

RibbonBar::RibbonBar(QWidget* parent) : QTabWidget(parent)
{
    setFixedHeight(kRibbonHeight);
    setDocumentMode(true);
    add_tab("File");
    add_tab("Geometry");
    add_tab("Mesh");
    add_tab("Physics");
    add_tab("Materials");
    add_tab("Boundary Conditions");
    add_tab("Solver");
    add_tab("Initialization");
    add_tab("Run");
    add_tab("Results");
    add_tab("Automation");
    add_tab("HPC");
    add_tab("Plugins");
}

void RibbonBar::add_tab(const QString& name)
{
    addTab(new RibbonTab(this), name);
}

QAction* RibbonBar::add_button(const QString& tab, const QString& label, const QString& iconPath)
{
    for (int i = 0; i < count(); ++i) {
        if (tabText(i) == tab) {
            auto* t = qobject_cast<RibbonTab*>(widget(i));
            if (!t)
                return nullptr;
            QAction* a = t->bar()->addAction(label);
            if (!iconPath.isEmpty())
                a->setIcon(QIcon(iconPath));
            return a;
        }
    }
    return nullptr;
}

} // namespace simall::gui
