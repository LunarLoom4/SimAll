// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/ConsolePanel.cpp
// =============================================================================
#include "gui/ConsolePanel.hpp"

#include <QComboBox>
#include <QDateTime>
#include <QHBoxLayout>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QToolButton>
#include <QVBoxLayout>
#include <QTextCursor>
#include <QScrollBar>

namespace simall::gui {

ConsolePanel::ConsolePanel(QWidget* parent) : QWidget(parent) {
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(4, 4, 4, 4);
    root->setSpacing(4);

    auto* toolbar = new QHBoxLayout();
    toolbar->setSpacing(4);

    levelCombo_ = new QComboBox(this);
    levelCombo_->addItems({"Trace", "Debug", "Info", "Warning", "Error", "Fatal"});
    levelCombo_->setCurrentIndex(2);
    toolbar->addWidget(levelCombo_);

    search_ = new QLineEdit(this);
    search_->setPlaceholderText("Filter (substring)…");
    toolbar->addWidget(search_, /*stretch*/1);

    autoscrollBtn_ = new QToolButton(this);
    autoscrollBtn_->setText("Auto-scroll");
    autoscrollBtn_->setCheckable(true);
    autoscrollBtn_->setChecked(true);
    toolbar->addWidget(autoscrollBtn_);

    clearBtn_ = new QToolButton(this);
    clearBtn_->setText("Clear");
    toolbar->addWidget(clearBtn_);

    root->addLayout(toolbar);

    view_ = new QPlainTextEdit(this);
    view_->setReadOnly(true);
    view_->setMaximumBlockCount(kMaxEntries);
    view_->setStyleSheet(
        "QPlainTextEdit{background:#1B1B1B;color:#D4D4D4;"
        "font-family:Consolas,Menlo,'Courier New',monospace;font-size:10pt;}");
    root->addWidget(view_, /*stretch*/1);

    connect(this, &ConsolePanel::newEntry, this, &ConsolePanel::rebuild_view,
            Qt::QueuedConnection);
    connect(levelCombo_, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, [this](int i){ levelThreshold_ = i; rebuild_view(); });
    connect(search_, &QLineEdit::textChanged, this, [this]{ rebuild_view(); });
    connect(autoscrollBtn_, &QToolButton::toggled, this,
            [this](bool on){ autoscroll_ = on; });
    connect(clearBtn_, &QToolButton::clicked, this, &ConsolePanel::clear_history);
}

void ConsolePanel::post(int level, const QString& category, const QString& msg) {
    {
        std::lock_guard lk(mu_);
        entries_.push_back({level, category, msg,
                            QDateTime::currentDateTime().toString("hh:mm:ss.zzz")});
        while (entries_.size() > kMaxEntries) entries_.pop_front();
    }
    emit newEntry();
}

void ConsolePanel::set_level_threshold(int level) {
    levelThreshold_ = level;
    if (levelCombo_) levelCombo_->setCurrentIndex(level);
    rebuild_view();
}

void ConsolePanel::clear_history() {
    {
        std::lock_guard lk(mu_);
        entries_.clear();
    }
    view_->clear();
}

bool ConsolePanel::passes_filter(const Entry& e) const {
    if (e.level < levelThreshold_) return false;
    const QString needle = search_ ? search_->text() : QString();
    if (needle.isEmpty()) return true;
    return e.text.contains(needle, Qt::CaseInsensitive) ||
           e.category.contains(needle, Qt::CaseInsensitive);
}

QString ConsolePanel::colour_for(int level) const {
    switch (level) {
        case 0: return "#888888";  // Trace
        case 1: return "#6E8FE8";  // Debug
        case 2: return "#D4D4D4";  // Info
        case 3: return "#E0A600";  // Warning
        case 4: return "#E52B50";  // Error
        case 5: return "#FF55FF";  // Fatal
        default: return "#D4D4D4";
    }
}

QString ConsolePanel::format_line(const Entry& e) const {
    return QString("<span style='color:%1'>[%2] [%3] %4</span>")
        .arg(colour_for(e.level), e.timestamp, e.category, e.text.toHtmlEscaped());
}

void ConsolePanel::append_entry_to_view(const Entry& e) {
    view_->appendHtml(format_line(e));
}

void ConsolePanel::rebuild_view() {
    std::deque<Entry> snapshot;
    {
        std::lock_guard lk(mu_);
        snapshot = entries_;
    }
    view_->clear();
    for (const auto& e : snapshot)
        if (passes_filter(e)) append_entry_to_view(e);
    if (autoscroll_) {
        auto* sb = view_->verticalScrollBar();
        sb->setValue(sb->maximum());
    }
}

}  // namespace simall::gui
