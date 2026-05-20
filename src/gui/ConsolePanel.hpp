// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/ConsolePanel.hpp
//
// Log console with level filter, substring search, auto-scroll toggle, and
// a clear button.  Receives entries from `core::Logger` via a `LogSink` that
// marshals into the GUI thread.  Holds an in-memory ring buffer so filter
// toggles re-render without losing history.
// =============================================================================
#pragma once

#include <deque>
#include <mutex>
#include <QWidget>

class QPlainTextEdit;
class QLineEdit;
class QComboBox;
class QToolButton;

namespace simall::gui
{

class ConsolePanel : public QWidget
{
    Q_OBJECT
public:
    explicit ConsolePanel(QWidget* parent = nullptr);

    // Thread-safe entry point — safe to call from worker threads.
    void post(int level, const QString& category, const QString& msg);

    void set_level_threshold(int level);
    void clear_history();

signals:
    void newEntry(); // queued connection trigger so render runs on UI thread

private:
    struct Entry
    {
        int level;
        QString category;
        QString text;
        QString timestamp;
    };
    void rebuild_view();
    void append_entry_to_view(const Entry& e);
    [[nodiscard]] bool passes_filter(const Entry& e) const;
    [[nodiscard]] QString format_line(const Entry& e) const;
    [[nodiscard]] QString colour_for(int level) const;

    QPlainTextEdit* view_ = nullptr;
    QLineEdit* search_ = nullptr;
    QComboBox* levelCombo_ = nullptr;
    QToolButton* autoscrollBtn_ = nullptr;
    QToolButton* clearBtn_ = nullptr;

    mutable std::mutex mu_;
    std::deque<Entry> entries_; // bounded
    int levelThreshold_ = 0;
    bool autoscroll_ = true;
    static constexpr int kMaxEntries = 4096;
};

} // namespace simall::gui
