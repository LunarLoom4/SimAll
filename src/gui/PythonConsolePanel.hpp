// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/PythonConsolePanel.hpp
//
// Read-Eval-Print interface for the embedded scripting layer (W16 hooks).
// In W15 the panel is a fully functional shell: history, multi-line input,
// banner, and a pluggable `Executor` callback that returns the printed
// output for a typed command.  The pybind11-backed executor is registered
// at startup in W16 — until then we ship a `NullExecutor` that returns
// "Python not initialised".
// =============================================================================
#pragma once

#include <functional>
#include <QWidget>
#include <vector>

class QPlainTextEdit;
class QLineEdit;
class QPushButton;

namespace simall::gui
{

using PythonExecutor = std::function<QString(const QString&)>;

class PythonConsolePanel : public QWidget
{
    Q_OBJECT
public:
    explicit PythonConsolePanel(QWidget* parent = nullptr);

    void set_executor(PythonExecutor exec);
    void execute(const QString& cmd);
    void clear_console();
    // Append text from an out-of-band source (the embedded interpreter's
    // captured stdout / stderr).  Safe to invoke via QueuedConnection.
    void post_output(const QString& text);
    [[nodiscard]] QStringList history() const { return history_; }

private slots:
    void on_return_pressed();
    void on_key_up();
    void on_key_down();

private:
    void print_banner();
    void append(const QString& text, const QString& colour = QString());

    QPlainTextEdit* output_ = nullptr;
    QLineEdit* prompt_ = nullptr;
    QPushButton* runBtn_ = nullptr;
    QPushButton* clearBtn_ = nullptr;

    PythonExecutor executor_;
    QStringList history_;
    int historyIdx_ = -1;
};

} // namespace simall::gui
