// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/PythonConsolePanel.cpp
// =============================================================================
#include "gui/PythonConsolePanel.hpp"

#include <QHBoxLayout>
#include <QKeyEvent>
#include <QLineEdit>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QScrollBar>
#include <QVBoxLayout>

namespace simall::gui
{

namespace
{
class HistoryLineEdit : public QLineEdit
{
public:
    using QLineEdit::QLineEdit;
    std::function<void()> upCb, downCb;

protected:
    void keyPressEvent(QKeyEvent* ev) override
    {
        if (ev->key() == Qt::Key_Up) {
            if (upCb)
                upCb();
            return;
        }
        if (ev->key() == Qt::Key_Down) {
            if (downCb)
                downCb();
            return;
        }
        QLineEdit::keyPressEvent(ev);
    }
};
} // namespace

PythonConsolePanel::PythonConsolePanel(QWidget* parent) : QWidget(parent)
{
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(4, 4, 4, 4);
    root->setSpacing(4);

    output_ = new QPlainTextEdit(this);
    output_->setReadOnly(true);
    output_->setStyleSheet("QPlainTextEdit{background:#101010;color:#D4D4D4;"
                           "font-family:Consolas,Menlo,'Courier New',monospace;font-size:10pt;}");
    root->addWidget(output_, 1);

    auto* row = new QHBoxLayout();
    auto* prompt = new HistoryLineEdit(this);
    prompt->setPlaceholderText(">>> Python command");
    prompt->upCb = [this] { on_key_up(); };
    prompt->downCb = [this] { on_key_down(); };
    prompt_ = prompt;
    runBtn_ = new QPushButton("Run", this);
    clearBtn_ = new QPushButton("Clear", this);
    row->addWidget(prompt_, 1);
    row->addWidget(runBtn_);
    row->addWidget(clearBtn_);
    root->addLayout(row);

    setMinimumHeight(180);
    print_banner();

    connect(prompt_, &QLineEdit::returnPressed, this, &PythonConsolePanel::on_return_pressed);
    connect(runBtn_, &QPushButton::clicked, this, &PythonConsolePanel::on_return_pressed);
    connect(clearBtn_, &QPushButton::clicked, this, &PythonConsolePanel::clear_console);

    executor_ = [](const QString&) {
        return QString("(Python not initialised — register an executor in W16)");
    };
}

void PythonConsolePanel::set_executor(PythonExecutor exec)
{
    executor_ = std::move(exec);
}

void PythonConsolePanel::print_banner()
{
    append("SimAll Beta — Python Console", "#9CDCFE");
    append("Type a command and press Enter.  ↑/↓ navigate history.", "#888888");
}

void PythonConsolePanel::execute(const QString& cmd)
{
    if (cmd.trimmed().isEmpty())
        return;
    history_ << cmd;
    historyIdx_ = history_.size();
    append(">>> " + cmd, "#9CDCFE");
    if (cmd.trimmed() == "clear") {
        clear_console();
        return;
    }
    QString out;
    try {
        out = executor_ ? executor_(cmd) : QString();
    } catch (const std::exception& ex) {
        out = QString("Exception: ") + ex.what();
    }
    if (!out.isEmpty())
        append(out, "#D4D4D4");
}

void PythonConsolePanel::on_return_pressed()
{
    const QString cmd = prompt_->text();
    prompt_->clear();
    execute(cmd);
}

void PythonConsolePanel::on_key_up()
{
    if (history_.isEmpty())
        return;
    if (historyIdx_ > 0)
        --historyIdx_;
    prompt_->setText(history_.value(historyIdx_));
}

void PythonConsolePanel::on_key_down()
{
    if (history_.isEmpty())
        return;
    if (historyIdx_ < history_.size() - 1) {
        ++historyIdx_;
        prompt_->setText(history_.value(historyIdx_));
    } else {
        historyIdx_ = history_.size();
        prompt_->clear();
    }
}

void PythonConsolePanel::clear_console()
{
    output_->clear();
    print_banner();
}

void PythonConsolePanel::post_output(const QString& text)
{
    if (text.isEmpty())
        return;
    append(text, "#B5CEA8");
}

void PythonConsolePanel::append(const QString& text, const QString& colour)
{
    QString html;
    if (colour.isEmpty())
        html = text.toHtmlEscaped();
    else
        html = QString("<span style='color:%1'>%2</span>").arg(colour, text.toHtmlEscaped());
    output_->appendHtml(html);
    auto* sb = output_->verticalScrollBar();
    sb->setValue(sb->maximum());
}

} // namespace simall::gui
