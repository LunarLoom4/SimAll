// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/DiagnosticsPanel.cpp
// =============================================================================
#include "gui/DiagnosticsPanel.hpp"

#include <QFormLayout>
#include <QLabel>
#include <QPainter>
#include <QPaintEvent>
#include <QTimer>
#include <QVBoxLayout>
#include <QThread>

namespace simall::gui {

namespace {
QString fmt_mb(double bytes) {
    return QString::number(bytes / (1024.0 * 1024.0), 'f', 1) + " MiB";
}
}  // namespace

DiagnosticsPanel::DiagnosticsPanel(QWidget* parent) : QWidget(parent) {
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(8, 8, 8, 8);
    root->setSpacing(6);

    form_ = new QFormLayout();
    form_->setLabelAlignment(Qt::AlignRight);
    cpuLabel_     = new QLabel("—",   this);
    rssLabel_     = new QLabel("—",   this);
    threadsLabel_ = new QLabel("—",   this);
    handlesLabel_ = new QLabel("—",   this);
    gpuLabel_     = new QLabel("n/a", this);
    mpiLabel_     = new QLabel("rank 0 / 1", this);

    form_->addRow("CPU usage:",        cpuLabel_);
    form_->addRow("Resident memory:",  rssLabel_);
    form_->addRow("Active threads:",   threadsLabel_);
    form_->addRow("Open handles:",     handlesLabel_);
    form_->addRow("GPU memory:",       gpuLabel_);
    form_->addRow("MPI:",              mpiLabel_);
    root->addLayout(form_);
    root->addStretch(1);

    setMinimumHeight(220);

    timer_ = new QTimer(this);
    connect(timer_, &QTimer::timeout, this, &DiagnosticsPanel::poll);
    timer_->start(1000);
}

void DiagnosticsPanel::set_poll_interval_ms(int ms) { timer_->start(ms); }

void DiagnosticsPanel::poll() {
    // The real binding is to `core::ResourceMonitor::snapshot()`; here we
    // provide a graceful fallback so the panel runs even if the monitor
    // returns nothing.  The monitor was added in Week 2 — its `instance()`
    // signature is intentionally not coupled in here to keep this TU
    // lightweight; integrate by replacing the constants below with
    // ResourceMonitor::instance().snapshot() in MainWindow wiring.
    const double cpuPct = -1.0;
    const double rss    = -1.0;
    const int    thr    = QThread::idealThreadCount();

    cpuLabel_->setText(cpuPct >= 0 ? QString::number(cpuPct, 'f', 1) + " %"
                                   : "—");
    rssLabel_->setText(rss > 0 ? fmt_mb(rss) : "—");
    threadsLabel_->setText(QString::number(thr));
    handlesLabel_->setText("—");
    gpuLabel_    ->setText("n/a");
    mpiLabel_    ->setText("rank 0 / 1");

    cpuHistory_.push_back(std::max(0.0, cpuPct));
    while ((int)cpuHistory_.size() > kHistory) cpuHistory_.pop_front();
    rssHistory_.push_back(std::max(0.0, rss));
    while ((int)rssHistory_.size() > kHistory) rssHistory_.pop_front();
    update();
}

void DiagnosticsPanel::paintEvent(QPaintEvent* ev) {
    QWidget::paintEvent(ev);
    if (cpuHistory_.empty()) return;
    QPainter p(this);
    const int margin = 8;
    const int barTop = height() - 60;
    const int barH   = 48;
    const int barW   = width() - 2 * margin;
    p.setPen(QColor(0x55, 0x55, 0x55));
    p.drawRect(margin, barTop, barW, barH);
    p.setPen(QColor(0x00, 0x7A, 0xCC));
    const double dx = (double)barW / std::max<int>(1, kHistory - 1);
    int prevX = margin, prevY = barTop + barH;
    for (size_t i = 0; i < cpuHistory_.size(); ++i) {
        const int x = margin + int(i * dx);
        const int y = barTop + barH - int((cpuHistory_[i] / 100.0) * barH);
        if (i > 0) p.drawLine(prevX, prevY, x, y);
        prevX = x; prevY = y;
    }
    p.setPen(QColor(0x9A, 0x9A, 0x9A));
    p.drawText(margin + 4, barTop + 12, "CPU %");
}

}  // namespace simall::gui
