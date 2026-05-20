// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/SolverMonitorPanel.cpp
// =============================================================================
#include "gui/SolverMonitorPanel.hpp"

#include <QTableWidget>

#include <cmath>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPainter>
#include <QPaintEvent>
#include <QPushButton>
#include <QSplitter>
#include <QVBoxLayout>

namespace simall::gui
{

namespace
{
QString fmt_time(double s)
{
    if (!std::isfinite(s) || s < 0)
        return "—";
    const int total = int(s);
    const int hh = total / 3600;
    const int mm = (total / 60) % 60;
    const int ss = total % 60;
    return QString("%1:%2:%3")
        .arg(hh, 2, 10, QChar('0'))
        .arg(mm, 2, 10, QChar('0'))
        .arg(ss, 2, 10, QChar('0'));
}
} // namespace

SolverMonitorPanel::SolverMonitorPanel(QWidget* parent) : QWidget(parent)
{
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(4, 4, 4, 4);
    root->setSpacing(4);

    auto* head = new QHBoxLayout();
    stateLbl_ = new QLabel("State: Idle", this);
    etaLbl_ = new QLabel("ETA: —", this);
    btnPause_ = new QPushButton("Pause", this);
    btnStop_ = new QPushButton("Stop", this);
    btnSnap_ = new QPushButton("Snapshot", this);
    head->addWidget(stateLbl_);
    head->addSpacing(12);
    head->addWidget(etaLbl_);
    head->addStretch(1);
    head->addWidget(btnPause_);
    head->addWidget(btnStop_);
    head->addWidget(btnSnap_);
    root->addLayout(head);

    summary_ = new QTableWidget(0, 4, this);
    summary_->setHorizontalHeaderLabels({"Variable", "Current", "Initial", "Drop"});
    summary_->horizontalHeader()->setStretchLastSection(true);
    summary_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    summary_->verticalHeader()->setVisible(false);
    summary_->setMaximumHeight(140);
    root->addWidget(summary_);

    setMinimumHeight(280);

    connect(btnPause_, &QPushButton::clicked, this, &SolverMonitorPanel::pauseRequested);
    connect(btnStop_, &QPushButton::clicked, this, &SolverMonitorPanel::stopRequested);
    connect(btnSnap_, &QPushButton::clicked, this, &SolverMonitorPanel::snapshotRequested);
}

void SolverMonitorPanel::register_variable(const QString& name, const QColor& c)
{
    if (residuals_.find(name) == residuals_.end()) {
        residuals_[name] = {c, {}};
        ordering_.push_back(name);
        rebuild_summary();
    } else {
        residuals_[name].colour = c;
    }
}

void SolverMonitorPanel::append_residual(const QString& name, int it, double r)
{
    auto& s = residuals_[name];
    if (s.colour.alpha() == 0)
        s.colour = QColor(0xE0, 0xE0, 0xE0);
    const double logr = (r > 0.0) ? std::log10(r) : -30.0;
    s.samples.push_back({double(it), logr});
    while ((int) s.samples.size() > kMaxSamples)
        s.samples.pop_front();
    if (std::find(ordering_.begin(), ordering_.end(), name) == ordering_.end()) {
        ordering_.push_back(name);
        rebuild_summary();
    }
    update();
}

void SolverMonitorPanel::append_imbalance(const QString& k, double pct)
{
    imbalances_[k] = pct;
    rebuild_summary();
}

void SolverMonitorPanel::set_run_state(const QString& s)
{
    state_ = s;
    stateLbl_->setText("State: " + s);
}

void SolverMonitorPanel::set_eta(double s)
{
    eta_ = s;
    etaLbl_->setText("ETA: " + fmt_time(s));
}

void SolverMonitorPanel::set_iteration(int it, double sim, double wall)
{
    iter_ = it;
    simTime_ = sim;
    wallTime_ = wall;
}

void SolverMonitorPanel::clear_history()
{
    for (auto& [_, s] : residuals_)
        s.samples.clear();
    imbalances_.clear();
    summary_->setRowCount(0);
    update();
}

void SolverMonitorPanel::rebuild_summary()
{
    summary_->setRowCount(int(ordering_.size()));
    for (int i = 0; i < int(ordering_.size()); ++i) {
        const auto& name = ordering_[i];
        const auto& s = residuals_[name];
        const double cur = s.samples.empty() ? 0.0 : std::pow(10.0, s.samples.back().y());
        const double init = s.samples.empty() ? 0.0 : std::pow(10.0, s.samples.front().y());
        const double drop = (init > 0 && cur > 0) ? init / cur : 0.0;
        summary_->setItem(i, 0, new QTableWidgetItem(name));
        summary_->setItem(i, 1, new QTableWidgetItem(QString::number(cur, 'e', 3)));
        summary_->setItem(i, 2, new QTableWidgetItem(QString::number(init, 'e', 3)));
        summary_->setItem(i, 3, new QTableWidgetItem(QString::number(drop, 'f', 2)));
        for (int c = 0; c < 4; ++c) {
            if (auto* it = summary_->item(i, c)) {
                if (c == 0)
                    it->setForeground(s.colour);
                it->setFlags(it->flags() & ~Qt::ItemIsEditable);
            }
        }
    }
}

void SolverMonitorPanel::paintEvent(QPaintEvent* ev)
{
    QWidget::paintEvent(ev);
    if (residuals_.empty())
        return;
    QPainter p(this);
    const int margin = 8;
    const int plotTop = 60 + (summary_ ? summary_->height() : 0);
    const int plotBottom = height() - margin - 24;
    if (plotBottom <= plotTop + 20)
        return;
    const int plotLeft = margin + 40;
    const int plotRight = width() - margin;

    // Determine axis ranges.
    double xMin = 1e300, xMax = -1e300, yMin = 1e300, yMax = -1e300;
    for (auto& [_, s] : residuals_) {
        for (const auto& pt : s.samples) {
            xMin = std::min(xMin, pt.x());
            xMax = std::max(xMax, pt.x());
            yMin = std::min(yMin, pt.y());
            yMax = std::max(yMax, pt.y());
        }
    }
    if (xMin >= xMax)
        return;
    if (yMin >= yMax) {
        yMin -= 1.0;
        yMax += 1.0;
    }

    p.setPen(QColor(0x55, 0x55, 0x55));
    p.drawRect(plotLeft, plotTop, plotRight - plotLeft, plotBottom - plotTop);

    auto X = [&](double x) {
        return plotLeft + int((x - xMin) / (xMax - xMin) * (plotRight - plotLeft));
    };
    auto Y = [&](double y) {
        return plotBottom - int((y - yMin) / (yMax - yMin) * (plotBottom - plotTop));
    };

    p.setPen(QColor(0x80, 0x80, 0x80));
    for (int g = 0; g <= 4; ++g) {
        const double ly = yMin + (yMax - yMin) * g / 4.0;
        const int yy = Y(ly);
        p.drawLine(plotLeft, yy, plotRight, yy);
        p.drawText(margin, yy + 4, QString::number(std::pow(10.0, ly), 'e', 1));
    }

    for (auto& [_, s] : residuals_) {
        if (s.samples.size() < 2)
            continue;
        QPen pen(s.colour);
        pen.setWidth(2);
        p.setPen(pen);
        for (size_t i = 1; i < s.samples.size(); ++i)
            p.drawLine(X(s.samples[i - 1].x()),
                       Y(s.samples[i - 1].y()),
                       X(s.samples[i].x()),
                       Y(s.samples[i].y()));
    }

    p.setPen(QColor(0xAA, 0xAA, 0xAA));
    p.drawText(plotLeft,
               plotBottom + 16,
               QString("iter %1   sim %2 s   wall %3   ETA %4")
                   .arg(iter_)
                   .arg(simTime_, 0, 'f', 4)
                   .arg(fmt_time(wallTime_))
                   .arg(fmt_time(eta_)));
}

} // namespace simall::gui
