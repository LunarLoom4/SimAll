// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/DiagnosticsPanel.hpp
//
// Live readout of process/system resource use: CPU%, RSS, GPU mem, thread
// count, open file handles, MPI rank info.  Polls `core::ResourceMonitor`
// on a QTimer and renders into a sparse 2-column grid plus a small history
// sparkline for CPU% and RSS.
// =============================================================================
#pragma once

#include <QWidget>
#include <deque>

class QTimer;
class QLabel;
class QFormLayout;

namespace simall::gui {

class DiagnosticsPanel : public QWidget {
    Q_OBJECT
public:
    explicit DiagnosticsPanel(QWidget* parent = nullptr);

    void set_poll_interval_ms(int ms);

protected:
    void paintEvent(QPaintEvent*) override;

private slots:
    void poll();

private:
    QFormLayout*        form_         = nullptr;
    QLabel*             cpuLabel_     = nullptr;
    QLabel*             rssLabel_     = nullptr;
    QLabel*             threadsLabel_ = nullptr;
    QLabel*             handlesLabel_ = nullptr;
    QLabel*             gpuLabel_     = nullptr;
    QLabel*             mpiLabel_     = nullptr;
    QTimer*             timer_        = nullptr;

    std::deque<double>  cpuHistory_;
    std::deque<double>  rssHistory_;
    static constexpr int kHistory = 120;     // ~2 min at 1 Hz
};

}  // namespace simall::gui
