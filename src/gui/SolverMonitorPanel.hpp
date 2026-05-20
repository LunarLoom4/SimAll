// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/SolverMonitorPanel.hpp
//
// Run-time solver dashboard.
//   * Residuals plot for arbitrary number of variables (each colour-coded).
//   * Imbalance plot (mass, momentum, energy) bottom strip.
//   * Convergence summary: iteration, sim time, wall time, ETA, drop ratio.
//   * Controls: Pause / Stop / Snapshot. Buttons emit signals; the
//     application layer translates those into command bus calls.
//
// The plot uses Qt-native QPainter; no external charting dependency.
// =============================================================================
#pragma once

#include <deque>
#include <QWidget>
#include <unordered_map>
#include <vector>

class QTableWidget;
class QPushButton;
class QLabel;
class QSplitter;

namespace simall::gui
{

class SolverMonitorPanel : public QWidget
{
    Q_OBJECT
public:
    explicit SolverMonitorPanel(QWidget* parent = nullptr);

    void register_variable(const QString& name, const QColor& colour);
    void append_residual(const QString& name, int iteration, double residual);
    void append_imbalance(const QString& kind, double percent);
    void set_run_state(const QString& state); // "Running", "Paused", "Idle"
    void set_eta(double secondsRemaining);
    void set_iteration(int iter, double simTimeSec, double wallTimeSec);
    void clear_history();

signals:
    void pauseRequested();
    void stopRequested();
    void snapshotRequested();

protected:
    void paintEvent(QPaintEvent*) override;

private:
    void rebuild_summary();

    struct Series
    {
        QColor colour;
        std::deque<QPointF> samples; // {iteration, log10(residual)}
    };
    std::unordered_map<QString, Series> residuals_;
    std::unordered_map<QString, double> imbalances_;
    std::vector<QString> ordering_; // insertion order
    QString state_ = "Idle";
    int iter_ = 0;
    double simTime_ = 0.0;
    double wallTime_ = 0.0;
    double eta_ = -1.0;

    QTableWidget* summary_ = nullptr;
    QPushButton* btnPause_ = nullptr;
    QPushButton* btnStop_ = nullptr;
    QPushButton* btnSnap_ = nullptr;
    QLabel* stateLbl_ = nullptr;
    QLabel* etaLbl_ = nullptr;

    static constexpr int kMaxSamples = 8192;
};

} // namespace simall::gui
