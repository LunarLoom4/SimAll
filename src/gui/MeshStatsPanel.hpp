// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/MeshStatsPanel.hpp
//
// Tabular mesh statistics: cell/face/node counts per zone, min/avg/max
// quality (skewness, aspect ratio, orthogonality), and a histogram of one
// chosen metric.  Population is fed by the caller via `set_stats()`.
// =============================================================================
#pragma once

#include <QWidget>
#include <vector>

class QTableWidget;
class QComboBox;
class QLabel;

namespace simall::gui {

struct MeshZoneStat {
    QString zoneName;
    int     cells   = 0;
    int     faces   = 0;
    int     nodes   = 0;
    double  minSkew = 0.0;
    double  avgSkew = 0.0;
    double  maxSkew = 0.0;
    double  minAspect = 0.0;
    double  avgAspect = 0.0;
    double  maxAspect = 0.0;
    double  minOrtho  = 0.0;
    double  avgOrtho  = 0.0;
    double  maxOrtho  = 0.0;
};

class MeshStatsPanel : public QWidget {
    Q_OBJECT
public:
    explicit MeshStatsPanel(QWidget* parent = nullptr);

    void set_stats(const std::vector<MeshZoneStat>& zones);
    void set_histogram(const std::vector<double>& bins, double xMin, double xMax);
    void clear_stats();

protected:
    void paintEvent(QPaintEvent*) override;

private:
    QTableWidget*           table_       = nullptr;
    QComboBox*              metricCombo_ = nullptr;
    QLabel*                 totalLabel_  = nullptr;
    std::vector<double>     histBins_;
    double                  histMin_ = 0.0, histMax_ = 1.0;
};

}  // namespace simall::gui
