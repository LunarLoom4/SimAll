// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/ResidualPlot.hpp
// Phase  : 2.1 (Dock panels → Residual Plot)
//
// Lightweight residual viewer: subscribes to core::events::SolverIteration
// and plots log-scale residuals. Pure QPainter, zero external chart dep.
// =============================================================================
#pragma once
#include <deque>
#include <QWidget>

namespace simall::gui
{

class ResidualPlot : public QWidget
{
    Q_OBJECT
public:
    explicit ResidualPlot(QWidget* parent = nullptr);
    void append(int iteration, double residualMax);
    void clear();

protected:
    void paintEvent(QPaintEvent*) override;

private:
    std::deque<std::pair<int, double>> history_;
};

} // namespace simall::gui
