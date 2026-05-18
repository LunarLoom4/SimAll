#include "gui/ResidualPlot.hpp"
#include <QPainter>
#include <QPaintEvent>
#include <algorithm>
#include <cmath>

namespace simall::gui {

ResidualPlot::ResidualPlot(QWidget* p) : QWidget(p) {
    setMinimumHeight(160);
    setAutoFillBackground(true);
    QPalette pal = palette();
    pal.setColor(QPalette::Window, QColor(0x1E, 0x1E, 0x1E));
    setPalette(pal);
}

void ResidualPlot::append(int it, double r) {
    history_.emplace_back(it, std::max(r, 1e-30));
    if (history_.size() > 4096) history_.pop_front();
    update();
}

void ResidualPlot::clear() { history_.clear(); update(); }

void ResidualPlot::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);
    const QRect r = rect().adjusted(36, 12, -12, -24);
    p.setPen(QColor(60, 60, 60));
    p.drawRect(r);
    if (history_.size() < 2) return;
    const double minL = -10.0, maxL = 2.0;
    auto xToPx = [&](int i){
        double t = double(i) / double(history_.back().first);
        return r.left() + t * r.width();
    };
    auto yToPx = [&](double v){
        double l = std::log10(v);
        double t = (l - minL) / (maxL - minL);
        t = std::clamp(t, 0.0, 1.0);
        return r.bottom() - t * r.height();
    };
    p.setPen(QPen(QColor(0x00, 0x7A, 0xCC), 1.5));
    for (std::size_t i = 1; i < history_.size(); ++i) {
        p.drawLine(QPointF(xToPx(history_[i-1].first), yToPx(history_[i-1].second)),
                   QPointF(xToPx(history_[i].first),   yToPx(history_[i].second)));
    }
    p.setPen(QColor(200, 200, 200));
    p.drawText(rect().adjusted(4, 4, -4, 0), Qt::AlignLeft|Qt::AlignTop, "Residuals (log10)");
}

}  // namespace simall::gui
