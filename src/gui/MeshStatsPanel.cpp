// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/MeshStatsPanel.cpp
// =============================================================================
#include "gui/MeshStatsPanel.hpp"

#include <QTableWidget>

#include <QComboBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QPainter>
#include <QPaintEvent>
#include <QVBoxLayout>

namespace simall::gui
{

MeshStatsPanel::MeshStatsPanel(QWidget* parent) : QWidget(parent)
{
    auto* root = new QVBoxLayout(this);
    root->setContentsMargins(4, 4, 4, 4);
    root->setSpacing(4);

    auto* row = new QHBoxLayout();
    row->addWidget(new QLabel("Histogram metric:"));
    metricCombo_ = new QComboBox(this);
    metricCombo_->addItems({"Skewness", "Aspect ratio", "Orthogonality"});
    row->addWidget(metricCombo_);
    row->addStretch(1);
    totalLabel_ = new QLabel("0 cells / 0 faces / 0 nodes");
    row->addWidget(totalLabel_);
    root->addLayout(row);

    table_ = new QTableWidget(this);
    table_->setColumnCount(11);
    table_->setHorizontalHeaderLabels({"Zone",
                                       "Cells",
                                       "Faces",
                                       "Nodes",
                                       "Skew min",
                                       "Skew avg",
                                       "Skew max",
                                       "Asp avg",
                                       "Asp max",
                                       "Ortho min",
                                       "Ortho avg"});
    table_->horizontalHeader()->setStretchLastSection(true);
    table_->setEditTriggers(QAbstractItemView::NoEditTriggers);
    table_->setSelectionBehavior(QAbstractItemView::SelectRows);
    table_->verticalHeader()->setVisible(false);
    root->addWidget(table_, /*stretch*/ 1);

    setMinimumHeight(220);
}

void MeshStatsPanel::set_stats(const std::vector<MeshZoneStat>& zones)
{
    table_->setRowCount(int(zones.size()));
    int totC = 0, totF = 0, totN = 0;
    for (int i = 0; i < int(zones.size()); ++i) {
        const auto& z = zones[i];
        totC += z.cells;
        totF += z.faces;
        totN += z.nodes;
        auto setCell = [&](int col, const QString& txt) {
            auto* it = new QTableWidgetItem(txt);
            it->setTextAlignment(col == 0 ? Qt::AlignLeft | Qt::AlignVCenter
                                          : Qt::AlignRight | Qt::AlignVCenter);
            table_->setItem(i, col, it);
        };
        setCell(0, z.zoneName);
        setCell(1, QString::number(z.cells));
        setCell(2, QString::number(z.faces));
        setCell(3, QString::number(z.nodes));
        setCell(4, QString::number(z.minSkew, 'f', 3));
        setCell(5, QString::number(z.avgSkew, 'f', 3));
        setCell(6, QString::number(z.maxSkew, 'f', 3));
        setCell(7, QString::number(z.avgAspect, 'f', 3));
        setCell(8, QString::number(z.maxAspect, 'f', 3));
        setCell(9, QString::number(z.minOrtho, 'f', 3));
        setCell(10, QString::number(z.avgOrtho, 'f', 3));
    }
    totalLabel_->setText(QString("%1 cells / %2 faces / %3 nodes").arg(totC).arg(totF).arg(totN));
    table_->resizeColumnsToContents();
}

void MeshStatsPanel::set_histogram(const std::vector<double>& bins, double xMin, double xMax)
{
    histBins_ = bins;
    histMin_ = xMin;
    histMax_ = xMax;
    update();
}

void MeshStatsPanel::clear_stats()
{
    table_->setRowCount(0);
    histBins_.clear();
    totalLabel_->setText("0 cells / 0 faces / 0 nodes");
    update();
}

void MeshStatsPanel::paintEvent(QPaintEvent* ev)
{
    QWidget::paintEvent(ev);
    if (histBins_.empty())
        return;
    QPainter p(this);
    const int margin = 6;
    const int h = 60;
    const int top = height() - h - margin;
    const int w = width() - 2 * margin;
    p.setPen(QColor(0x55, 0x55, 0x55));
    p.drawRect(margin, top, w, h);
    const double bw = double(w) / histBins_.size();
    double mx = 0.0;
    for (double v : histBins_)
        mx = std::max(mx, v);
    if (mx <= 0.0)
        return;
    p.setBrush(QColor(0x00, 0x7A, 0xCC));
    p.setPen(Qt::NoPen);
    for (size_t i = 0; i < histBins_.size(); ++i) {
        const int bh = int((histBins_[i] / mx) * (h - 2));
        p.drawRect(int(margin + i * bw), top + h - bh, int(std::max(1.0, bw - 1.0)), bh);
    }
}

} // namespace simall::gui
