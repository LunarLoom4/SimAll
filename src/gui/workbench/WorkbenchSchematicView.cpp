// =============================================================================
// SimAll Beta — GUI / Workbench panel
// File   : src/gui/workbench/WorkbenchSchematicView.cpp
// Phase  : 22 Pass 22.4
// =============================================================================
#include "gui/workbench/WorkbenchSchematicView.hpp"

#include "workbench/Cell.hpp"
#include "workbench/CellPort.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/WorkflowEngine.hpp"

#include <algorithm>
#include <QAction>
#include <QBrush>
#include <QColor>
#include <QContextMenuEvent>
#include <QFont>
#include <QFontMetrics>
#include <QGraphicsItem>
#include <QGraphicsLineItem>
#include <QGraphicsPathItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsSimpleTextItem>
#include <QInputDialog>
#include <QKeyEvent>
#include <QLineF>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <unordered_set>

namespace simall::gui::workbench
{

using simall::workbench::Cell;
using simall::workbench::CellId;
using simall::workbench::CellKind;
using simall::workbench::CellLink;
using simall::workbench::CellPort;
using simall::workbench::CellState;
using simall::workbench::kInvalidCellId;
using simall::workbench::kInvalidPortId;
using simall::workbench::PortDirection;
using simall::workbench::PortId;
using simall::workbench::Schematic;
using simall::workbench::StateMachine;
using simall::workbench::WorkflowEngine;

namespace
{

// ---- visual constants (kept local so theming can swap them out later) ----
constexpr qreal kCellW = 180.0;
constexpr qreal kCellH = 92.0;
constexpr qreal kHeaderH = 22.0;
constexpr qreal kPortR = 5.0;
constexpr qreal kColumnGap = 240.0;
constexpr qreal kRowGap = 140.0;

QColor state_color(CellState s)
{
    switch (s) {
    case CellState::UpToDate:
        return QColor(76, 175, 80); // green
    case CellState::RefreshRequired:
        return QColor(255, 193, 7); // amber
    case CellState::Unfulfilled:
        return QColor(120, 144, 156); // slate
    case CellState::Failed:
        return QColor(244, 67, 54); // red
    }
    return QColor(120, 144, 156);
}

QString kind_initial(CellKind k)
{
    switch (k) {
    case CellKind::Geometry:
        return "G";
    case CellKind::Mesh:
        return "M";
    case CellKind::Setup:
        return "S";
    case CellKind::Solution:
        return "R"; // Run
    case CellKind::Results:
        return "V"; // View
    case CellKind::Custom:
        return "C";
    }
    return "?";
}

// Stringify CellState for status tooltips / menu enables.
QString state_text(CellState s)
{
    switch (s) {
    case CellState::UpToDate:
        return "Up-to-date";
    case CellState::RefreshRequired:
        return "Refresh required";
    case CellState::Unfulfilled:
        return "Unfulfilled";
    case CellState::Failed:
        return "Failed";
    }
    return "?";
}

} // namespace


// ===========================================================================
// CellNodeItem -- one rounded card per Cell.  Holds child PortHandleItems
// laid out on the left / right edges so hit-testing is trivial.
// ===========================================================================
class CellNodeItem final : public QGraphicsItem
{
public:
    CellNodeItem(
        CellId id, CellKind kind, QString label, CellState state, std::vector<CellPort> ports)
        : id_(id), kind_(kind), label_(std::move(label)), state_(state), ports_(std::move(ports))
    {
        setFlag(ItemIsSelectable, true);
        setFlag(ItemIsMovable, true);
        setAcceptHoverEvents(false);
        setData(0, QVariant::fromValue<uint>(id_));
    }

    [[nodiscard]] QRectF boundingRect() const override
    {
        return QRectF(-kPortR, 0, kCellW + 2 * kPortR, kCellH);
    }

    void paint(QPainter* p, const QStyleOptionGraphicsItem*, QWidget*) override
    {
        const QRectF body(0, 0, kCellW, kCellH);
        QPainterPath frame;
        frame.addRoundedRect(body, 8, 8);

        // Body
        p->setRenderHint(QPainter::Antialiasing, true);
        p->fillPath(frame, QColor(40, 44, 52));

        // Header strip in state colour.  We fill a plain rectangle on
        // top of the rounded body; the body's corner radius hides the
        // top corners, leaving a flat bottom edge against the body --
        // exactly the look the Workbench shipped product uses.
        p->fillRect(QRectF(0, 0, kCellW, kHeaderH), state_color(state_));

        // Outline (selected -> bright blue)
        QPen pen(isSelected() ? QColor(33, 150, 243) : QColor(70, 75, 85), 1.5);
        p->setPen(pen);
        p->drawPath(frame);

        // Kind initial badge
        p->setPen(Qt::white);
        QFont f = p->font();
        f.setBold(true);
        f.setPointSize(10);
        p->setFont(f);
        p->drawText(QRectF(4, 0, 20, kHeaderH), Qt::AlignCenter, kind_initial(kind_));

        // Label
        p->setPen(QColor(200, 210, 220));
        QFont lf = p->font();
        lf.setBold(false);
        lf.setPointSize(9);
        p->setFont(lf);
        p->drawText(QRectF(28, 0, kCellW - 32, kHeaderH), Qt::AlignLeft | Qt::AlignVCenter, label_);

        // State text + small status circle bottom-right
        p->setPen(QColor(150, 160, 170));
        QFont sf = p->font();
        sf.setPointSize(8);
        p->setFont(sf);
        p->drawText(QRectF(8, kCellH - 16, kCellW - 16, 14),
                    Qt::AlignRight | Qt::AlignVCenter,
                    state_text(state_));

        // Port dots
        int in_idx = 0, out_idx = 0;
        const int n_in = num_ports(PortDirection::Input);
        const int n_out = num_ports(PortDirection::Output);
        for (const CellPort& port : ports_) {
            const QPointF pos = port_local_pos(port, in_idx, out_idx, n_in, n_out);
            if (port.direction == PortDirection::Input)
                ++in_idx;
            else
                ++out_idx;
            p->setPen(QPen(QColor(220, 220, 220), 1.0));
            p->setBrush(port.direction == PortDirection::Input ? QColor(33, 150, 243)
                                                               : QColor(255, 152, 0));
            p->drawEllipse(pos, kPortR, kPortR);
        }
    }

    [[nodiscard]] CellId id() const noexcept { return id_; }
    [[nodiscard]] CellKind kind() const noexcept { return kind_; }
    [[nodiscard]] CellState state() const noexcept { return state_; }
    [[nodiscard]] const QString& label() const noexcept { return label_; }
    [[nodiscard]] const std::vector<CellPort>& ports() const noexcept { return ports_; }

    // Scene-space position of a given port (used by LinkEdgeItem and by
    // drag-link hit-testing).  Returns top-left+offset; the disc is drawn
    // centred on that point.
    [[nodiscard]] QPointF port_scene_pos(PortId pid) const
    {
        int in_idx = 0, out_idx = 0;
        const int n_in = num_ports(PortDirection::Input);
        const int n_out = num_ports(PortDirection::Output);
        for (const CellPort& port : ports_) {
            const QPointF local = port_local_pos(port, in_idx, out_idx, n_in, n_out);
            if (port.direction == PortDirection::Input)
                ++in_idx;
            else
                ++out_idx;
            if (port.id == pid)
                return mapToScene(local);
        }
        return mapToScene(QPointF(0, 0));
    }

private:
    int num_ports(PortDirection d) const
    {
        int n = 0;
        for (const CellPort& p : ports_)
            if (p.direction == d)
                ++n;
        return n;
    }

    static QPointF port_local_pos(
        const CellPort& port, int in_idx, int out_idx, int n_in, int n_out)
    {
        const qreal usable = kCellH - kHeaderH - 12.0;
        if (port.direction == PortDirection::Input) {
            const qreal y = kHeaderH + 6.0 + usable * (in_idx + 0.5) / std::max(1, n_in);
            return QPointF(0.0, y);
        } else {
            const qreal y = kHeaderH + 6.0 + usable * (out_idx + 0.5) / std::max(1, n_out);
            return QPointF(kCellW, y);
        }
    }

    CellId id_;
    CellKind kind_;
    QString label_;
    CellState state_;
    std::vector<CellPort> ports_;
};


// ===========================================================================
// LinkEdgeItem -- cubic bezier between two scene points; updates on demand.
// ===========================================================================
class LinkEdgeItem final : public QGraphicsPathItem
{
public:
    LinkEdgeItem(CellLink link, QPointF from, QPointF to) : link_(link)
    {
        setPen(QPen(QColor(180, 190, 200), 1.6));
        update_geometry(from, to);
        setData(0, QStringLiteral("link"));
    }

    void update_geometry(QPointF a, QPointF b)
    {
        QPainterPath p(a);
        const qreal dx = std::max<qreal>(40.0, std::abs(b.x() - a.x()) * 0.5);
        p.cubicTo(QPointF(a.x() + dx, a.y()), QPointF(b.x() - dx, b.y()), b);
        setPath(p);
    }

    [[nodiscard]] CellLink link() const noexcept { return link_; }

private:
    CellLink link_;
};


// ===========================================================================
// WorkbenchSchematicView
// ===========================================================================
WorkbenchSchematicView::WorkbenchSchematicView(QWidget* parent)
    : QGraphicsView(parent), scene_(new QGraphicsScene(this))
{
    setScene(scene_);
    setRenderHint(QPainter::Antialiasing, true);
    setDragMode(QGraphicsView::RubberBandDrag);
    setBackgroundBrush(QColor(28, 30, 34));
    setMinimumSize(420, 320);
}

WorkbenchSchematicView::~WorkbenchSchematicView() = default;

void WorkbenchSchematicView::set_sources(Schematic* s, StateMachine* sm, WorkflowEngine* eng)
{
    schematic_ = s;
    state_ = sm;
    engine_ = eng;
    refresh();
}

void WorkbenchSchematicView::refresh()
{
    rebuild_scene();
}

void WorkbenchSchematicView::rebuild_scene()
{
    scene_->clear();
    cell_items_.clear();
    drag_preview_ = nullptr;
    dragging_link_ = false;
    if (!schematic_)
        return;

    // Create cell items
    for (const Cell& c : schematic_->cells()) {
        auto* item = new CellNodeItem(
            c.id(), c.kind(), QString::fromStdString(std::string(c.label())), c.state(), c.ports());
        scene_->addItem(item);
        cell_items_[c.id()] = item;
    }

    layout_cells();

    // Create link items
    for (const CellLink& l : schematic_->links()) {
        auto fi = cell_items_.find(l.from_cell);
        auto ti = cell_items_.find(l.to_cell);
        if (fi == cell_items_.end() || ti == cell_items_.end())
            continue;
        const QPointF a = fi->second->port_scene_pos(l.from_port);
        const QPointF b = ti->second->port_scene_pos(l.to_port);
        scene_->addItem(new LinkEdgeItem(l, a, b));
    }

    scene_->setSceneRect(scene_->itemsBoundingRect().adjusted(-60, -60, 60, 60));
}

void WorkbenchSchematicView::layout_cells()
{
    if (!schematic_)
        return;

    // Column = depth in topological order (longest path from a root).
    std::unordered_map<CellId, int> depth;
    for (const CellId id : schematic_->topological_order()) {
        int d = 0;
        for (const CellId up : schematic_->upstream(id)) {
            auto it = depth.find(up);
            if (it != depth.end())
                d = std::max(d, it->second + 1);
        }
        depth[id] = d;
    }

    // Bucket cells by column, sort each column by CellId so layout is
    // deterministic across runs.
    std::unordered_map<int, std::vector<CellId>> by_col;
    for (const auto& [id, d] : depth)
        by_col[d].push_back(id);
    for (auto& [_, v] : by_col)
        std::sort(v.begin(), v.end());

    for (const auto& [col, ids] : by_col) {
        for (std::size_t row = 0; row < ids.size(); ++row) {
            auto it = cell_items_.find(ids[row]);
            if (it == cell_items_.end())
                continue;
            it->second->setPos(col * kColumnGap, static_cast<qreal>(row) * kRowGap);
        }
    }
}


// ---------------------------------------------------------------------------
// Hit-testing helpers
// ---------------------------------------------------------------------------
std::optional<WorkbenchSchematicView::PortHit> WorkbenchSchematicView::hit_port(
    const QPointF& scenePos) const
{
    constexpr qreal kPortPickR = 8.0;
    for (const auto& [id, item] : cell_items_) {
        for (const CellPort& p : item->ports()) {
            const QPointF c = item->port_scene_pos(p.id);
            if (QLineF(c, scenePos).length() <= kPortPickR) {
                return PortHit{id, p.id, p.direction, c};
            }
        }
    }
    return std::nullopt;
}

CellNodeItem* WorkbenchSchematicView::hit_cell(const QPointF& scenePos) const
{
    for (QGraphicsItem* it : scene_->items(scenePos)) {
        if (auto* c = dynamic_cast<CellNodeItem*>(it))
            return c;
    }
    return nullptr;
}


// ---------------------------------------------------------------------------
// Drag-link interaction
// ---------------------------------------------------------------------------
void WorkbenchSchematicView::mousePressEvent(QMouseEvent* ev)
{
    const QPointF s = mapToScene(ev->pos());
    if (ev->button() == Qt::LeftButton) {
        if (auto hit = hit_port(s); hit && hit->direction == PortDirection::Output) {
            dragging_link_ = true;
            drag_origin_ = *hit;
            drag_preview_ = scene_->addLine(QLineF(hit->scene_pos, s),
                                            QPen(QColor(33, 150, 243), 1.4, Qt::DashLine));
            ev->accept();
            return;
        }
    }
    QGraphicsView::mousePressEvent(ev);
}

void WorkbenchSchematicView::mouseMoveEvent(QMouseEvent* ev)
{
    if (dragging_link_ && drag_preview_) {
        drag_preview_->setLine(QLineF(drag_origin_.scene_pos, mapToScene(ev->pos())));
        ev->accept();
        return;
    }
    QGraphicsView::mouseMoveEvent(ev);
}

void WorkbenchSchematicView::mouseReleaseEvent(QMouseEvent* ev)
{
    if (dragging_link_) {
        const QPointF s = mapToScene(ev->pos());
        if (auto hit = hit_port(s); hit && hit->direction == PortDirection::Input) {
            emit connectRequested(
                CellLink{drag_origin_.cell, drag_origin_.port, hit->cell, hit->port});
        }
        if (drag_preview_) {
            scene_->removeItem(drag_preview_);
            delete drag_preview_;
        }
        drag_preview_ = nullptr;
        dragging_link_ = false;
        ev->accept();
        return;
    }
    QGraphicsView::mouseReleaseEvent(ev);
}


// ---------------------------------------------------------------------------
// Context menu
// ---------------------------------------------------------------------------
void WorkbenchSchematicView::contextMenuEvent(QContextMenuEvent* ev)
{
    if (!schematic_)
        return;
    const QPointF s = mapToScene(ev->pos());

    // 1) Link under cursor?
    for (QGraphicsItem* it : scene_->items(s)) {
        if (auto* e = dynamic_cast<LinkEdgeItem*>(it)) {
            QMenu menu(this);
            QAction* del = menu.addAction(tr("Disconnect"));
            if (menu.exec(ev->globalPos()) == del) {
                emit disconnectRequested(e->link());
            }
            return;
        }
    }

    // 2) Cell under cursor?
    if (CellNodeItem* c = hit_cell(s)) {
        QMenu menu(this);
        QAction* aRun = menu.addAction(tr("Run cell"));
        QAction* aBind = menu.addAction(tr("Bind adapter..."));
        menu.addSeparator();
        QAction* aRename = menu.addAction(tr("Rename..."));
        QAction* aRemove = menu.addAction(tr("Remove cell"));
        menu.addSeparator();
        QMenu* state_menu = menu.addMenu(tr("Mark state"));
        QAction* sUp = state_menu->addAction(tr("Up-to-date"));
        QAction* sRR = state_menu->addAction(tr("Refresh required"));
        QAction* sUnf = state_menu->addAction(tr("Unfulfilled"));
        QAction* sFail = state_menu->addAction(tr("Failed"));

        QAction* picked = menu.exec(ev->globalPos());
        if (!picked)
            return;
        if (picked == aRun) {
            emit cellRunRequested(c->id());
        } else if (picked == aBind) {
            emit cellBindAdapterRequested(c->id());
        } else if (picked == aRename) {
            bool ok = false;
            const QString fresh = QInputDialog::getText(
                this, tr("Rename cell"), tr("Label:"), QLineEdit::Normal, c->label(), &ok);
            if (ok && !fresh.isEmpty())
                emit cellRenameRequested(c->id(), fresh);
        } else if (picked == aRemove) {
            emit cellRemoveRequested(c->id());
        } else if (picked == sUp)
            emit cellStateChangeRequested(c->id(), CellState::UpToDate);
        else if (picked == sRR)
            emit cellStateChangeRequested(c->id(), CellState::RefreshRequired);
        else if (picked == sUnf)
            emit cellStateChangeRequested(c->id(), CellState::Unfulfilled);
        else if (picked == sFail)
            emit cellStateChangeRequested(c->id(), CellState::Failed);
        return;
    }

    // 3) Empty canvas -> Add cell submenu by kind.
    QMenu menu(this);
    QMenu* add_menu = menu.addMenu(tr("Add cell"));
    QAction* aG = add_menu->addAction(tr("Geometry"));
    QAction* aM = add_menu->addAction(tr("Mesh"));
    QAction* aS = add_menu->addAction(tr("Setup"));
    QAction* aR = add_menu->addAction(tr("Solution"));
    QAction* aV = add_menu->addAction(tr("Results"));
    QAction* aC = add_menu->addAction(tr("Custom"));
    QAction* picked = menu.exec(ev->globalPos());
    if (picked == aG)
        emit cellAddRequested(CellKind::Geometry);
    else if (picked == aM)
        emit cellAddRequested(CellKind::Mesh);
    else if (picked == aS)
        emit cellAddRequested(CellKind::Setup);
    else if (picked == aR)
        emit cellAddRequested(CellKind::Solution);
    else if (picked == aV)
        emit cellAddRequested(CellKind::Results);
    else if (picked == aC)
        emit cellAddRequested(CellKind::Custom);
}


// ---------------------------------------------------------------------------
// Keyboard: Delete removes the selected cell or link, F2 renames.
// ---------------------------------------------------------------------------
void WorkbenchSchematicView::keyPressEvent(QKeyEvent* ev)
{
    const auto selected = scene_->selectedItems();
    if (ev->key() == Qt::Key_Delete && !selected.isEmpty()) {
        for (QGraphicsItem* it : selected) {
            if (auto* c = dynamic_cast<CellNodeItem*>(it)) {
                emit cellRemoveRequested(c->id());
            } else if (auto* e = dynamic_cast<LinkEdgeItem*>(it)) {
                emit disconnectRequested(e->link());
            }
        }
        ev->accept();
        return;
    }
    if (ev->key() == Qt::Key_F2 && selected.size() == 1) {
        if (auto* c = dynamic_cast<CellNodeItem*>(selected.front())) {
            bool ok = false;
            const QString fresh = QInputDialog::getText(
                this, tr("Rename cell"), tr("Label:"), QLineEdit::Normal, c->label(), &ok);
            if (ok && !fresh.isEmpty())
                emit cellRenameRequested(c->id(), fresh);
            ev->accept();
            return;
        }
    }
    QGraphicsView::keyPressEvent(ev);
}

} // namespace simall::gui::workbench
