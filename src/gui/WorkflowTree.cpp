#include "gui/WorkflowTree.hpp"

namespace simall::gui
{

WorkflowTree::WorkflowTree(QWidget* p) : QTreeWidget(p)
{
    setHeaderHidden(true);
    setIndentation(14);
    setAnimated(true);
    setMinimumWidth(260); // Section 3.4
    setUniformRowHeights(true);
}

void WorkflowTree::rebuild_from(io::ProjectNode* root)
{
    clear();
    if (!root)
        return;
    auto* item = new QTreeWidgetItem(this, {QString::fromStdString(root->name())});
    item->setForeground(0, state_color(root->state()));
    for (const auto& c : root->children())
        add_node(item, c.get());
    expandAll();
}

void WorkflowTree::add_node(QTreeWidgetItem* parent, io::ProjectNode* node)
{
    auto* item = new QTreeWidgetItem(parent, {QString::fromStdString(node->name())});
    item->setForeground(0, state_color(node->state()));
    for (const auto& c : node->children())
        add_node(item, c.get());
}

QColor WorkflowTree::state_color(io::NodeState s) const
{
    switch (s) {
    case io::NodeState::Empty:
        return QColor(0x80, 0x80, 0x80);
    case io::NodeState::Valid:
        return QColor(0x6A, 0xB0, 0x4A);
    case io::NodeState::Stale:
        return QColor(0xD7, 0xBA, 0x7D);
    case io::NodeState::Computing:
        return QColor(0x00, 0x7A, 0xCC);
    case io::NodeState::Failed:
        return QColor(0xF4, 0x47, 0x47);
    }
    return Qt::white;
}

} // namespace simall::gui
