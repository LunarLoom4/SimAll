// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/WorkflowTree.hpp
// Phase  : 2.1 / Section 3.9 — workflow tree mirrors io::ProjectNode graph.
// Tree-node colour reflects NodeState (Empty/Valid/Stale/Computing/Failed).
// =============================================================================
#pragma once
#include <QTreeWidget>
#include "io/ProjectNode.hpp"

namespace simall::gui {

class WorkflowTree : public QTreeWidget {
    Q_OBJECT
public:
    explicit WorkflowTree(QWidget* parent = nullptr);
    void rebuild_from(io::ProjectNode* root);
private:
    void add_node(QTreeWidgetItem* parent, io::ProjectNode* node);
    QColor state_color(io::NodeState s) const;
};

}  // namespace simall::gui
