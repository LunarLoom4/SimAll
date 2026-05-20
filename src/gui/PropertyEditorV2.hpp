// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/PropertyEditorV2.hpp
//
// Metadata-driven editor that consumes `gui_core::PropertyDescriptor` and
// builds the correct widget per type, wires validation feedback, and
// re-evaluates visibility / enablement after every commit.
//
// Differences from the legacy `PropertyEditor`:
//   * Groups properties by `category` into collapsible sections.
//   * Supports Vec3, FilePath, Color, Range editor types.
//   * Shows inline error message under the offending widget when validation
//     fails (red micro-label) — does NOT abort the edit (the value remains
//     uncommitted in the bag).
//   * Honours `advanced` flag with a master "Show Advanced" toggle.
// =============================================================================
#pragma once

#include "gui_core/PropertyDescriptor.hpp"

#include <QWidget>
#include <unordered_map>

class QVBoxLayout;
class QFormLayout;
class QCheckBox;
class QGroupBox;

namespace simall::gui
{

class PropertyEditorV2 : public QWidget
{
    Q_OBJECT
public:
    explicit PropertyEditorV2(QWidget* parent = nullptr);

    // Take ownership of the descriptor bag and rebuild the editor.
    void set_bag(std::vector<gui_core::PropertyDescriptor> bag);
    [[nodiscard]] const std::vector<gui_core::PropertyDescriptor>& bag() const { return bag_; }

    void set_show_advanced(bool on);
    void refresh_visibility(); // re-evaluates predicates for all widgets

signals:
    void valueCommitted(const QString& propertyName);
    void validationFailed(const QString& propertyName, const QString& reason);

private:
    QWidget* build_widget_for(gui_core::PropertyDescriptor& d);
    void wire_commit(QWidget* w,
                     gui_core::PropertyDescriptor& d,
                     std::function<gui_core::Variant()> reader);

    QVBoxLayout* root_ = nullptr;
    QCheckBox* showAdvanced_ = nullptr;
    std::vector<gui_core::PropertyDescriptor> bag_;
    std::unordered_map<std::string, QWidget*> widgets_;     // propertyName → widget
    std::unordered_map<std::string, QWidget*> rowParents_;  // propertyName → row container
    std::unordered_map<std::string, QWidget*> errorLabels_; // propertyName → micro-label
    bool advanced_ = false;
};

} // namespace simall::gui
