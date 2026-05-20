// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/PropertyEditorV2.cpp
// =============================================================================
#include "gui/PropertyEditorV2.hpp"

#include <QToolButton>

#include <QCheckBox>
#include <QColorDialog>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QScrollArea>
#include <QSpinBox>
#include <QVBoxLayout>

namespace simall::gui
{

namespace
{

QString to_qstring(const gui_core::Variant& v)
{
    QString out;
    std::visit(
        [&](const auto& x) {
            using T = std::decay_t<decltype(x)>;
            if constexpr (std::is_same_v<T, double>)
                out = QString::number(x, 'g', 8);
            else if constexpr (std::is_same_v<T, int>)
                out = QString::number(x);
            else if constexpr (std::is_same_v<T, bool>)
                out = x ? "true" : "false";
            else if constexpr (std::is_same_v<T, std::string>)
                out = QString::fromStdString(x);
            else if constexpr (std::is_same_v<T, std::array<double, 3>>)
                out = QString("%1, %2, %3").arg(x[0]).arg(x[1]).arg(x[2]);
        },
        v);
    return out;
}

} // namespace

PropertyEditorV2::PropertyEditorV2(QWidget* parent) : QWidget(parent)
{
    auto* outer = new QVBoxLayout(this);
    outer->setContentsMargins(2, 2, 2, 2);
    outer->setSpacing(4);

    auto* topRow = new QHBoxLayout();
    showAdvanced_ = new QCheckBox("Show advanced", this);
    topRow->addWidget(showAdvanced_);
    topRow->addStretch(1);
    outer->addLayout(topRow);

    auto* scroll = new QScrollArea(this);
    scroll->setWidgetResizable(true);
    auto* inner = new QWidget(scroll);
    root_ = new QVBoxLayout(inner);
    root_->setContentsMargins(4, 4, 4, 4);
    root_->setSpacing(6);
    scroll->setWidget(inner);
    outer->addWidget(scroll, 1);

    connect(showAdvanced_, &QCheckBox::toggled, this, &PropertyEditorV2::set_show_advanced);
}

void PropertyEditorV2::set_show_advanced(bool on)
{
    advanced_ = on;
    refresh_visibility();
}

void PropertyEditorV2::set_bag(std::vector<gui_core::PropertyDescriptor> bag)
{
    bag_ = std::move(bag);
    widgets_.clear();
    rowParents_.clear();
    errorLabels_.clear();

    // Clear current rows.
    QLayoutItem* item;
    while ((item = root_->takeAt(0)) != nullptr) {
        if (auto* w = item->widget())
            w->deleteLater();
        delete item;
    }

    // Group by category (preserve insertion order).
    std::vector<QString> order;
    std::unordered_map<QString, QFormLayout*> formByCat;
    for (auto& d : bag_) {
        const QString cat = QString::fromStdString(d.category.empty() ? "General" : d.category);
        if (formByCat.find(cat) == formByCat.end()) {
            auto* gb = new QGroupBox(cat, this);
            auto* form = new QFormLayout(gb);
            form->setLabelAlignment(Qt::AlignRight);
            form->setFormAlignment(Qt::AlignLeft | Qt::AlignTop);
            root_->addWidget(gb);
            formByCat[cat] = form;
            order.push_back(cat);
        }
        QWidget* w = build_widget_for(d);
        if (!w)
            continue;

        auto* rowHolder = new QWidget(this);
        auto* rowL = new QVBoxLayout(rowHolder);
        rowL->setContentsMargins(0, 0, 0, 0);
        rowL->setSpacing(1);
        rowL->addWidget(w);
        auto* err = new QLabel(rowHolder);
        err->setStyleSheet("color:#E52B50;font-size:9pt;");
        err->hide();
        rowL->addWidget(err);

        const QString label =
            QString::fromStdString(d.display())
            + (d.units.empty() ? "" : "  [" + QString::fromStdString(d.units) + "]");
        formByCat[cat]->addRow(label, rowHolder);

        widgets_[d.propertyName] = w;
        rowParents_[d.propertyName] = rowHolder;
        errorLabels_[d.propertyName] = err;
        if (!d.tooltip.empty())
            w->setToolTip(QString::fromStdString(d.tooltip));
    }
    root_->addStretch(1);
    refresh_visibility();
}

void PropertyEditorV2::refresh_visibility()
{
    for (auto& d : bag_) {
        auto it = rowParents_.find(d.propertyName);
        if (it == rowParents_.end())
            continue;
        const bool vis = d.is_visible() && (!d.advanced || advanced_);
        it->second->setVisible(vis);
        if (auto wit = widgets_.find(d.propertyName); wit != widgets_.end())
            wit->second->setEnabled(d.is_enabled());
    }
}

void PropertyEditorV2::wire_commit(QWidget* w,
                                   gui_core::PropertyDescriptor& d,
                                   std::function<gui_core::Variant()> reader)
{
    const std::string name = d.propertyName;
    auto trigger = [this, name, reader] {
        auto* descPtr = gui_core::find(bag_, name);
        if (!descPtr)
            return;
        auto candidate = reader();
        auto v = descPtr->validate(candidate);
        QWidget* errLbl = errorLabels_[name];
        if (!v.ok) {
            if (auto* lbl = qobject_cast<QLabel*>(errLbl)) {
                lbl->setText(QString::fromStdString(v.message));
                lbl->show();
            }
            emit validationFailed(QString::fromStdString(name), QString::fromStdString(v.message));
            return;
        }
        descPtr->commit(std::move(candidate));
        if (auto* lbl = qobject_cast<QLabel*>(errLbl))
            lbl->hide();
        emit valueCommitted(QString::fromStdString(name));
        refresh_visibility();
    };

    // Hook up the appropriate signal based on the QWidget's runtime type.
    if (auto* sb = qobject_cast<QDoubleSpinBox*>(w))
        QObject::connect(sb, QOverload<double>::of(&QDoubleSpinBox::valueChanged), w, trigger);
    else if (auto* sb = qobject_cast<QSpinBox*>(w))
        QObject::connect(sb, QOverload<int>::of(&QSpinBox::valueChanged), w, trigger);
    else if (auto* cb = qobject_cast<QCheckBox*>(w))
        QObject::connect(cb, &QCheckBox::toggled, w, trigger);
    else if (auto* cmb = qobject_cast<QComboBox*>(w))
        QObject::connect(cmb, QOverload<int>::of(&QComboBox::currentIndexChanged), w, trigger);
    else if (auto* le = qobject_cast<QLineEdit*>(w))
        QObject::connect(le, &QLineEdit::editingFinished, w, trigger);
}

QWidget* PropertyEditorV2::build_widget_for(gui_core::PropertyDescriptor& d)
{
    using gui_core::PropertyType;
    using gui_core::Variant;

    switch (d.type) {
    case PropertyType::Double: {
        auto* sb = new QDoubleSpinBox(this);
        sb->setRange(d.minimum, d.maximum);
        sb->setDecimals(6);
        sb->setSingleStep(d.step > 0 ? d.step : 0.1);
        sb->setValue(gui_core::variant_as_double(d.currentValue));
        wire_commit(sb, d, [sb] { return Variant{double{sb->value()}}; });
        return sb;
    }
    case PropertyType::Int: {
        auto* sb = new QSpinBox(this);
        sb->setRange(int(std::max<double>(d.minimum, INT_MIN)),
                     int(std::min<double>(d.maximum, INT_MAX)));
        sb->setSingleStep(int(d.step > 0 ? d.step : 1));
        sb->setValue(int(gui_core::variant_as_double(d.currentValue)));
        wire_commit(sb, d, [sb] { return Variant{int{sb->value()}}; });
        return sb;
    }
    case PropertyType::Bool: {
        auto* cb = new QCheckBox(this);
        cb->setChecked(std::holds_alternative<bool>(d.currentValue) ? std::get<bool>(d.currentValue)
                                                                    : false);
        wire_commit(cb, d, [cb] { return Variant{bool{cb->isChecked()}}; });
        return cb;
    }
    case PropertyType::Enum: {
        auto* cmb = new QComboBox(this);
        for (const auto& s : d.enumOptions)
            cmb->addItem(QString::fromStdString(s));
        int cur = std::holds_alternative<int>(d.currentValue) ? std::get<int>(d.currentValue) : 0;
        cmb->setCurrentIndex(std::clamp(cur, 0, cmb->count() - 1));
        wire_commit(cmb, d, [cmb] { return Variant{int{cmb->currentIndex()}}; });
        return cmb;
    }
    case PropertyType::String: {
        auto* le = new QLineEdit(this);
        if (std::holds_alternative<std::string>(d.currentValue))
            le->setText(QString::fromStdString(std::get<std::string>(d.currentValue)));
        wire_commit(le, d, [le] { return Variant{le->text().toStdString()}; });
        return le;
    }
    case PropertyType::Vec3:
    case PropertyType::Range: {
        auto* host = new QWidget(this);
        auto* row = new QHBoxLayout(host);
        row->setContentsMargins(0, 0, 0, 0);
        auto* x = new QDoubleSpinBox(host);
        auto* y = new QDoubleSpinBox(host);
        auto* z = new QDoubleSpinBox(host);
        for (auto* sb : {x, y, z}) {
            sb->setRange(d.minimum, d.maximum);
            sb->setDecimals(6);
            sb->setSingleStep(d.step > 0 ? d.step : 0.1);
            row->addWidget(sb);
        }
        std::array<double, 3> v = std::holds_alternative<std::array<double, 3>>(d.currentValue)
                                      ? std::get<std::array<double, 3>>(d.currentValue)
                                      : std::array<double, 3>{0, 0, 0};
        x->setValue(v[0]);
        y->setValue(v[1]);
        z->setValue(v[2]);
        auto reader = [x, y, z] {
            return Variant{std::array<double, 3>{x->value(), y->value(), z->value()}};
        };
        for (auto* sb : {x, y, z})
            wire_commit(sb, d, reader);
        return host;
    }
    case PropertyType::Color: {
        auto* host = new QWidget(this);
        auto* row = new QHBoxLayout(host);
        row->setContentsMargins(0, 0, 0, 0);
        auto* swatch = new QPushButton(host);
        swatch->setFixedSize(48, 22);
        auto* le = new QLineEdit(host);
        row->addWidget(swatch);
        row->addWidget(le, 1);
        std::array<double, 3> v = std::holds_alternative<std::array<double, 3>>(d.currentValue)
                                      ? std::get<std::array<double, 3>>(d.currentValue)
                                      : std::array<double, 3>{1, 1, 1};
        auto applyColour = [swatch, le](std::array<double, 3> c) {
            QColor q(int(c[0] * 255), int(c[1] * 255), int(c[2] * 255));
            swatch->setStyleSheet(QString("background:%1;border:1px solid #555;").arg(q.name()));
            le->setText(q.name());
        };
        applyColour(v);
        QObject::connect(swatch,
                         &QPushButton::clicked,
                         host,
                         [this, host, le, swatch, applyColour, name = d.propertyName] {
                             auto* dsc = gui_core::find(bag_, name);
                             if (!dsc)
                                 return;
                             QColor cur(le->text());
                             QColor c =
                                 QColorDialog::getColor(cur.isValid() ? cur : Qt::white, host);
                             if (!c.isValid())
                                 return;
                             std::array<double, 3> arr{c.redF(), c.greenF(), c.blueF()};
                             applyColour(arr);
                             auto r = dsc->commit(gui_core::Variant{arr});
                             if (r.ok)
                                 emit valueCommitted(QString::fromStdString(name));
                         });
        return host;
    }
    case PropertyType::FilePath: {
        auto* host = new QWidget(this);
        auto* row = new QHBoxLayout(host);
        row->setContentsMargins(0, 0, 0, 0);
        auto* le = new QLineEdit(host);
        auto* btn = new QPushButton("…", host);
        btn->setFixedWidth(28);
        row->addWidget(le, 1);
        row->addWidget(btn);
        if (std::holds_alternative<std::string>(d.currentValue))
            le->setText(QString::fromStdString(std::get<std::string>(d.currentValue)));
        wire_commit(le, d, [le] { return Variant{le->text().toStdString()}; });
        QObject::connect(btn, &QPushButton::clicked, host, [host, le, filter = d.fileFilter] {
            QString f = QFileDialog::getOpenFileName(
                host, "Select file", {}, QString::fromStdString(filter));
            if (!f.isEmpty()) {
                le->setText(f);
                emit le->editingFinished();
            }
        });
        return host;
    }
    }
    return nullptr;
}

} // namespace simall::gui
