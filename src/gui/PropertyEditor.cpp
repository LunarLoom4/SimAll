#include "gui/PropertyEditor.hpp"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>

namespace simall::gui
{

PropertyEditor::PropertyEditor(QWidget* p) : QWidget(p)
{
    form_ = new QFormLayout(this);
    form_->setLabelAlignment(Qt::AlignLeft);
    form_->setFormAlignment(Qt::AlignTop);
    setMinimumWidth(280); // Section 3.5
}

void PropertyEditor::set_properties(std::vector<PropertyDescriptor> ps)
{
    props_ = std::move(ps);
    rebuild();
}

void PropertyEditor::clear()
{
    props_.clear();
    rebuild();
}

void PropertyEditor::rebuild()
{
    while (form_->rowCount())
        form_->removeRow(0);
    for (auto& p : props_) {
        if (p.visibilityCondition && !p.visibilityCondition())
            continue;
        QWidget* w = nullptr;
        switch (p.type) {
        case PropertyType::Double: {
            auto* sb = new QDoubleSpinBox(this);
            sb->setRange(p.minimum, p.maximum);
            sb->setDecimals(6);
            sb->setValue(std::get<double>(p.defaultValue));
            if (!p.units.empty())
                sb->setSuffix(QString(" %1").arg(p.units.c_str()));
            QObject::connect(sb, &QDoubleSpinBox::valueChanged, this, [&](double v) {
                if (p.onChanged)
                    p.onChanged(v);
            });
            w = sb;
            break;
        }
        case PropertyType::Int: {
            auto* sb = new QSpinBox(this);
            sb->setRange(int(p.minimum), int(p.maximum));
            sb->setValue(std::get<int>(p.defaultValue));
            QObject::connect(sb, &QSpinBox::valueChanged, this, [&](int v) {
                if (p.onChanged)
                    p.onChanged(v);
            });
            w = sb;
            break;
        }
        case PropertyType::Bool: {
            auto* cb = new QCheckBox(this);
            cb->setChecked(std::get<bool>(p.defaultValue));
            QObject::connect(cb, &QCheckBox::toggled, this, [&](bool v) {
                if (p.onChanged)
                    p.onChanged(v);
            });
            w = cb;
            break;
        }
        case PropertyType::Enum: {
            auto* combo = new QComboBox(this);
            for (auto& opt : p.enumOptions)
                combo->addItem(QString::fromStdString(opt));
            if (auto* s = std::get_if<std::string>(&p.defaultValue))
                combo->setCurrentText(QString::fromStdString(*s));
            QObject::connect(combo, &QComboBox::currentTextChanged, this, [&](const QString& v) {
                if (p.onChanged)
                    p.onChanged(v.toStdString());
                rebuild(); // re-evaluate visibility
            });
            w = combo;
            break;
        }
        case PropertyType::String: {
            auto* e = new QLineEdit(this);
            if (auto* s = std::get_if<std::string>(&p.defaultValue))
                e->setText(QString::fromStdString(*s));
            QObject::connect(e, &QLineEdit::editingFinished, this, [&, e] {
                if (p.onChanged)
                    p.onChanged(e->text().toStdString());
            });
            w = e;
            break;
        }
        case PropertyType::Vec3:
            // Composed widget would be created here.
            w = new QLabel("(vec3)", this);
            break;
        }
        if (w) {
            w->setToolTip(QString::fromStdString(p.tooltip));
            form_->addRow(QString::fromStdString(p.propertyName), w);
        }
    }
}

} // namespace simall::gui
