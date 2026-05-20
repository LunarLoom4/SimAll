// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/Dialogs.cpp
// =============================================================================
#include "gui/Dialogs.hpp"

#include "gui/PropertyEditorV2.hpp"

#include <QThread>

#include <QApplication>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QSizePolicy>
#include <QVBoxLayout>

namespace simall::gui
{

namespace
{

using gui_core::PropertyDescriptor;
using gui_core::PropertyType;
using gui_core::Variant;

PropertyDescriptor pd_double(std::string name,
                             std::string disp,
                             double def,
                             double mn,
                             double mx,
                             std::string units = {},
                             std::string cat = "General")
{
    PropertyDescriptor d;
    d.propertyName = std::move(name);
    d.displayName = std::move(disp);
    d.type = PropertyType::Double;
    d.defaultValue = def;
    d.currentValue = def;
    d.minimum = mn;
    d.maximum = mx;
    d.units = std::move(units);
    d.category = std::move(cat);
    return d;
}

PropertyDescriptor pd_int(
    std::string name, std::string disp, int def, int mn, int mx, std::string cat = "General")
{
    PropertyDescriptor d;
    d.propertyName = std::move(name);
    d.displayName = std::move(disp);
    d.type = PropertyType::Int;
    d.defaultValue = def;
    d.currentValue = def;
    d.minimum = mn;
    d.maximum = mx;
    d.category = std::move(cat);
    return d;
}

PropertyDescriptor pd_bool(std::string name,
                           std::string disp,
                           bool def,
                           std::string cat = "General")
{
    PropertyDescriptor d;
    d.propertyName = std::move(name);
    d.displayName = std::move(disp);
    d.type = PropertyType::Bool;
    d.defaultValue = def;
    d.currentValue = def;
    d.category = std::move(cat);
    return d;
}

PropertyDescriptor pd_enum(std::string name,
                           std::string disp,
                           std::vector<std::string> opts,
                           int def,
                           std::string cat = "General")
{
    PropertyDescriptor d;
    d.propertyName = std::move(name);
    d.displayName = std::move(disp);
    d.type = PropertyType::Enum;
    d.enumOptions = std::move(opts);
    d.defaultValue = def;
    d.currentValue = def;
    d.minimum = 0;
    d.maximum = double(d.enumOptions.size() - 1);
    d.category = std::move(cat);
    return d;
}

PropertyDescriptor pd_path(std::string name,
                           std::string disp,
                           std::string filter,
                           bool mustExist = false,
                           std::string cat = "General")
{
    PropertyDescriptor d;
    d.propertyName = std::move(name);
    d.displayName = std::move(disp);
    d.type = PropertyType::FilePath;
    d.fileFilter = std::move(filter);
    d.pathMustExist = mustExist;
    d.defaultValue = std::string{};
    d.currentValue = std::string{};
    d.category = std::move(cat);
    return d;
}

QDialog* shell(QWidget* parent,
               const QString& title,
               PropertyEditorV2*& editor,
               std::vector<PropertyDescriptor> bag,
               std::function<void()> onAccept = {})
{
    auto* dlg = new QDialog(parent);
    dlg->setWindowTitle(title);
    dlg->resize(560, 480);
    auto* lay = new QVBoxLayout(dlg);
    editor = new PropertyEditorV2(dlg);
    editor->set_bag(std::move(bag));
    lay->addWidget(editor, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, dlg);
    lay->addWidget(btns);
    QObject::connect(btns, &QDialogButtonBox::accepted, dlg, [dlg, onAccept] {
        if (onAccept)
            onAccept();
        dlg->accept();
    });
    QObject::connect(btns, &QDialogButtonBox::rejected, dlg, &QDialog::reject);
    return dlg;
}

template <class T>
T read_or(const std::vector<PropertyDescriptor>& bag, std::string_view name, T fallback)
{
    auto* d = gui_core::find(bag, name);
    if (!d)
        return fallback;
    if constexpr (std::is_same_v<T, double>)
        return gui_core::variant_as_double(d->currentValue);
    else if constexpr (std::is_same_v<T, int>) {
        if (std::holds_alternative<int>(d->currentValue))
            return std::get<int>(d->currentValue);
        return int(gui_core::variant_as_double(d->currentValue));
    } else if constexpr (std::is_same_v<T, bool>) {
        if (std::holds_alternative<bool>(d->currentValue))
            return std::get<bool>(d->currentValue);
        return fallback;
    } else if constexpr (std::is_same_v<T, QString>) {
        if (std::holds_alternative<std::string>(d->currentValue))
            return QString::fromStdString(std::get<std::string>(d->currentValue));
        return fallback;
    }
    return fallback;
}

} // namespace

// ----------------- PreferencesDialog -------------------------------------
PreferencesDialog::PreferencesDialog(QWidget* parent) : QDialog(parent)
{
    std::vector<PropertyDescriptor> bag;
    bag.push_back(
        pd_enum("ui.theme", "Theme", {"Dark", "Light", "High contrast"}, 0, "Appearance"));
    bag.push_back(pd_int("ui.font", "Editor font size", 10, 7, 24, "Appearance"));
    bag.push_back(pd_bool("ui.autosave", "Autosave project", true, "Workspace"));
    bag.push_back(pd_int("ui.autoFreq", "Autosave every (s)", 60, 5, 3600, "Workspace"));
    bag.push_back(pd_enum("units.length", "Length units", {"m", "mm", "in"}, 0, "Units"));
    bag.push_back(
        pd_enum("units.pressure", "Pressure units", {"Pa", "bar", "atm", "psi"}, 0, "Units"));
    bag.push_back(pd_enum("units.temperature", "Temperature units", {"K", "°C", "°F"}, 0, "Units"));
    bag.push_back(pd_int(
        "perf.threads", "Worker threads", QThread::idealThreadCount(), 1, 256, "Performance"));
    bag.push_back(pd_bool("perf.gpu", "Enable GPU acceleration", false, "Performance"));
    setWindowTitle("Preferences");
    resize(560, 520);
    auto* lay = new QVBoxLayout(this);
    editor_ = new PropertyEditorV2(this);
    editor_->set_bag(std::move(bag));
    lay->addWidget(editor_, 1);
    buttons_ = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(buttons_);
    connect(buttons_, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttons_, &QDialogButtonBox::rejected, this, &QDialog::reject);
}
const std::vector<PropertyDescriptor>& PreferencesDialog::bag() const
{
    return editor_->bag();
}

// ----------------- AboutDialog -------------------------------------------
AboutDialog::AboutDialog(QWidget* parent) : QDialog(parent)
{
    setWindowTitle("About SimAll Beta");
    setFixedSize(480, 320);
    auto* lay = new QVBoxLayout(this);
    auto* title = new QLabel("<h2>SimAll Beta</h2>", this);
    title->setAlignment(Qt::AlignCenter);
    auto* body = new QLabel(
        "<p style='text-align:center;'>Industrial multiphysics platform.<br>"
        "Build " __DATE__ " " __TIME__ "</p>"
        "<p style='text-align:center;'>"
        "Qt 6 · VTK · OpenCASCADE · Eigen · spdlog · METIS · Scotch</p>"
        "<p style='text-align:center;color:#888'>© SimAll project — for research use.</p>",
        this);
    body->setWordWrap(true);
    lay->addWidget(title);
    lay->addWidget(body, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Close, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
}

// ----------------- ImportCadDialog ---------------------------------------
ImportCadDialog::ImportCadDialog(QWidget* parent) : QDialog(parent)
{
    std::vector<PropertyDescriptor> bag;
    bag.push_back(pd_path("cad.path",
                          "CAD file",
                          "CAD files (*.step *.stp *.iges *.igs *.brep *.stl)",
                          true,
                          "Source"));
    bag.push_back(pd_double("cad.scale", "Unit scale", 1.0, 1e-6, 1e6, "", "Geometry"));
    bag.push_back(pd_double("cad.sewTol", "Sewing tolerance", 1e-4, 1e-9, 1e-1, "m", "Healing"));
    bag.push_back(
        pd_double("cad.deflection", "Tessellation deflection", 1e-3, 1e-9, 1e-1, "m", "Healing"));
    bag.push_back(pd_bool("cad.heal", "Heal on import", true, "Healing"));
    bag.push_back(pd_bool("cad.defeature", "Defeature small features", false, "Healing"));
    setWindowTitle("Import CAD");
    resize(560, 420);
    auto* lay = new QVBoxLayout(this);
    editor_ = new PropertyEditorV2(this);
    editor_->set_bag(std::move(bag));
    lay->addWidget(editor_, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}
ImportCadDialog::Options ImportCadDialog::options() const
{
    Options o;
    const auto& b = editor_->bag();
    o.path = read_or<QString>(b, "cad.path", {});
    o.scale = read_or<double>(b, "cad.scale", 1.0);
    o.sewTolerance = read_or<double>(b, "cad.sewTol", 1e-4);
    o.deflection = read_or<double>(b, "cad.deflection", 1e-3);
    o.healOnImport = read_or<bool>(b, "cad.heal", true);
    o.defeature = read_or<bool>(b, "cad.defeature", false);
    return o;
}

// ----------------- MeshSettingsDialog ------------------------------------
MeshSettingsDialog::MeshSettingsDialog(QWidget* parent) : QDialog(parent)
{
    std::vector<PropertyDescriptor> bag;
    bag.push_back(
        pd_double("mesh.surface", "Target surface size", 0.01, 1e-6, 1e3, "m", "Surface"));
    bag.push_back(pd_double("mesh.min", "Minimum size", 1e-4, 1e-9, 1e3, "m", "Surface"));
    bag.push_back(pd_double("mesh.growth", "Growth rate", 1.2, 1.0, 3.0, "", "Surface"));
    bag.push_back(pd_int("mesh.blLayers", "Boundary layers", 8, 0, 100, "Boundary layer"));
    bag.push_back(pd_double(
        "mesh.blFirst", "First-layer thickness", 1e-5, 1e-10, 1e-1, "m", "Boundary layer"));
    bag.push_back(
        pd_double("mesh.blGrowth", "BL growth rate", 1.15, 1.0, 3.0, "", "Boundary layer"));
    bag.push_back(pd_bool("mesh.conformal", "Conformal interfaces", true, "Topology"));
    setWindowTitle("Mesh settings");
    resize(560, 440);
    auto* lay = new QVBoxLayout(this);
    editor_ = new PropertyEditorV2(this);
    editor_->set_bag(std::move(bag));
    lay->addWidget(editor_, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}
MeshSettingsDialog::Settings MeshSettingsDialog::settings() const
{
    Settings s;
    const auto& b = editor_->bag();
    s.surfaceSize = read_or<double>(b, "mesh.surface", 0.01);
    s.minSize = read_or<double>(b, "mesh.min", 1e-4);
    s.growthRate = read_or<double>(b, "mesh.growth", 1.2);
    s.boundaryLayers = read_or<int>(b, "mesh.blLayers", 8);
    s.firstLayerThickness = read_or<double>(b, "mesh.blFirst", 1e-5);
    s.blGrowthRate = read_or<double>(b, "mesh.blGrowth", 1.15);
    s.conformal = read_or<bool>(b, "mesh.conformal", true);
    return s;
}

// ----------------- SolverSetupDialog -------------------------------------
SolverSetupDialog::SolverSetupDialog(QWidget* parent) : QDialog(parent)
{
    std::vector<PropertyDescriptor> bag;
    bag.push_back(pd_enum("solver.algo",
                          "Algorithm",
                          {"SIMPLE", "SIMPLEC", "PISO", "Coupled"},
                          0,
                          "Pressure-velocity"));
    bag.push_back(pd_enum("solver.conv",
                          "Convection scheme",
                          {"Upwind", "SOU", "QUICK", "MUSCL", "CDS"},
                          1,
                          "Discretisation"));
    bag.push_back(pd_double("solver.relP", "Pressure relaxation", 0.3, 0.0, 1.0, "", "Relaxation"));
    bag.push_back(pd_double("solver.relU", "Velocity relaxation", 0.7, 0.0, 1.0, "", "Relaxation"));
    bag.push_back(pd_enum(
        "solver.lin", "Linear solver", {"CG", "BiCGSTAB", "GMRES", "TFQMR", "AMG"}, 1, "Linear"));
    bag.push_back(pd_double("solver.tol", "Residual target", 1e-5, 1e-15, 1e-1, "", "Convergence"));
    bag.push_back(pd_int("solver.iter", "Max iterations", 2000, 1, 1000000, "Convergence"));
    setWindowTitle("Solver setup");
    resize(560, 460);
    auto* lay = new QVBoxLayout(this);
    editor_ = new PropertyEditorV2(this);
    editor_->set_bag(std::move(bag));
    lay->addWidget(editor_, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}
SolverSetupDialog::Setup SolverSetupDialog::setup() const
{
    Setup s;
    const auto& b = editor_->bag();
    s.algorithmIdx = read_or<int>(b, "solver.algo", 0);
    s.convectionSchemeIdx = read_or<int>(b, "solver.conv", 1);
    s.relaxationP = read_or<double>(b, "solver.relP", 0.3);
    s.relaxationU = read_or<double>(b, "solver.relU", 0.7);
    s.linearSolverIdx = read_or<int>(b, "solver.lin", 1);
    s.residualTarget = read_or<double>(b, "solver.tol", 1e-5);
    s.maxIterations = read_or<int>(b, "solver.iter", 2000);
    return s;
}

// ----------------- MaterialsBrowserDialog --------------------------------
MaterialsBrowserDialog::MaterialsBrowserDialog(const QStringList& available, QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle("Materials browser");
    resize(420, 480);
    auto* lay = new QVBoxLayout(this);
    auto* list = new QListWidget(this);
    list->addItems(available);
    lay->addWidget(list, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::accepted, this, [this, list] {
        if (auto* it = list->currentItem()) {
            chosen_ = it->text();
            accept();
        } else
            reject();
    });
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
    connect(list, &QListWidget::itemDoubleClicked, this, [this](QListWidgetItem* it) {
        chosen_ = it->text();
        accept();
    });
}

// ----------------- BcEditorDialog ----------------------------------------
namespace
{
std::vector<PropertyDescriptor> make_bc_bag(BcEditorDialog::Kind k)
{
    using K = BcEditorDialog::Kind;
    std::vector<PropertyDescriptor> b;
    switch (k) {
    case K::Wall:
        b.push_back(pd_enum("wall.type", "Type", {"No-slip", "Slip", "Moving"}, 0, "Mechanics"));
        b.push_back(pd_double("wall.rough", "Roughness height", 0.0, 0.0, 1.0, "m", "Mechanics"));
        b.push_back(pd_enum("wall.heat",
                            "Thermal",
                            {"Adiabatic", "Fixed temperature", "Heat flux", "HTC"},
                            0,
                            "Heat"));
        b.push_back(pd_double("wall.value", "Thermal value", 300.0, -1e6, 1e6, "", "Heat"));
        break;
    case K::Inlet:
        b.push_back(pd_enum("inlet.type", "Type", {"Velocity", "Mass flow", "Total pressure"}, 0));
        b.push_back(pd_double("inlet.value", "Magnitude", 1.0, -1e6, 1e6));
        b.push_back(pd_double("inlet.T", "Temperature", 300.0, 0, 1e5, "K"));
        b.push_back(pd_double("inlet.k", "Turbulence k", 1e-3, 0, 1e6, "m²/s²"));
        b.push_back(pd_double("inlet.eps", "Turbulence ε", 1e-3, 0, 1e9, "m²/s³"));
        break;
    case K::Outlet:
        b.push_back(pd_enum("outlet.type", "Type", {"Pressure", "Outflow", "Mass flow"}, 0));
        b.push_back(pd_double("outlet.p", "Pressure", 0.0, -1e9, 1e9, "Pa"));
        b.push_back(pd_bool("outlet.backflow", "Backflow corrections", true));
        break;
    case K::Symmetry:
        b.push_back(pd_bool("sym.enable", "Enabled", true));
        break;
    case K::Periodic:
        b.push_back(pd_enum("per.kind", "Kind", {"Translational", "Rotational"}, 0));
        b.push_back(pd_double("per.angle", "Rotation angle (deg)", 0.0, -360.0, 360.0));
        break;
    case K::PorousJump:
        b.push_back(pd_double("pj.permeability", "Permeability", 1e-10, 0.0, 1.0, "m²"));
        b.push_back(pd_double("pj.thickness", "Thickness", 0.01, 0.0, 10.0, "m"));
        b.push_back(pd_double("pj.coefC2", "Inertial C2", 0.0, 0.0, 1e9, "1/m"));
        break;
    case K::Fan:
        b.push_back(pd_double("fan.dp", "Δp (Pa)", 0.0, -1e9, 1e9));
        b.push_back(pd_double("fan.rpm", "Rotational speed", 0.0, -1e6, 1e6, "rpm"));
        break;
    case K::Interface:
        b.push_back(pd_enum("if.kind", "Coupling", {"Conformal", "Mortar", "Overset"}, 0));
        b.push_back(pd_bool("if.cons", "Conservative interpolation", true));
        break;
    }
    return b;
}
} // namespace

BcEditorDialog::BcEditorDialog(Kind k, QWidget* parent) : QDialog(parent), kind_(k)
{
    const char* titles[] = {
        "Wall", "Inlet", "Outlet", "Symmetry", "Periodic", "Porous Jump", "Fan", "Interface"};
    setWindowTitle(QString("Boundary condition — %1").arg(titles[int(k)]));
    resize(540, 420);
    auto* lay = new QVBoxLayout(this);
    editor_ = new PropertyEditorV2(this);
    editor_->set_bag(make_bc_bag(k));
    lay->addWidget(editor_, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}
const std::vector<PropertyDescriptor>& BcEditorDialog::bag() const
{
    return editor_->bag();
}

// ----------------- RunSetupDialog ----------------------------------------
RunSetupDialog::RunSetupDialog(QWidget* parent) : QDialog(parent)
{
    std::vector<PropertyDescriptor> bag;
    bag.push_back(pd_bool("run.steady", "Steady state", true, "Time"));
    bag.push_back(pd_int("run.iter", "Iterations", 1000, 1, 1000000000, "Time"));
    bag.push_back(pd_double("run.endTime", "End time", 1.0, 0.0, 1e15, "s", "Time"));
    bag.push_back(pd_double("run.dt", "Time step", 1e-3, 1e-15, 1e6, "s", "Time"));
    bag.push_back(pd_int("run.mpi", "MPI ranks", 1, 1, 65536, "Parallel"));
    bag.push_back(pd_int("run.omp", "OpenMP threads", 8, 1, 1024, "Parallel"));
    bag.push_back(pd_bool("run.gpu", "Enable GPU", false, "Parallel"));
    // visibility: time step / end time only when transient.
    auto* steadyDesc = &bag[0];
    for (auto& d : bag) {
        if (d.propertyName == "run.endTime" || d.propertyName == "run.dt") {
            d.visibilityCondition = [steadyDesc] {
                return std::holds_alternative<bool>(steadyDesc->currentValue)
                       && !std::get<bool>(steadyDesc->currentValue);
            };
        }
        if (d.propertyName == "run.iter") {
            d.visibilityCondition = [steadyDesc] {
                return std::holds_alternative<bool>(steadyDesc->currentValue)
                       && std::get<bool>(steadyDesc->currentValue);
            };
        }
    }
    setWindowTitle("Run setup");
    resize(540, 440);
    auto* lay = new QVBoxLayout(this);
    editor_ = new PropertyEditorV2(this);
    editor_->set_bag(std::move(bag));
    lay->addWidget(editor_, 1);
    auto* btns = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    lay->addWidget(btns);
    connect(btns, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(btns, &QDialogButtonBox::rejected, this, &QDialog::reject);
}
RunSetupDialog::Run RunSetupDialog::run() const
{
    Run r;
    const auto& b = editor_->bag();
    r.steady = read_or<bool>(b, "run.steady", true);
    r.iterations = read_or<int>(b, "run.iter", 1000);
    r.endTime = read_or<double>(b, "run.endTime", 1.0);
    r.timeStep = read_or<double>(b, "run.dt", 1e-3);
    r.mpiRanks = read_or<int>(b, "run.mpi", 1);
    r.ompThreads = read_or<int>(b, "run.omp", 8);
    r.gpuAcceleration = read_or<bool>(b, "run.gpu", false);
    return r;
}

} // namespace simall::gui
