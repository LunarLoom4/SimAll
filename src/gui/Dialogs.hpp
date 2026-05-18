// =============================================================================
// SimAll Beta - GUI Subsystem
// File   : src/gui/Dialogs.hpp
//
// Modal dialogs collected here for the W15 GUI build-out.  Each dialog
// owns a `PropertyEditorV2` populated from a domain-specific descriptor
// bag — that gives us validation, units, advanced toggles, and headless
// reuse for free.
//
//   * PreferencesDialog       — UI / performance / units
//   * AboutDialog             — version / build / acknowledgements
//   * ImportCadDialog         — tolerance, healing, scale, deflection
//   * MeshSettingsDialog      — surface mesh + boundary layer params
//   * SolverSetupDialog       — algorithm / schemes / relaxation
//   * MaterialsBrowserDialog  — chooser into the materials registry
//   * BcEditorDialog          — per-BC parameter editor
//   * RunSetupDialog          — iterations / dt / parallel options
// =============================================================================
#pragma once

#include "gui_core/PropertyDescriptor.hpp"

#include <QDialog>
#include <QString>
#include <vector>

class QDialogButtonBox;

namespace simall::gui {

class PropertyEditorV2;

// -----------------------------------------------------------------------------
class PreferencesDialog : public QDialog {
    Q_OBJECT
public:
    explicit PreferencesDialog(QWidget* parent = nullptr);
    [[nodiscard]] const std::vector<gui_core::PropertyDescriptor>& bag() const;
private:
    PropertyEditorV2*  editor_ = nullptr;
    QDialogButtonBox*  buttons_= nullptr;
};

// -----------------------------------------------------------------------------
class AboutDialog : public QDialog {
    Q_OBJECT
public:
    explicit AboutDialog(QWidget* parent = nullptr);
};

// -----------------------------------------------------------------------------
class ImportCadDialog : public QDialog {
    Q_OBJECT
public:
    explicit ImportCadDialog(QWidget* parent = nullptr);
    struct Options {
        QString path;
        double  scale          = 1.0;
        double  sewTolerance   = 1e-4;
        double  deflection     = 1e-3;
        bool    healOnImport   = true;
        bool    defeature      = false;
    };
    [[nodiscard]] Options options() const;
private:
    PropertyEditorV2*  editor_ = nullptr;
};

// -----------------------------------------------------------------------------
class MeshSettingsDialog : public QDialog {
    Q_OBJECT
public:
    explicit MeshSettingsDialog(QWidget* parent = nullptr);
    struct Settings {
        double surfaceSize        = 0.01;
        double minSize            = 1e-4;
        double growthRate         = 1.2;
        int    boundaryLayers     = 8;
        double firstLayerThickness= 1e-5;
        double blGrowthRate       = 1.15;
        bool   conformal          = true;
    };
    [[nodiscard]] Settings settings() const;
private:
    PropertyEditorV2*  editor_ = nullptr;
};

// -----------------------------------------------------------------------------
class SolverSetupDialog : public QDialog {
    Q_OBJECT
public:
    explicit SolverSetupDialog(QWidget* parent = nullptr);
    struct Setup {
        int     algorithmIdx        = 0;  // SIMPLE / SIMPLEC / PISO / Coupled
        int     convectionSchemeIdx = 0;  // Upwind / SOU / QUICK / MUSCL
        double  relaxationP         = 0.3;
        double  relaxationU         = 0.7;
        int     linearSolverIdx     = 0;  // CG / BiCGSTAB / GMRES / AMG
        double  residualTarget      = 1e-5;
        int     maxIterations       = 2000;
    };
    [[nodiscard]] Setup setup() const;
private:
    PropertyEditorV2*  editor_ = nullptr;
};

// -----------------------------------------------------------------------------
class MaterialsBrowserDialog : public QDialog {
    Q_OBJECT
public:
    explicit MaterialsBrowserDialog(const QStringList& available,
                                    QWidget* parent = nullptr);
    [[nodiscard]] QString selected() const { return chosen_; }
private:
    QString chosen_;
};

// -----------------------------------------------------------------------------
class BcEditorDialog : public QDialog {
    Q_OBJECT
public:
    enum class Kind { Wall, Inlet, Outlet, Symmetry, Periodic, PorousJump, Fan, Interface };
    explicit BcEditorDialog(Kind kind, QWidget* parent = nullptr);
    [[nodiscard]] const std::vector<gui_core::PropertyDescriptor>& bag() const;
    [[nodiscard]] Kind kind() const { return kind_; }
private:
    Kind               kind_;
    PropertyEditorV2*  editor_ = nullptr;
};

// -----------------------------------------------------------------------------
class RunSetupDialog : public QDialog {
    Q_OBJECT
public:
    explicit RunSetupDialog(QWidget* parent = nullptr);
    struct Run {
        bool   steady          = true;
        int    iterations      = 1000;
        double endTime         = 1.0;
        double timeStep        = 1e-3;
        int    mpiRanks        = 1;
        int    ompThreads      = 8;
        bool   gpuAcceleration = false;
    };
    [[nodiscard]] Run run() const;
private:
    PropertyEditorV2*  editor_ = nullptr;
};

}  // namespace simall::gui
