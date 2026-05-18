// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/Viewport.hpp
// Phase  : 3 (VTK VISUALIZATION ARCHITECTURE)
//
// QVTKOpenGLNativeWidget host. Owns:
//   - vtkRenderer + vtkRenderWindow
//   - camera state (CameraState struct mirrors Section 4.3)
//   - actor registry (CAD edges, CAD surfaces, mesh, scalar postprocessing)
//   - picking adapter producing core::events::SelectionChanged
//
// VTK NEVER owns exact geometry — it owns ONLY vtkPolyData built from a
// cad::TriangleMesh. Picking translates vtkCellId → cad::PersistentId via
// the triangleFaceId table embedded in cad::TriangleMesh.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"
#include "meshing/MeshStorage.hpp"
#include "utilities/MathTypes.hpp"

#include <QVTKOpenGLNativeWidget.h>
#include <vtkSmartPointer.h>

#include <unordered_map>

class vtkRenderer;
class vtkRenderWindow;
class vtkActor;
class vtkPolyData;
class vtkCellPicker;
class vtkScalarBarActor;
class vtkLookupTable;
class vtkOrientationMarkerWidget;

namespace simall::visualization {

struct CameraState {
    util::Vec3d position{0, 0, 10};
    util::Vec3d focalPoint{0, 0, 0};
    util::Vec3d upVector{0, 1, 0};
    double      fieldOfView = 30.0;
    double      nearPlane   = 0.1;
    double      farPlane    = 1000.0;
    bool        orthographic = false;
};

class Viewport : public QVTKOpenGLNativeWidget {
    Q_OBJECT
public:
    explicit Viewport(QWidget* parent = nullptr);
    ~Viewport() override;

    // CAD display
    void show_cad_mesh(const cad::TriangleMesh& tri);
    void clear_cad();

    // Mesh display (Phase 5 visualisation)
    void show_mesh(const meshing::Mesh& mesh);
    void clear_mesh();

    // Scalar visualisation
    void show_scalar_field(const std::string& name,
                           const std::vector<double>& cell_values,
                           double clip_min, double clip_max);

    // Camera
    void          set_camera(const CameraState&);
    CameraState   camera() const;
    void          fit_view();
    void          set_orthographic(bool on);

    // Display toggles
    void          set_edge_overlay(bool on);
    void          set_transparency(double alpha);   // [0,1]

signals:
    void selectionChanged(std::uint64_t topologyId);

protected:
    void mousePressEvent(QMouseEvent* e) override;

private:
    void initialize_renderer();
    void install_axes();
    void pick_at(int x, int y);

    vtkSmartPointer<vtkRenderer>                renderer_;
    vtkSmartPointer<vtkRenderWindow>            window_;
    vtkSmartPointer<vtkActor>                   cad_actor_;
    vtkSmartPointer<vtkActor>                   mesh_actor_;
    vtkSmartPointer<vtkActor>                   scalar_actor_;
    vtkSmartPointer<vtkScalarBarActor>          scalar_bar_;
    vtkSmartPointer<vtkLookupTable>             lut_;
    vtkSmartPointer<vtkCellPicker>              picker_;
    vtkSmartPointer<vtkOrientationMarkerWidget> axes_;

    // VTK triangle index → CAD persistent topology id (Phase 3.3)
    std::vector<util::PersistentId> tri_to_face_;
};

}  // namespace simall::visualization
