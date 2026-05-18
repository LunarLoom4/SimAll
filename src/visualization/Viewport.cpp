// =============================================================================
// SimAll Beta - Visualization Subsystem
// File   : src/visualization/Viewport.cpp
// =============================================================================
#include "visualization/Viewport.hpp"
#include "core/EventBus.hpp"

#include <QMouseEvent>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkCellPicker.h>
#include <vtkCamera.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkGenericOpenGLRenderWindow.h>

namespace simall::visualization {

Viewport::Viewport(QWidget* parent) : QVTKOpenGLNativeWidget(parent) {
    initialize_renderer();
}

Viewport::~Viewport() = default;

void Viewport::initialize_renderer() {
    window_   = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
    setRenderWindow(window_);
    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderer_->SetBackground(0.117, 0.117, 0.117);     // #1E1E1E (Section 3.6)
    renderer_->SetBackground2(0.145, 0.145, 0.149);
    renderer_->GradientBackgroundOn();
    window_->AddRenderer(renderer_);

    auto style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    window_->GetInteractor()->SetInteractorStyle(style);

    picker_ = vtkSmartPointer<vtkCellPicker>::New();
    picker_->SetTolerance(0.0005);

    install_axes();
}

void Viewport::install_axes() {
    auto axes = vtkSmartPointer<vtkAxesActor>::New();
    axes_     = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    axes_->SetOrientationMarker(axes);
    axes_->SetInteractor(window_->GetInteractor());
    axes_->SetViewport(0.0, 0.0, 0.18, 0.22);
    axes_->SetEnabled(true);
    axes_->InteractiveOff();
}

void Viewport::show_cad_mesh(const cad::TriangleMesh& tri) {
    clear_cad();

    auto pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetNumberOfPoints(static_cast<vtkIdType>(tri.points.size()));
    for (vtkIdType i = 0; i < pts->GetNumberOfPoints(); ++i)
        pts->SetPoint(i, tri.points[i].x, tri.points[i].y, tri.points[i].z);

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->AllocateEstimate(tri.triangles.size(), 3);
    for (const auto& t : tri.triangles) {
        vtkIdType ids[3] = {t[0], t[1], t[2]};
        cells->InsertNextCell(3, ids);
    }

    auto poly = vtkSmartPointer<vtkPolyData>::New();
    poly->SetPoints(pts);
    poly->SetPolys(cells);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(poly);

    cad_actor_ = vtkSmartPointer<vtkActor>::New();
    cad_actor_->SetMapper(mapper);
    cad_actor_->GetProperty()->SetColor(0.75, 0.78, 0.82);
    cad_actor_->GetProperty()->SetAmbient(0.25);
    cad_actor_->GetProperty()->SetDiffuse(0.75);
    cad_actor_->GetProperty()->SetSpecular(0.3);
    cad_actor_->GetProperty()->EdgeVisibilityOn();
    cad_actor_->GetProperty()->SetEdgeColor(0.3, 0.55, 0.85);

    renderer_->AddActor(cad_actor_);
    tri_to_face_ = tri.triangleFaceId;

    fit_view();
    window_->Render();
}

void Viewport::clear_cad() {
    if (cad_actor_) { renderer_->RemoveActor(cad_actor_); cad_actor_ = nullptr; }
    tri_to_face_.clear();
}

void Viewport::show_mesh(const meshing::Mesh&) { /* Phase 5 viz */ }
void Viewport::clear_mesh() { if (mesh_actor_) { renderer_->RemoveActor(mesh_actor_); mesh_actor_ = nullptr; } }

void Viewport::show_scalar_field(const std::string&,
                                 const std::vector<double>&,
                                 double, double) {
    // Implemented in Phase 19 with proper LUT mapping + scalar-bar legend.
}

void Viewport::set_camera(const CameraState& c) {
    auto cam = renderer_->GetActiveCamera();
    cam->SetPosition  (c.position.x,   c.position.y,   c.position.z);
    cam->SetFocalPoint(c.focalPoint.x, c.focalPoint.y, c.focalPoint.z);
    cam->SetViewUp    (c.upVector.x,   c.upVector.y,   c.upVector.z);
    cam->SetViewAngle (c.fieldOfView);
    cam->SetClippingRange(c.nearPlane, c.farPlane);
    cam->SetParallelProjection(c.orthographic ? 1 : 0);
    window_->Render();
}

CameraState Viewport::camera() const {
    CameraState c;
    auto cam = renderer_->GetActiveCamera();
    double p[3], f[3], u[3];
    cam->GetPosition(p); cam->GetFocalPoint(f); cam->GetViewUp(u);
    c.position   = {p[0],p[1],p[2]};
    c.focalPoint = {f[0],f[1],f[2]};
    c.upVector   = {u[0],u[1],u[2]};
    c.fieldOfView = cam->GetViewAngle();
    c.orthographic = cam->GetParallelProjection() != 0;
    return c;
}

void Viewport::fit_view()                    { renderer_->ResetCamera(); window_->Render(); }
void Viewport::set_orthographic(bool on)     { renderer_->GetActiveCamera()->SetParallelProjection(on); window_->Render(); }
void Viewport::set_edge_overlay(bool on)     { if (cad_actor_) cad_actor_->GetProperty()->SetEdgeVisibility(on); window_->Render(); }
void Viewport::set_transparency(double alpha){ if (cad_actor_) cad_actor_->GetProperty()->SetOpacity(alpha); window_->Render(); }

void Viewport::mousePressEvent(QMouseEvent* e) {
    QVTKOpenGLNativeWidget::mousePressEvent(e);
    if (e->button() == Qt::LeftButton && (e->modifiers() & Qt::ControlModifier))
        pick_at(e->pos().x(), height() - e->pos().y());
}

void Viewport::pick_at(int x, int y) {
    if (!picker_->Pick(x, y, 0, renderer_)) return;
    const vtkIdType cellId = picker_->GetCellId();
    if (cellId < 0 || cellId >= static_cast<vtkIdType>(tri_to_face_.size())) return;
    const auto topoId = tri_to_face_[cellId];
    emit selectionChanged(topoId);
    core::EventBus::instance().publish(
        core::events::SelectionChanged{ std::vector<std::uint64_t>{topoId} });
}

}  // namespace simall::visualization
