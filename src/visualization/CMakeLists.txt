# -----------------------------------------------------------------------------
# Visualization & post-processing engine (Week 14 expands this module from
# the original VTK-only viewport into a full domain-pure pipeline:
#
#   Scene management : SceneGraph, ActorRegistry, CameraController
#   Filters          : ContourFilter, IsoSurface, StreamlineRk4, VectorGlyph,
#                      ClippingPlane, SectionCut, VolumeRaycast, LicFilter
#   Overlays         : ScalarBar, Annotation, OrientationGizmo
#   Recording        : ScreenshotRecorder, AnimationRecorder
#   Interaction      : PickingBridge
#
# Every new translation unit is VTK-free.  Viewport remains the only TU that
# depends on QVTKOpenGLNativeWidget / vtkRenderer, so headless CI builds and
# the unit tests link the rest without a display server.
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Visualization & post-processing engine (Week 14).
#
# The module is split into two libraries so unit tests and headless CI nodes
# do not need Qt or VTK at link time:
#
#   SimAll::VisualizationCore — VTK-free.  Scene management, all filters,
#       overlays, recording, picking.  Operates purely on POD geometry views.
#
#   SimAll::Visualization     — VTK + Qt adapter (Viewport).  Links against
#       VisualizationCore so application code can hand SurfaceMesh / LineSet
#       / GlyphSet straight to the renderer.
#
# Contents:
#   Scene management : SceneGraph, ActorRegistry, CameraController
#   Filters          : ContourFilter, IsoSurface, StreamlineRk4, VectorGlyph,
#                      ClippingPlane, SectionCut, VolumeRaycast, LicFilter
#   Overlays         : ScalarBar, Annotation, OrientationGizmo
#   Recording        : ScreenshotRecorder, AnimationRecorder
#   Interaction      : PickingBridge
# -----------------------------------------------------------------------------
simall_add_module(visualization_core
    SOURCES
        SceneGraph.cpp
        ActorRegistry.cpp
        CameraController.cpp
        ContourFilter.cpp
        IsoSurface.cpp
        StreamlineRk4.cpp
        VectorGlyph.cpp
        ClippingPlane.cpp
        SectionCut.cpp
        VolumeRaycast.cpp
        LicFilter.cpp
        ScalarBar.cpp
        Annotation.cpp
        OrientationGizmo.cpp
        ScreenshotRecorder.cpp
        AnimationRecorder.cpp
        PickingBridge.cpp
    HEADERS
        VisualizationTypes.hpp
        SceneGraph.hpp
        ActorRegistry.hpp
        CameraController.hpp
        ContourFilter.hpp
        IsoSurface.hpp
        StreamlineRk4.hpp
        VectorGlyph.hpp
        ClippingPlane.hpp
        SectionCut.hpp
        VolumeRaycast.hpp
        LicFilter.hpp
        ScalarBar.hpp
        Annotation.hpp
        OrientationGizmo.hpp
        ScreenshotRecorder.hpp
        AnimationRecorder.hpp
        PickingBridge.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities)

simall_add_module(visualization
    SOURCES Viewport.cpp
    HEADERS Viewport.hpp
    PUBLIC_DEPS SimAll::VisualizationCore SimAll::Cad SimAll::Meshing
                Qt6::Core Qt6::Widgets Qt6::OpenGLWidgets
                ${VTK_LIBRARIES})

vtk_module_autoinit(TARGETS simall_visualization MODULES ${VTK_LIBRARIES})
