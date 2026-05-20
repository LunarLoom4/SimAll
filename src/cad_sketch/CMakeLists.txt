# =============================================================================
# SimAll Beta -- Parametric CAD / Sketcher subsystem
# Phase 23 Pass 23.1
#
# The sketcher is OpenCASCADE-free on purpose: it operates on pure POD
# parameter / entity / constraint data and links only against Eigen3 (for
# the dense Gauss-Newton linear algebra) and SimAll::Utilities.  This lets
# tests run without an OCCT install.
# =============================================================================
simall_add_module(cad_sketch
    SOURCES
        Sketch.cpp
        ConstraintSolver.cpp
    HEADERS
        Vec2.hpp
        Sketch.hpp
        ConstraintSolver.hpp
    PUBLIC_DEPS
        SimAll::Utilities
    PRIVATE_DEPS
        Eigen3::Eigen)
