# =============================================================================
# SimAll Beta - Industrial CFD/CAE Platform
# Root CMake Build Script
#
# Spec references:
#   - Master Blueprint: Phase 1.2 (CMAKE INFRASTRUCTURE)
#   - Ultra-Detailed Spec: Section 2 (CORE TECHNOLOGY STACK)
#
# Build modes: Debug | Release | RelWithDebInfo | ASan | UBSan
# Required external packages: Qt6, VTK (>=9.2), OpenCASCADE (>=7.6),
#   Eigen3, optionally MPI, OpenMP, TBB, CUDA, PETSc.
# =============================================================================

cmake_minimum_required(VERSION 3.22)

project(SimAll_Beta
    VERSION      1.0.0
    DESCRIPTION  "SimAll Beta - Industrial Multiphysics CFD/CAE Platform"
    LANGUAGES    CXX C
)

# -----------------------------------------------------------------------------
# Global standards
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD          20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS        OFF)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS   ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(SimAllCompilerFlags)
include(SimAllModule)
include(SimAllDependencies)
include(SimAllInstall)

# Standard CTest module (provides BUILD_TESTING, Dart/CDash hooks).
include(CTest)

# -----------------------------------------------------------------------------
# Project-wide options
# -----------------------------------------------------------------------------
option(SIMALL_BUILD_TESTS         "Build unit + integration tests"        ON)
option(SIMALL_BUILD_APPLICATIONS  "Build executable applications"         ON)
option(SIMALL_BUILD_PLUGINS       "Build bundled plugins"                 ON)
option(SIMALL_ENABLE_MPI          "Enable MPI parallel solver backend"    OFF)
option(SIMALL_ENABLE_OPENMP       "Enable OpenMP shared-memory parallel"  ON)
option(SIMALL_ENABLE_TBB          "Enable Intel TBB task scheduler"       ON)
option(SIMALL_ENABLE_CUDA         "Enable CUDA GPU backend"               OFF)
option(SIMALL_ENABLE_SCRIPTING    "Enable embedded Python / pybind11"     ON)
option(SIMALL_ENABLE_SANITIZERS   "Enable ASan/UBSan in Debug builds"     OFF)
option(SIMALL_USE_AVX2            "Enable AVX2 codegen"                   ON)
option(SIMALL_USE_AVX512          "Enable AVX512 codegen (override AVX2)" OFF)
option(SIMALL_USE_FETCHCONTENT    "Fetch header-only third-party deps via FetchContent" ON)

# -----------------------------------------------------------------------------
# Output layout
# -----------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# -----------------------------------------------------------------------------
# External dependencies (resolved lazily; missing optional deps are skipped)
# -----------------------------------------------------------------------------
if(SIMALL_USE_FETCHCONTENT AND EXISTS "${CMAKE_SOURCE_DIR}/third_party/submodules.cmake")
    include("${CMAKE_SOURCE_DIR}/third_party/submodules.cmake")
endif()

simall_find_required_dependencies()
simall_find_optional_dependencies()

# -----------------------------------------------------------------------------
# Subsystem libraries — STRICT ORDER (Master Blueprint, "IMPLEMENTATION ORDER")
# Each subsystem compiles to its own static/shared library and depends ONLY
# upon already-declared subsystems. Cyclic dependencies are FORBIDDEN.
# -----------------------------------------------------------------------------
add_subdirectory(src/utilities)        # Phase 24: memory, math, SoA primitives
add_subdirectory(src/core)             # Phase 1.3: service locator, events, log, undo
add_subdirectory(src/io)               # Phase 22: serializers (depends: core, utilities)
add_subdirectory(src/parallel)         # Phase 17: MPI / OpenMP / TBB facades
add_subdirectory(src/gpu)              # Phase 18: device abstraction
add_subdirectory(src/materials)        # Phase 15: material database
add_subdirectory(src/cad)              # Phase 4 : OpenCASCADE wrapper
add_subdirectory(src/meshing)          # Phase 5 : mesh data + meshers
add_subdirectory(src/solver)           # Phase 6 : finite-volume infrastructure
add_subdirectory(src/turbulence)       # Phase 8 : RANS/LES/Hybrid plugins
add_subdirectory(src/heat_transfer)    # Phase 9
add_subdirectory(src/radiation)        # Phase 10
add_subdirectory(src/combustion)       # Phase 11
add_subdirectory(src/multiphase)       # Phase 12
add_subdirectory(src/particles)        # Phase 13: Lagrangian DPM
add_subdirectory(src/morphing)         # Phase 14: RBF mesh morphing (FSI)
add_subdirectory(src/ibm)              # Phase 5.7: Immersed Boundary Method
add_subdirectory(src/porous)           # Phase 12.5: Darcy-Forchheimer porous media
add_subdirectory(src/rotating)         # Phase 15: MRF rotating-frame source
add_subdirectory(src/dynamics)         # Phase 14.5: Six-DOF rigid body
add_subdirectory(src/acoustics)        # Phase 18: FW-H aeroacoustics
add_subdirectory(src/emag)             # Phase 19: low-Rm MHD
add_subdirectory(src/rom)              # Phase 21: POD/DMD reduced-order modelling
add_subdirectory(src/zones)            # Phase 22: Fan / heat-exchanger lumped zones
add_subdirectory(src/amr)              # Phase 16: adaptive mesh refinement
add_subdirectory(src/scripting)        # Phase 20
add_subdirectory(src/optimization)     # (future)
add_subdirectory(src/adjoint)          # (future)
add_subdirectory(src/visualization)    # Phase 3 : VTK viewport + actors
add_subdirectory(src/gui_core)         # Phase 2 : Qt-free property + dock + panel registry
add_subdirectory(src/gui)              # Phase 2 : Qt6 shell (depends on everything above)

if(SIMALL_BUILD_PLUGINS)
    add_subdirectory(plugins)
endif()

if(SIMALL_BUILD_APPLICATIONS)
    add_subdirectory(applications)
endif()

if(SIMALL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# -----------------------------------------------------------------------------
# Install + export every subsystem under a single SimAllTargets export set.
# -----------------------------------------------------------------------------
set(SIMALL_SUBSYSTEMS
    utilities core io parallel gpu materials cad meshing solver
    turbulence heat_transfer radiation combustion multiphase particles
    morphing ibm porous rotating dynamics acoustics emag rom zones amr
    scripting optimization adjoint visualization_core visualization gui_core gui)

foreach(_mod IN LISTS SIMALL_SUBSYSTEMS)
    if(TARGET simall_${_mod})
        simall_install_module(${_mod})
    endif()
endforeach()

simall_export_targets()

# -----------------------------------------------------------------------------
# Doxygen target (opt-in: cmake --build build --target docs)
# -----------------------------------------------------------------------------
find_package(Doxygen QUIET)
if(Doxygen_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/docs/Doxyfile")
    add_custom_target(docs
        COMMAND Doxygen::doxygen "${CMAKE_SOURCE_DIR}/docs/Doxyfile"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
    set_target_properties(docs PROPERTIES FOLDER "SimAll/Tools")
endif()

# -----------------------------------------------------------------------------
# CPack — runtime / devel / docs components.
# -----------------------------------------------------------------------------
simall_setup_cpack()

# -----------------------------------------------------------------------------
# Build summary
# -----------------------------------------------------------------------------
message(STATUS "")
message(STATUS "================ SimAll Beta build configuration ================")
message(STATUS "  Build type        : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C++ standard      : ${CMAKE_CXX_STANDARD}")
message(STATUS "  Qt6               : ${Qt6_FOUND}")
message(STATUS "  VTK               : ${VTK_FOUND}  (${VTK_VERSION})")
message(STATUS "  OpenCASCADE       : ${OpenCASCADE_FOUND}")
message(STATUS "  Eigen3            : ${Eigen3_FOUND}")
message(STATUS "  MPI    enabled    : ${SIMALL_ENABLE_MPI}")
message(STATUS "  OpenMP enabled    : ${SIMALL_ENABLE_OPENMP}")
message(STATUS "  TBB    enabled    : ${SIMALL_ENABLE_TBB}")
message(STATUS "  CUDA   enabled    : ${SIMALL_ENABLE_CUDA}")
message(STATUS "  Sanitizers        : ${SIMALL_ENABLE_SANITIZERS}")
message(STATUS "=================================================================")
message(STATUS "")
