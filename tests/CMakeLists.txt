# =============================================================================
# Tests — Phase 23. Unit tests use Catch2 (FetchContent), integration tests
# drop into tests/integration/, verification cases into tests/verification/
# (lid cavity, backward-facing step, Taylor-Green vortex, etc. per spec §23).
# =============================================================================
include(FetchContent)
FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.4)
FetchContent_MakeAvailable(Catch2)

add_executable(simall_core_tests
    unit/test_event_bus.cpp
    unit/test_undo.cpp
    unit/test_aligned.cpp
    unit/core/test_thread_pool.cpp
    unit/core/test_task_graph.cpp
    unit/core/test_config_store.cpp
    unit/core/test_unit_system.cpp
    unit/core/test_composite_command.cpp
    unit/core/test_rotating_sink.cpp
    unit/cad/test_topology_graph.cpp
    ${CMAKE_SOURCE_DIR}/src/cad/TopologyGraph.cpp)
target_link_libraries(simall_core_tests PRIVATE
    SimAll::CompilerFlags SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Materials unit tests (Week 4). Links SimAll::Materials directly so we get
# Janaf / SolidsCatalog / FluidsCatalog / MaterialDatabase.cpp.
add_executable(simall_materials_tests
    unit/materials/test_polynomial_fit.cpp
    unit/materials/test_tabulated_property.cpp
    unit/materials/test_anisotropic_conductivity.cpp
    unit/materials/test_catalogs.cpp
    unit/materials/test_janaf.cpp)
target_link_libraries(simall_materials_tests PRIVATE
    SimAll::CompilerFlags SimAll::Materials SimAll::Core SimAll::Utilities
    Catch2::Catch2WithMain)

# Solver unit tests (Week 5). Links SimAll::Solver for SIMPLE / PISO /
# SIMPLEC / Coupled / TFQMR / SSOR / FCT / Rhie-Chow / dual-time / overset.
add_executable(simall_solver_tests
    unit/solver/test_ssor_preconditioner.cpp
    unit/solver/test_tfqmr_solver.cpp
    unit/solver/test_fct_zalesak.cpp)
target_link_libraries(simall_solver_tests PRIVATE
    SimAll::CompilerFlags SimAll::Solver SimAll::Meshing SimAll::Materials
    SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Turbulence unit tests (Week 6). Links SimAll::Turbulence for RNG / Realizable /
# k-ω Std / dynamic Smag / SA-IDDES / SAS-SST / LRR-RSM / k-kL-ω.
add_executable(simall_turbulence_tests
    unit/turbulence/test_turbulence_models.cpp)
target_link_libraries(simall_turbulence_tests PRIVATE
    SimAll::CompilerFlags SimAll::Turbulence SimAll::Solver SimAll::Meshing
    SimAll::Materials SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Multiphase unit tests (Week 7). Links SimAll::Multiphase for Kunz/Zwart
# cavitation, MixtureModel, DriftFlux, EulerianMultifluid, IsoAdvector,
# RpiWallBoiling.
add_executable(simall_multiphase_tests
    unit/multiphase/test_multiphase_models.cpp)
target_link_libraries(simall_multiphase_tests PRIVATE
    SimAll::CompilerFlags SimAll::Multiphase SimAll::Solver SimAll::Meshing
    SimAll::Materials SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Combustion unit tests (Week 8). Links SimAll::Combustion for Chemkin/NASA
# parsers, FGM, SLFM, transported PDF, NOx (thermal/prompt/fuel), and the
# two soot models (Moss-Brookes + MOMIC).
add_executable(simall_combustion_tests
    unit/combustion/test_combustion_models.cpp
    unit/combustion/test_cantera_parser.cpp)
target_link_libraries(simall_combustion_tests PRIVATE
    SimAll::CompilerFlags SimAll::Combustion SimAll::Solver SimAll::Meshing
    SimAll::Materials SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Zones unit tests — persistent NamedSelections + lumped FanHexZone.
add_executable(simall_zones_tests
    unit/zones/test_selection_manager.cpp)
target_link_libraries(simall_zones_tests PRIVATE
    SimAll::CompilerFlags SimAll::Zones SimAll::Core SimAll::Utilities
    Catch2::Catch2WithMain)

# Particles unit tests (Week 9). Links SimAll::Particles for the eight
# new modules: Rosin-Rammler / blob / cone injectors, evaporation, heating,
# coalescence, wall film, charged-particle field.
add_executable(simall_particles_tests
    unit/particles/test_particles_models.cpp)
target_link_libraries(simall_particles_tests PRIVATE
    SimAll::CompilerFlags SimAll::Particles SimAll::Solver SimAll::Meshing
    SimAll::Materials SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Radiation unit tests (Week 10). Links SimAll::Radiation for Rosseland,
# Monte-Carlo RTE, surface-to-surface, solar load, WSGG spectral model.
add_executable(simall_radiation_tests
    unit/radiation/test_radiation_models.cpp)
target_link_libraries(simall_radiation_tests PRIVATE
    SimAll::CompilerFlags SimAll::Radiation SimAll::Solver SimAll::Meshing
    SimAll::Materials SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Meshing-advanced unit tests (Week 11). Links SimAll::Meshing for the
# advancing-front surface, hex sweep, multi-block hex, polyhedral
# agglomeration, boundary recovery, medial-axis decomposition, METIS /
# Scotch partitioners, and AMR adaptor modules.
add_executable(simall_meshing_advanced_tests
    unit/meshing/test_meshing_advanced.cpp)
target_link_libraries(simall_meshing_advanced_tests PRIVATE
    SimAll::CompilerFlags SimAll::Meshing SimAll::Solver SimAll::Materials
    SimAll::Core SimAll::Utilities SimAll::Cad Catch2::Catch2WithMain)

# Parallel / HPC unit tests (Week 12). MpiContext, GhostExchange,
# DomainPartition, FieldReducer, LoadBalancer, NumaPinning.  All
# exercise the serial fall-back path (no MPI required at test time).
add_executable(simall_parallel_tests
    unit/parallel/test_parallel_models.cpp)
target_link_libraries(simall_parallel_tests PRIVATE
    SimAll::CompilerFlags SimAll::Parallel SimAll::Meshing
    SimAll::Core SimAll::Utilities SimAll::Cad Catch2::Catch2WithMain)

# GPU subsystem unit tests (Week 13). CudaContext, DeviceMemoryPool,
# StreamScheduler, HostDeviceMirror, GpuKernels, CgGpu, BicgstabGpu,
# RenderingBridge.  Exercises the SIMALL_NO_CUDA fall-back so the
# whole suite runs on CI nodes without a GPU.
add_executable(simall_gpu_tests
    unit/gpu/test_gpu_models.cpp)
target_link_libraries(simall_gpu_tests PRIVATE
    SimAll::CompilerFlags SimAll::Gpu
    SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Visualization & post-processing engine (Week 14). SceneGraph, ActorRegistry,
# CameraController, all eight filters, three overlays, two recorders, and the
# picking bridge.  Links SimAll::VisualizationCore (the VTK-free half of the
# subsystem) so the suite runs headless.
add_executable(simall_visualization_tests
    unit/visualization/test_visualization_models.cpp)
target_link_libraries(simall_visualization_tests PRIVATE
    SimAll::CompilerFlags SimAll::VisualizationCore
    SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# GUI Core (Week 15).  PropertyDescriptor validation, DockPerspective text
# round-trip, PanelRegistry factory dispatch, ThemePaths discovery — all
# Qt-free, runs on any CI worker.
add_executable(simall_gui_core_tests
    unit/gui_core/test_gui_core_models.cpp)
target_link_libraries(simall_gui_core_tests PRIVATE
    SimAll::CompilerFlags SimAll::GuiCore
    SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Scripting (Week 16). Repl continuation analysis, ScriptContext bag,
# CommandBus dispatch, MacroRecorder serialization, UdfHost registry.
# Python is NOT required: the stub branch of PyBindings.cpp keeps the
# library link-complete on hosts without pybind11.
add_executable(simall_scripting_tests
    unit/scripting/test_scripting_models.cpp)
target_link_libraries(simall_scripting_tests PRIVATE
    SimAll::CompilerFlags SimAll::Scripting
    SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# IO subsystem (Week 17).  Fluent .msh + Plot3D readers, CGNS native
# round-trip, Tecplot ASCII writer, native SRS result store, project
# graph hash-commit, .simall chunked container, APNG movie encoder.
# All backends ship a native fallback so the suite runs without libhdf5,
# libcgns or ffmpeg on the host.
add_executable(simall_io_tests
    unit/io/test_io_models.cpp)
target_link_libraries(simall_io_tests PRIVATE
    SimAll::CompilerFlags SimAll::Io
    SimAll::Core SimAll::Utilities Catch2::Catch2WithMain)

# Long-tail modules (Week 18).  Adjoint (discrete GMRES on Aᵀ + continuous
# shape sensitivity), Optimization (MMA, CMA-ES, LHS / Sobol / full-factorial
# DOE, Kriging surrogate), Morphing (FFD lattice), Dynamics (6DOF coupler),
# Rotating (sliding-mesh driver), AMR (orchestrator), IBM (direct-forcing),
# EMag (A-V tet4 assembler), Acoustics (Curle dipole), ROM (operator
# inference + ECSW hyper-reduction).
add_executable(simall_longtail_tests
    unit/longtail/test_longtail_modules.cpp)
target_link_libraries(simall_longtail_tests PRIVATE
    SimAll::CompilerFlags
    SimAll::Adjoint SimAll::Optimization SimAll::Morphing
    SimAll::Dynamics SimAll::Rotating  SimAll::Amr
    SimAll::Ibm     SimAll::Emag      SimAll::Acoustics
    SimAll::Rom
    SimAll::Solver  SimAll::Meshing   SimAll::Materials
    SimAll::Core    SimAll::Utilities
    Catch2::Catch2WithMain)

# Regression / verification / performance / golden suites (Week 19).  These
# do not link any SimAll subsystem library — they are header-only reference
# packs whose role is to lock in the canonical benchmark data and the FD
# Laplacian convergence rate.  Each executable carries a Catch2 tag that we
# expose to ctest via labels for selective CI execution.
add_executable(simall_regression_tests
    regression/test_regression_suite.cpp)
target_include_directories(simall_regression_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simall_regression_tests PRIVATE
    SimAll::CompilerFlags Catch2::Catch2WithMain)

add_executable(simall_verification_tests
    verification/test_mms.cpp)
target_include_directories(simall_verification_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simall_verification_tests PRIVATE
    SimAll::CompilerFlags Catch2::Catch2WithMain)

add_executable(simall_performance_tests
    performance/test_perf_harness.cpp)
target_include_directories(simall_performance_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simall_performance_tests PRIVATE
    SimAll::CompilerFlags Catch2::Catch2WithMain)

add_executable(simall_golden_tests
    golden/test_golden_hashes.cpp)
target_include_directories(simall_golden_tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(simall_golden_tests PRIVATE
    SimAll::CompilerFlags Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(simall_core_tests)
catch_discover_tests(simall_materials_tests)
catch_discover_tests(simall_solver_tests)
catch_discover_tests(simall_turbulence_tests)
catch_discover_tests(simall_multiphase_tests)
catch_discover_tests(simall_combustion_tests)
catch_discover_tests(simall_zones_tests)
catch_discover_tests(simall_particles_tests)
catch_discover_tests(simall_radiation_tests)
catch_discover_tests(simall_meshing_advanced_tests)
catch_discover_tests(simall_parallel_tests)
catch_discover_tests(simall_gpu_tests)
catch_discover_tests(simall_visualization_tests)
catch_discover_tests(simall_gui_core_tests)
catch_discover_tests(simall_scripting_tests)
catch_discover_tests(simall_io_tests)
catch_discover_tests(simall_longtail_tests)
catch_discover_tests(simall_regression_tests
    PROPERTIES LABELS "regression")
catch_discover_tests(simall_verification_tests
    PROPERTIES LABELS "verification")
catch_discover_tests(simall_performance_tests
    PROPERTIES LABELS "performance")
catch_discover_tests(simall_golden_tests
    PROPERTIES LABELS "golden")
