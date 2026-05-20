add_executable(simall_beta WIN32 main.cpp)
target_link_libraries(simall_beta PRIVATE
    SimAll::CompilerFlags
    SimAll::Core
    SimAll::Io
    SimAll::Parallel
    SimAll::Cad
    SimAll::Meshing
    SimAll::Solver
    SimAll::Turbulence
    SimAll::HeatTransfer
    SimAll::Radiation
    SimAll::Combustion
    SimAll::Multiphase
    SimAll::Visualization
    SimAll::Gui
    Qt6::Core Qt6::Gui Qt6::Widgets)

vtk_module_autoinit(TARGETS simall_beta MODULES ${VTK_LIBRARIES})

set_target_properties(simall_beta PROPERTIES
    FOLDER "SimAll/Applications"
    WIN32_EXECUTABLE TRUE)

# Install the GUI exe into bin/ (runtime component) and bundle every
# third-party DLL it transitively needs (Qt6 via windeployqt, the rest via
# file(GET_RUNTIME_DEPENDENCIES)) so the resulting MSI / NSIS installer is
# self-contained.  See cmake/SimAllInstall.cmake.
install(TARGETS simall_beta
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    COMPONENT runtime)
simall_install_runtime_deps(simall_beta)
