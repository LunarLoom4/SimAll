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
