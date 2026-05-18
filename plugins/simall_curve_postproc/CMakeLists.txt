cmake_minimum_required(VERSION 3.22)

add_library(simall_curve_postproc MODULE CurvePostprocPlugin.cpp)
target_include_directories(simall_curve_postproc PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(simall_curve_postproc PRIVATE SimAll::CompilerFlags)
set_target_properties(simall_curve_postproc PROPERTIES PREFIX "")

install(TARGETS simall_curve_postproc
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/simall/plugins COMPONENT runtime)
