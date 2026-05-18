cmake_minimum_required(VERSION 3.22)

add_library(simall_hello_solver MODULE HelloSolverPlugin.cpp)
target_include_directories(simall_hello_solver PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(simall_hello_solver PRIVATE SimAll::CompilerFlags)
set_target_properties(simall_hello_solver PROPERTIES PREFIX "")

install(TARGETS simall_hello_solver
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/simall/plugins COMPONENT runtime)
