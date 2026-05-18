cmake_minimum_required(VERSION 3.22)

add_library(simall_csv_writer MODULE CsvWriterPlugin.cpp)
target_include_directories(simall_csv_writer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(simall_csv_writer PRIVATE SimAll::CompilerFlags)
set_target_properties(simall_csv_writer PROPERTIES PREFIX "")

install(TARGETS simall_csv_writer
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/simall/plugins COMPONENT runtime)
