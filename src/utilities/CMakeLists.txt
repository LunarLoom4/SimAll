# utilities — header-only foundation library (Phase 24)
simall_add_module(utilities INTERFACE)
target_include_directories(simall_utilities INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>)
