# =============================================================================
# SimAllModule.cmake
# Helper for declaring SimAll subsystem libraries with a consistent layout:
#   simall_add_module(<short_name>
#       SOURCES   <list>
#       HEADERS   <list>
#       PUBLIC_DEPS  <targets>
#       PRIVATE_DEPS <targets>
#       [STATIC|SHARED|INTERFACE])
# Creates target `simall_<short_name>` and alias `SimAll::<ShortName>`.
# =============================================================================

function(simall_add_module name)
    set(options STATIC SHARED INTERFACE)
    set(one_value_args)
    set(multi_value_args SOURCES HEADERS PUBLIC_DEPS PRIVATE_DEPS)
    cmake_parse_arguments(ARG "${options}" "${one_value_args}" "${multi_value_args}" ${ARGN})

    include(GNUInstallDirs)
    set(_target "simall_${name}")

    if(ARG_INTERFACE)
        add_library(${_target} INTERFACE)
        target_include_directories(${_target} INTERFACE
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
        target_link_libraries(${_target} INTERFACE
            SimAll::CompilerFlags ${ARG_PUBLIC_DEPS})
    else()
        if(ARG_SHARED)
            add_library(${_target} SHARED ${ARG_SOURCES} ${ARG_HEADERS})
        else()
            add_library(${_target} STATIC ${ARG_SOURCES} ${ARG_HEADERS})
        endif()

        target_include_directories(${_target}
            PUBLIC  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
                    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
            PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

        target_link_libraries(${_target}
            PUBLIC  SimAll::CompilerFlags ${ARG_PUBLIC_DEPS}
            PRIVATE ${ARG_PRIVATE_DEPS})

        set_target_properties(${_target} PROPERTIES
            CXX_VISIBILITY_PRESET hidden
            VISIBILITY_INLINES_HIDDEN ON
            FOLDER "SimAll/Subsystems")
    endif()

    # Convert snake_case → CamelCase alias.
    string(REPLACE "_" ";" _parts "${name}")
    set(_alias "")
    foreach(p ${_parts})
        string(SUBSTRING ${p} 0 1 _first)
        string(TOUPPER ${_first} _first)
        string(SUBSTRING ${p} 1 -1 _rest)
        set(_alias "${_alias}${_first}${_rest}")
    endforeach()
    add_library(SimAll::${_alias} ALIAS ${_target})
endfunction()
