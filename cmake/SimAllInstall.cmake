# =============================================================================
# SimAllInstall.cmake
#
# Provides:
#   simall_install_module(<short_name>)
#       Install a simall_<short_name> target (created by simall_add_module)
#       under the SimAllTargets export set, including its public headers.
#
#   simall_export_targets()
#       Emit the SimAllConfig.cmake / SimAllConfigVersion.cmake package
#       files so downstream projects can `find_package(SimAll CONFIG)`.
#
#   simall_setup_cpack()
#       Configure NSIS (Windows), DEB (Linux), and TGZ generators.
# =============================================================================

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(SIMALL_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/SimAll"
    CACHE STRING "CMake package install dir")

# -----------------------------------------------------------------------------
function(simall_install_module name)
    set(_target "simall_${name}")
    if(NOT TARGET ${_target})
        message(WARNING "simall_install_module: no such target '${_target}'")
        return()
    endif()

    install(TARGETS ${_target}
        EXPORT  SimAllTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
                COMPONENT runtime
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                COMPONENT runtime
                NAMELINK_COMPONENT devel
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                COMPONENT devel
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    # Install the module's public headers preserving the source layout
    # (consumers use `#include "<name>/Foo.hpp"` both in-build and after install).
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/${name}/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${name}
        COMPONENT   devel
        FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "CMakeFiles" EXCLUDE
        PATTERN ".git"       EXCLUDE)
endfunction()

# -----------------------------------------------------------------------------
function(simall_export_targets)
    # Compiler-flag interface target is part of the export set.
    if(TARGET SimAll_CompilerFlags)
        install(TARGETS SimAll_CompilerFlags EXPORT SimAllTargets)
    endif()

    install(EXPORT SimAllTargets
        FILE        SimAllTargets.cmake
        NAMESPACE   SimAll::
        DESTINATION ${SIMALL_INSTALL_CMAKEDIR}
        COMPONENT   devel)

    write_basic_package_version_file(
        "${CMAKE_BINARY_DIR}/SimAllConfigVersion.cmake"
        VERSION       ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion)

    # SimAllConfig.cmake.in is generated below if it does not exist.
    set(_cfg_in "${CMAKE_SOURCE_DIR}/cmake/SimAllConfig.cmake.in")
    if(NOT EXISTS "${_cfg_in}")
        file(WRITE "${_cfg_in}"
"@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(Eigen3 3.4 REQUIRED)

include(\"\${CMAKE_CURRENT_LIST_DIR}/SimAllTargets.cmake\")

check_required_components(SimAll)
")
    endif()

    configure_package_config_file(
        "${_cfg_in}"
        "${CMAKE_BINARY_DIR}/SimAllConfig.cmake"
        INSTALL_DESTINATION ${SIMALL_INSTALL_CMAKEDIR})

    install(FILES
        "${CMAKE_BINARY_DIR}/SimAllConfig.cmake"
        "${CMAKE_BINARY_DIR}/SimAllConfigVersion.cmake"
        DESTINATION ${SIMALL_INSTALL_CMAKEDIR}
        COMPONENT   devel)

    # Resource bundles (always installed when present).
    if(EXISTS "${CMAKE_SOURCE_DIR}/shaders")
        install(DIRECTORY "${CMAKE_SOURCE_DIR}/shaders/"
            DESTINATION ${CMAKE_INSTALL_DATADIR}/simall/shaders
            COMPONENT   runtime
            FILES_MATCHING PATTERN "*.vert" PATTERN "*.frag"
                           PATTERN "*.geom" PATTERN "*.glsl")
    endif()
    if(EXISTS "${CMAKE_SOURCE_DIR}/resources")
        install(DIRECTORY "${CMAKE_SOURCE_DIR}/resources/"
            DESTINATION ${CMAKE_INSTALL_DATADIR}/simall/resources
            COMPONENT   runtime
            FILES_MATCHING PATTERN "*.qss" PATTERN "*.svg" PATTERN "*.qrc")
    endif()
    if(EXISTS "${CMAKE_SOURCE_DIR}/docs")
        install(DIRECTORY "${CMAKE_SOURCE_DIR}/docs/"
            DESTINATION ${CMAKE_INSTALL_DOCDIR}
            COMPONENT   docs
            FILES_MATCHING PATTERN "*.md" PATTERN "*.puml" PATTERN "*.yml")
    endif()
endfunction()

# -----------------------------------------------------------------------------
function(simall_setup_cpack)
    set(CPACK_PACKAGE_NAME                "SimAll-Beta"        PARENT_SCOPE)
    set(CPACK_PACKAGE_VENDOR              "SimAll Project"     PARENT_SCOPE)
    set(CPACK_PACKAGE_VERSION_MAJOR       ${PROJECT_VERSION_MAJOR} PARENT_SCOPE)
    set(CPACK_PACKAGE_VERSION_MINOR       ${PROJECT_VERSION_MINOR} PARENT_SCOPE)
    set(CPACK_PACKAGE_VERSION_PATCH       ${PROJECT_VERSION_PATCH} PARENT_SCOPE)
    set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}"  PARENT_SCOPE)
    set(CPACK_PACKAGE_INSTALL_DIRECTORY   "SimAll Beta"        PARENT_SCOPE)
    set(CPACK_PACKAGE_CONTACT             "simall@example.org" PARENT_SCOPE)

    # Week 20 — embed the top-level legal / changelog assets into every
    # generator so installers can render the license-agreement page and
    # the documentation component carries a release-notes file.
    if(EXISTS "${CMAKE_SOURCE_DIR}/LICENSE")
        set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE" PARENT_SCOPE)
    endif()
    if(EXISTS "${CMAKE_SOURCE_DIR}/README.md")
        set(CPACK_RESOURCE_FILE_README  "${CMAKE_SOURCE_DIR}/README.md"  PARENT_SCOPE)
    endif()
    if(EXISTS "${CMAKE_SOURCE_DIR}/CHANGELOG.md")
        install(FILES
            "${CMAKE_SOURCE_DIR}/LICENSE"
            "${CMAKE_SOURCE_DIR}/CHANGELOG.md"
            "${CMAKE_SOURCE_DIR}/RELEASE_NOTES_1.0-rc1.md"
            DESTINATION ${CMAKE_INSTALL_DOCDIR}
            COMPONENT   docs
            OPTIONAL)
    endif()

    set(CPACK_COMPONENTS_ALL              runtime devel docs   PARENT_SCOPE)
    set(CPACK_COMPONENT_RUNTIME_DISPLAY_NAME "Runtime"          PARENT_SCOPE)
    set(CPACK_COMPONENT_DEVEL_DISPLAY_NAME   "SDK / Development" PARENT_SCOPE)
    set(CPACK_COMPONENT_DOCS_DISPLAY_NAME    "Documentation"    PARENT_SCOPE)
    set(CPACK_COMPONENT_DEVEL_DEPENDS        runtime            PARENT_SCOPE)

    if(WIN32)
        set(CPACK_GENERATOR "NSIS;ZIP" PARENT_SCOPE)
        set(CPACK_NSIS_PACKAGE_NAME       "SimAll Beta"       PARENT_SCOPE)
        set(CPACK_NSIS_DISPLAY_NAME       "SimAll Beta"       PARENT_SCOPE)
        set(CPACK_NSIS_HELP_LINK          "https://example.org/simall" PARENT_SCOPE)
        set(CPACK_NSIS_URL_INFO_ABOUT     "https://example.org/simall" PARENT_SCOPE)
        set(CPACK_NSIS_MODIFY_PATH        ON                  PARENT_SCOPE)
        set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON     PARENT_SCOPE)
    elseif(APPLE)
        set(CPACK_GENERATOR "TGZ;DragNDrop" PARENT_SCOPE)
    else()
        set(CPACK_GENERATOR "TGZ;DEB" PARENT_SCOPE)
        set(CPACK_DEBIAN_PACKAGE_MAINTAINER  "SimAll Project <simall@example.org>" PARENT_SCOPE)
        set(CPACK_DEBIAN_PACKAGE_SECTION     "science"             PARENT_SCOPE)
        set(CPACK_DEBIAN_PACKAGE_PRIORITY    "optional"            PARENT_SCOPE)
        set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS   ON                    PARENT_SCOPE)
        set(CPACK_DEBIAN_COMPONENT_INSTALL   ON                    PARENT_SCOPE)
    endif()

    include(CPack)
endfunction()
