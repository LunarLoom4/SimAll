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
# simall_install_runtime_deps(<target>)
#
# Bundles all third-party DLLs (Qt6, VTK, OCCT, HDF5, TBB, etc.) that the
# given Windows executable target depends on into the install tree's bin/
# folder, so the resulting MSI / NSIS installer is self-contained.
#
# Implementation:
#   1. Qt deployment uses windeployqt (ships with Qt6) which knows about
#      QML imports, translations, platform plugins, etc.  Found via the
#      Qt6::windeployqt imported tool target when available, falling back
#      to scanning Qt6Core's bin/ dir.
#   2. Everything else uses CMake's file(GET_RUNTIME_DEPENDENCIES) which
#      recursively walks PE-imports and copies resolved DLLs.  System DLLs
#      (api-ms-*, system32) are excluded.
# -----------------------------------------------------------------------------
function(simall_install_runtime_deps target)
    if(NOT WIN32)
        return()
    endif()
    if(NOT TARGET ${target})
        message(WARNING "simall_install_runtime_deps: target '${target}' does not exist; skipping.")
        return()
    endif()

    # --- Qt6 windeployqt ----------------------------------------------------
    if(TARGET Qt6::Core)
        get_target_property(_qt6_core_loc Qt6::Core LOCATION)
        get_filename_component(_qt6_bin_dir "${_qt6_core_loc}" DIRECTORY)
        find_program(SIMALL_WINDEPLOYQT
            NAMES windeployqt
            HINTS "${_qt6_bin_dir}" "${_qt6_bin_dir}/../tools/Qt6/bin"
            DOC "Qt6 windeployqt tool")
        if(SIMALL_WINDEPLOYQT)
            install(CODE "
                message(STATUS \"windeployqt: deploying Qt runtime for ${target}\")
                execute_process(COMMAND
                    \"${SIMALL_WINDEPLOYQT}\"
                    --no-translations
                    --no-system-d3d-compiler
                    --no-opengl-sw
                    --no-quick-import
                    --release
                    \"\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/$<TARGET_FILE_NAME:${target}>\")
            " COMPONENT runtime)
        else()
            message(WARNING "simall_install_runtime_deps: windeployqt not found near Qt6::Core; "
                "Qt DLLs will be picked up by GET_RUNTIME_DEPENDENCIES but plugins (platforms/, imageformats/, ...) "
                "will be missing.  Install Qt6 with the windeployqt tool to fix.")
        endif()
    endif()

    # --- Generic DLL closure via GET_RUNTIME_DEPENDENCIES -------------------
    # NB: we capture $<TARGET_FILE_NAME:${target}> at configure time so the
    # install script is a plain string by the time CMake runs it.  Variables
    # we want CMake to expand at *install* time (CMAKE_INSTALL_PREFIX, our
    # loop vars) are escaped with a leading backslash.
    set(_simall_exe_name "$<TARGET_FILE_NAME:${target}>")
    install(CODE "
        if(WIN32)
            file(GET_RUNTIME_DEPENDENCIES
                EXECUTABLES \"\${CMAKE_INSTALL_PREFIX}/bin/${_simall_exe_name}\"
                RESOLVED_DEPENDENCIES_VAR   _simall_resolved
                UNRESOLVED_DEPENDENCIES_VAR _simall_unresolved
                PRE_EXCLUDE_REGEXES  \"api-ms-\" \"ext-ms-\"
                POST_EXCLUDE_REGEXES \".*[\\\\\\\\/][Ss]ystem32[\\\\\\\\/].*\"
                                     \".*[\\\\\\\\/]SysWOW64[\\\\\\\\/].*\"
                                     \".*[\\\\\\\\/]WinSxS[\\\\\\\\/].*\")
            foreach(_dll IN LISTS _simall_resolved)
                file(INSTALL \"\${_dll}\"
                     DESTINATION \"\${CMAKE_INSTALL_PREFIX}/bin\"
                     FOLLOW_SYMLINK_CHAIN)
            endforeach()
            foreach(_u IN LISTS _simall_unresolved)
                message(STATUS \"simall_install_runtime_deps: unresolved '\${_u}' (delay-load/plugin; usually safe)\")
            endforeach()
        endif()
    " COMPONENT runtime)
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
        # NSIS + ZIP are the historical defaults; WIX (MSI) is added when
        # the WiX toolset (v3 candle.exe / v4 `wix.exe`) is on PATH, which
        # is what enterprise / Synopsys-style deployments expect.
        set(_simall_win_generators "NSIS;ZIP")
        find_program(_simall_wix4   NAMES wix.exe)
        find_program(_simall_wix3   NAMES candle.exe)
        if(_simall_wix4 OR _simall_wix3)
            list(PREPEND _simall_win_generators "WIX")
            # Pre-minted, never change across releases: Windows uses this to
            # recognise upgrades vs side-by-side installs.
            set(CPACK_WIX_UPGRADE_GUID
                "A2512926-3AAC-4CC2-BF04-46782C7E01FA" PARENT_SCOPE)
            set(CPACK_WIX_PRODUCT_GUID  "*"             PARENT_SCOPE)
            set(CPACK_WIX_PROGRAM_MENU_FOLDER "SimAll Beta" PARENT_SCOPE)
            if(EXISTS "${CMAKE_SOURCE_DIR}/assets/icon.ico")
                set(CPACK_WIX_PRODUCT_ICON
                    "${CMAKE_SOURCE_DIR}/assets/icon.ico" PARENT_SCOPE)
            endif()
            message(STATUS "SimAll: WiX toolset found -> MSI generator enabled.")
        else()
            message(STATUS "SimAll: WiX toolset not on PATH -> MSI generator skipped "
                "(install via `dotnet tool install --global wix` to enable).")
        endif()
        set(CPACK_GENERATOR "${_simall_win_generators}" PARENT_SCOPE)

        set(CPACK_NSIS_PACKAGE_NAME       "SimAll Beta"       PARENT_SCOPE)
        set(CPACK_NSIS_DISPLAY_NAME       "SimAll Beta"       PARENT_SCOPE)
        set(CPACK_NSIS_HELP_LINK          "https://example.org/simall" PARENT_SCOPE)
        set(CPACK_NSIS_URL_INFO_ABOUT     "https://example.org/simall" PARENT_SCOPE)
        set(CPACK_NSIS_MODIFY_PATH        ON                  PARENT_SCOPE)
        set(CPACK_NSIS_ENABLE_UNINSTALL_BEFORE_INSTALL ON     PARENT_SCOPE)

        # MSVC runtime DLLs (vcredist content) -- bundled into runtime component.
        set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_SKIP OFF        PARENT_SCOPE)
        set(CMAKE_INSTALL_SYSTEM_RUNTIME_COMPONENT runtime    PARENT_SCOPE)
        include(InstallRequiredSystemLibraries)

        # Post-build signing hook.  Activated only when SIMALL_SIGN_IDENTITY
        # is set (env var or cache var) -- typical CI flow.  Uses signtool's
        # `/a` switch to auto-pick the cert from the Windows cert store
        # (works with both file-based and HSM/cloud-stored certs).
        if(DEFINED ENV{SIMALL_SIGN_IDENTITY} OR SIMALL_SIGN_IDENTITY)
            set(CPACK_POST_BUILD_SCRIPTS
                "${CMAKE_SOURCE_DIR}/cmake/SignPackage.cmake" PARENT_SCOPE)
            message(STATUS "SimAll: code-signing enabled (SignPackage.cmake will run after CPack).")
        endif()
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
