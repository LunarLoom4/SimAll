# =============================================================================
# SignPackage.cmake -- CPack post-build hook
#
# Called automatically by CPack (via CPACK_POST_BUILD_SCRIPTS, set in
# SimAllInstall.cmake's simall_setup_cpack()) when the SIMALL_SIGN_IDENTITY
# env var (or CMake cache var) is defined.
#
# Signs:
#   1. The main GUI executable inside the staged install tree.
#   2. The final installer artefacts (*.msi, *.exe) produced by CPack.
#
# Conventions:
#   * Uses signtool.exe (Windows SDK).  If absent on PATH, the script logs
#     a warning and exits 0 -- the package still ships unsigned, which is
#     correct behaviour for dev / unsigned builds.
#   * Uses RFC 3161 timestamping against DigiCert (configurable via
#     SIMALL_SIGN_TIMESTAMP_URL).  Timestamping is REQUIRED -- without it
#     signatures expire with the cert.
#   * `signtool /a` auto-picks the cert from the current user's "My"
#     certificate store.  For HSM / cloud-stored EV certs (the modern
#     standard since June 2023), install the vendor's CSP / KSP first and
#     /a will Just Work.
#
# Variables CPack exposes when invoking this script:
#   CPACK_PACKAGE_FILES        -- list of paths to generated installers
#   CPACK_TEMPORARY_DIRECTORY  -- staged install tree (pre-packaging)
# =============================================================================

if(NOT WIN32)
    return()
endif()

# --- Resolve signtool ---------------------------------------------------------
find_program(SIGNTOOL_EXE
    NAMES signtool
    PATHS
        "$ENV{ProgramFiles\(x86\)}/Windows Kits/10/bin/x64"
        "$ENV{ProgramFiles\(x86\)}/Windows Kits/10/App Certification Kit"
    PATH_SUFFIXES "bin/x64" "x64"
    NO_DEFAULT_PATH)
if(NOT SIGNTOOL_EXE)
    find_program(SIGNTOOL_EXE NAMES signtool)
endif()

if(NOT SIGNTOOL_EXE)
    message(WARNING "SignPackage: signtool.exe not found -- skipping code signing. "
        "Install the Windows 10/11 SDK to enable.")
    return()
endif()

# --- Config -------------------------------------------------------------------
set(_ts_url "$ENV{SIMALL_SIGN_TIMESTAMP_URL}")
if(NOT _ts_url)
    set(_ts_url "http://timestamp.digicert.com")
endif()

set(_sign_args /tr "${_ts_url}" /td sha256 /fd sha256 /a)

# Optional: subject-name pinning for the cert (avoids picking wrong cert when
# multiple are in the store).  Set SIMALL_SIGN_SUBJECT="CN=Synopsys, Inc."
if(DEFINED ENV{SIMALL_SIGN_SUBJECT})
    list(APPEND _sign_args /n "$ENV{SIMALL_SIGN_SUBJECT}")
endif()

# --- 1. Sign binaries inside the staged install tree --------------------------
# (Belt-and-braces: most installers expect their payload already signed.)
file(GLOB_RECURSE _staged_exes
    "${CPACK_TEMPORARY_DIRECTORY}/*.exe"
    "${CPACK_TEMPORARY_DIRECTORY}/*.dll")
foreach(_bin IN LISTS _staged_exes)
    # Only sign our own outputs -- don't re-sign vendor DLLs (Qt/VTK/OCCT)
    # which are already signed by their publishers and would lose those sigs.
    get_filename_component(_name "${_bin}" NAME)
    if(_name MATCHES "^(simall_|SimAll)" )
        message(STATUS "SignPackage: signing staged ${_name}")
        execute_process(
            COMMAND "${SIGNTOOL_EXE}" sign ${_sign_args} "${_bin}"
            RESULT_VARIABLE _rc)
        if(NOT _rc EQUAL 0)
            message(WARNING "SignPackage: signtool returned ${_rc} for ${_bin}")
        endif()
    endif()
endforeach()

# --- 2. Sign the final installer artefacts ------------------------------------
foreach(_pkg IN LISTS CPACK_PACKAGE_FILES)
    if(_pkg MATCHES "\\.(msi|exe)$")
        message(STATUS "SignPackage: signing installer ${_pkg}")
        execute_process(
            COMMAND "${SIGNTOOL_EXE}" sign ${_sign_args} "${_pkg}"
            RESULT_VARIABLE _rc)
        if(NOT _rc EQUAL 0)
            message(WARNING "SignPackage: signtool returned ${_rc} for ${_pkg}")
        endif()
    endif()
endforeach()
