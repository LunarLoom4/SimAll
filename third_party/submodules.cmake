# ----------------------------------------------------------------------
# SimAll Beta — third-party fetch manifest.
# Included from the top-level CMakeLists.txt (after SIMALL_USE_FETCHCONTENT
# is honoured). Each block is guarded so the project also builds against
# a system-installed copy when one is present.
# ----------------------------------------------------------------------

include(FetchContent)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

# ---- Eigen -----------------------------------------------------------
# Prefer a packaged Eigen3 (e.g. from vcpkg). Fetching upstream Eigen and
# add_subdirectory'ing it pulls in eigen/blas/CMakeLists.txt, which calls
# enable_language(Fortran). On runners that have a non-functional gfortran
# on PATH (e.g. mingw gfortran with MSVC linker flags) this aborts the
# whole configure. We sidestep that entirely when a CONFIG package exists.
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.4 CONFIG QUIET)
endif()
if(NOT TARGET Eigen3::Eigen)
    FetchContent_Declare(eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
        GIT_SHALLOW    TRUE)
    set(EIGEN_BUILD_DOC      OFF CACHE BOOL "" FORCE)
    set(EIGEN_BUILD_TESTING  OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING        OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
endif()

# ---- fmt -------------------------------------------------------------
if(NOT TARGET fmt::fmt)
    FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG        10.2.1
        GIT_SHALLOW    TRUE)
    set(FMT_INSTALL ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(fmt)
endif()

# ---- spdlog ----------------------------------------------------------
if(NOT TARGET spdlog::spdlog)
    set(SPDLOG_FMT_EXTERNAL ON CACHE BOOL "" FORCE)
    set(SPDLOG_INSTALL      ON CACHE BOOL "" FORCE)
    FetchContent_Declare(spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.13.0
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(spdlog)
endif()

# ---- nlohmann/json ---------------------------------------------------
if(NOT TARGET nlohmann_json::nlohmann_json)
    set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
    set(JSON_Install    ON  CACHE BOOL "" FORCE)
    FetchContent_Declare(nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.11.3
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(nlohmann_json)
endif()

# ---- Catch2 (tests only) --------------------------------------------
if(SIMALL_BUILD_TESTS AND NOT TARGET Catch2::Catch2WithMain)
    FetchContent_Declare(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.5.4
        GIT_SHALLOW    TRUE)
    set(CATCH_INSTALL_DOCS OFF CACHE BOOL "" FORCE)
    set(CATCH_INSTALL_EXTRAS OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${Catch2_SOURCE_DIR}/extras)
endif()

# ---- pybind11 (scripting only) --------------------------------------
if(SIMALL_ENABLE_SCRIPTING AND NOT TARGET pybind11::pybind11)
    FetchContent_Declare(pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG        v2.12.0
        GIT_SHALLOW    TRUE)
    FetchContent_MakeAvailable(pybind11)
endif()

# ---- METIS (mesh partitioning) --------------------------------------
# METIS does not ship a modern CMake; we provide a thin wrapper under
# third_party/metis_wrapper/ that calls into the official sources.
# See third_party/README.md for licensing/versioning notes.
