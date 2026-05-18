# =============================================================================
# SimAllDependencies.cmake
# Centralized find_package() calls. NEVER hardcode include paths in modules.
# =============================================================================

function(simall_find_required_dependencies)
    # ---- Qt6 (GUI) ----------------------------------------------------------
    find_package(Qt6 6.5 REQUIRED COMPONENTS
        Core Gui Widgets OpenGL OpenGLWidgets Concurrent Network Xml Svg)
    set(CMAKE_AUTOMOC ON PARENT_SCOPE)
    set(CMAKE_AUTORCC ON PARENT_SCOPE)
    set(CMAKE_AUTOUIC ON PARENT_SCOPE)

    # ---- VTK (Visualization) ------------------------------------------------
    find_package(VTK 9.2 REQUIRED COMPONENTS
        CommonCore CommonDataModel CommonExecutionModel
        FiltersCore FiltersGeneral FiltersGeometry FiltersSources
        FiltersFlowPaths FiltersHybrid FiltersModeling
        RenderingCore RenderingOpenGL2 RenderingAnnotation
        RenderingVolume RenderingVolumeOpenGL2 RenderingFreeType
        GUISupportQt InteractionStyle InteractionWidgets
        IOXML IOLegacy IOGeometry IOImage
        ImagingCore)

    # ---- OpenCASCADE (CAD) --------------------------------------------------
    find_package(OpenCASCADE 7.6 REQUIRED)

    # ---- Eigen3 (linear algebra) --------------------------------------------
    find_package(Eigen3 3.4 REQUIRED NO_MODULE)
endfunction()

function(simall_find_optional_dependencies)
    if(SIMALL_ENABLE_OPENMP)
        find_package(OpenMP)
        if(NOT OpenMP_CXX_FOUND)
            message(WARNING "OpenMP requested but not found; disabling.")
            set(SIMALL_ENABLE_OPENMP OFF PARENT_SCOPE)
        endif()
    endif()

    if(SIMALL_ENABLE_TBB)
        find_package(TBB CONFIG)
        if(NOT TBB_FOUND)
            message(WARNING "TBB requested but not found; disabling.")
            set(SIMALL_ENABLE_TBB OFF PARENT_SCOPE)
        endif()
    endif()

    if(SIMALL_ENABLE_MPI)
        find_package(MPI)
        if(NOT MPI_CXX_FOUND)
            message(WARNING "MPI requested but not found; disabling.")
            set(SIMALL_ENABLE_MPI OFF PARENT_SCOPE)
        endif()
    endif()

    if(SIMALL_ENABLE_CUDA)
        include(CheckLanguage)
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
            enable_language(CUDA)
            set(CMAKE_CUDA_STANDARD 17 PARENT_SCOPE)
        else()
            message(WARNING "CUDA requested but not available; disabling.")
            set(SIMALL_ENABLE_CUDA OFF PARENT_SCOPE)
        endif()
    endif()

    # ---- Python + pybind11 (scripting) --------------------------------------
    # SIMALL_ENABLE_SCRIPTING gates the build-in attempt; if Python or
    # pybind11 are missing we fall back to the no-op stub in PyBindings.cpp
    # so the rest of the codebase always sees the same API.
    set(SIMALL_HAVE_PYTHON OFF PARENT_SCOPE)
    if(SIMALL_ENABLE_SCRIPTING)
        find_package(Python3 3.9 COMPONENTS Interpreter Development)
        if(Python3_FOUND)
            find_package(pybind11 2.10 CONFIG QUIET)
            if(NOT pybind11_FOUND AND SIMALL_USE_FETCHCONTENT)
                include(FetchContent)
                FetchContent_Declare(pybind11
                    GIT_REPOSITORY https://github.com/pybind/pybind11.git
                    GIT_TAG        v2.12.0)
                FetchContent_MakeAvailable(pybind11)
            endif()
            if(TARGET pybind11::embed)
                set(SIMALL_HAVE_PYTHON ON PARENT_SCOPE)
                message(STATUS "SimAll scripting: Python ${Python3_VERSION} + pybind11 enabled")
            else()
                message(WARNING "SimAll scripting: pybind11 not found; Python bindings disabled")
            endif()
        else()
            message(WARNING "SimAll scripting: Python3 dev headers not found; Python bindings disabled")
        endif()
    endif()

    # ---- HDF5 (results store) -----------------------------------------------
    # Optional.  Native SRS fallback in Hdf5ResultStore.cpp is always built so
    # the test-suite runs on hosts without libhdf5.
    set(SIMALL_HAVE_HDF5 OFF PARENT_SCOPE)
    find_package(HDF5 QUIET COMPONENTS C)
    if(HDF5_FOUND)
        set(SIMALL_HAVE_HDF5 ON PARENT_SCOPE)
        message(STATUS "SimAll IO: HDF5 ${HDF5_VERSION} backend enabled")
    else()
        message(STATUS "SimAll IO: HDF5 not found; native SRS fallback only")
    endif()

    # ---- CGNS (mesh / solution interchange) ---------------------------------
    set(SIMALL_HAVE_CGNS OFF PARENT_SCOPE)
    find_package(CGNS QUIET)
    if(CGNS_FOUND)
        set(SIMALL_HAVE_CGNS ON PARENT_SCOPE)
        message(STATUS "SimAll IO: libcgns backend enabled")
    else()
        message(STATUS "SimAll IO: libcgns not found; CGNS-native fallback only")
    endif()
endfunction()
