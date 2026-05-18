simall_add_module(gpu
    SOURCES
        Device.cpp
        CudaContext.cpp
        DeviceMemoryPool.cpp
        StreamScheduler.cpp
        GpuKernels.cpp
        CgGpu.cpp
        BicgstabGpu.cpp
        RenderingBridge.cpp
    HEADERS
        Device.hpp
        CudaContext.hpp
        DeviceMemoryPool.hpp
        StreamScheduler.hpp
        HostDeviceMirror.hpp
        GpuKernels.hpp
        CgGpu.hpp
        BicgstabGpu.hpp
        RenderingBridge.hpp
    PUBLIC_DEPS SimAll::Core SimAll::Utilities)

# -----------------------------------------------------------------------------
# Optional CUDA wiring.  When SIMALL_ENABLE_CUDA is ON and the CUDA toolkit
# was successfully located by simall_find_optional_dependencies(), expose the
# SIMALL_HAVE_CUDA compile definition and link cudart + cuBLAS + cuSPARSE so
# the device paths in CudaContext / DeviceMemoryPool / GpuKernels become
# active.  Otherwise the module compiles every translation unit under the
# host fall-back and SIMALL_NO_CUDA semantics apply transparently.
# -----------------------------------------------------------------------------
if(SIMALL_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    find_package(CUDAToolkit QUIET)
    if(CUDAToolkit_FOUND)
        target_compile_definitions(simall_gpu PUBLIC SIMALL_HAVE_CUDA)
        target_link_libraries(simall_gpu PUBLIC
            CUDA::cudart CUDA::cublas CUDA::cusparse)
        message(STATUS "  SimAll::Gpu  : CUDA toolkit "
            "${CUDAToolkit_VERSION} (cudart + cublas + cusparse)")
    else()
        message(WARNING "SIMALL_ENABLE_CUDA=ON but CUDAToolkit not found; "
            "SimAll::Gpu builds the serial fall-back.")
    endif()
endif()
