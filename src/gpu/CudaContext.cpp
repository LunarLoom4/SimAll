// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/CudaContext.cpp
// Phase  : 18.1 — CUDA device context (W13).
// =============================================================================
#include "gpu/CudaContext.hpp"

#include "core/Logger.hpp"

#include <cstdlib>
#include <mutex>

#ifdef SIMALL_HAVE_CUDA
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#endif

namespace simall::gpu
{

// -----------------------------------------------------------------------------
// Internal state.  All CUDA-specific members are gated; the serial fall-back
// keeps only the cached DeviceProperties object.
// -----------------------------------------------------------------------------
struct CudaContext::Impl
{
    int deviceId = -1;
    int deviceCount = 0;
    DeviceProperties props{};
    std::string lastErrorMessage = "ok";

#ifdef SIMALL_HAVE_CUDA
    cublasHandle_t cublas = nullptr;
    cusparseHandle_t cusparse = nullptr;
#endif
};

namespace
{

#ifdef SIMALL_HAVE_CUDA
/// Translate a cudaError_t into a SimAll-flavoured human string.
inline std::string cuda_err_str(cudaError_t e)
{
    if (e == cudaSuccess)
        return "ok";
    return std::string(cudaGetErrorName(e)) + ": " + cudaGetErrorString(e);
}

/// Fill `out` from cudaDeviceProp + memInfo for the currently selected device.
DeviceProperties query_device(int id)
{
    DeviceProperties out{};
    out.device_id = id;
    cudaDeviceProp p{};
    if (cudaGetDeviceProperties(&p, id) != cudaSuccess)
        return out;
    out.name = p.name;
    out.compute_capability_major = p.major;
    out.compute_capability_minor = p.minor;
    out.sm_count = p.multiProcessorCount;
    out.warp_size = p.warpSize;
    out.max_threads_per_block = p.maxThreadsPerBlock;
    out.total_memory_bytes = p.totalGlobalMem;
    out.shared_memory_per_block = p.sharedMemPerBlock;
    out.unified_addressing = p.unifiedAddressing != 0;
    out.managed_memory = p.managedMemory != 0;
    std::size_t free_b = 0, total_b = 0;
    if (cudaMemGetInfo(&free_b, &total_b) == cudaSuccess) {
        out.free_memory_bytes = free_b;
        // total_memory_bytes already populated from cudaDeviceProp.
    }
    return out;
}
#endif

} // namespace

// -----------------------------------------------------------------------------
// Singleton accessor
// -----------------------------------------------------------------------------
CudaContext& CudaContext::current()
{
    static CudaContext instance;
    return instance;
}

// -----------------------------------------------------------------------------
// Construction (lazy via static instance)
// -----------------------------------------------------------------------------
CudaContext::CudaContext() : impl_(std::make_unique<Impl>())
{
#ifdef SIMALL_HAVE_CUDA
    int n = 0;
    const cudaError_t qErr = cudaGetDeviceCount(&n);
    if (qErr != cudaSuccess || n == 0) {
        impl_->deviceCount = 0;
        impl_->deviceId = -1;
        impl_->lastErrorMessage = cuda_err_str(qErr);
        SIMALL_LOG_INFO("Gpu",
                        "CUDA toolkit linked but no usable device "
                        "(device_count=",
                        n,
                        ", err=",
                        impl_->lastErrorMessage,
                        ").  "
                        "Falling back to host kernels.");
        return;
    }
    impl_->deviceCount = n;
    impl_->deviceId = 0;
    cudaSetDevice(0);
    impl_->props = query_device(0);

    if (cublasCreate(&impl_->cublas) != CUBLAS_STATUS_SUCCESS) {
        impl_->cublas = nullptr;
        impl_->lastErrorMessage = "cublasCreate failed";
    } else {
        cublasSetPointerMode(impl_->cublas, CUBLAS_POINTER_MODE_HOST);
    }
    if (cusparseCreate(&impl_->cusparse) != CUSPARSE_STATUS_SUCCESS) {
        impl_->cusparse = nullptr;
        if (impl_->lastErrorMessage == "ok")
            impl_->lastErrorMessage = "cusparseCreate failed";
    }
    SIMALL_LOG_INFO("Gpu",
                    "CUDA context ready — device=",
                    impl_->props.name,
                    " sm=",
                    impl_->props.sm_count,
                    " cc=",
                    impl_->props.compute_capability_major,
                    ".",
                    impl_->props.compute_capability_minor,
                    " mem=",
                    (impl_->props.total_memory_bytes >> 20),
                    " MiB");
#else
    // Serial fall-back: device_count=0, properties default-initialised.
    SIMALL_LOG_INFO("Gpu", "SIMALL_NO_CUDA build — host kernels only.");
#endif
}

CudaContext::~CudaContext()
{
#ifdef SIMALL_HAVE_CUDA
    if (impl_->cusparse)
        cusparseDestroy(impl_->cusparse);
    if (impl_->cublas)
        cublasDestroy(impl_->cublas);
    // No cudaDeviceReset() — leave teardown to the runtime so other
    // singletons (e.g. cuBLAS in plugins) can still run.
#endif
}

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------
bool is_cuda_available() noexcept
{
    return CudaContext::current().device_count() > 0;
}

int CudaContext::device_count() const noexcept
{
    return impl_->deviceCount;
}
int CudaContext::device_id() const noexcept
{
    return impl_->deviceId;
}

int CudaContext::select_device(int id)
{
    const int prev = impl_->deviceId;
#ifdef SIMALL_HAVE_CUDA
    if (id < 0 || id >= impl_->deviceCount)
        return prev;
    if (id == prev)
        return prev;
    if (impl_->cusparse) {
        cusparseDestroy(impl_->cusparse);
        impl_->cusparse = nullptr;
    }
    if (impl_->cublas) {
        cublasDestroy(impl_->cublas);
        impl_->cublas = nullptr;
    }
    cudaSetDevice(id);
    impl_->deviceId = id;
    impl_->props = query_device(id);
    cublasCreate(&impl_->cublas);
    cublasSetPointerMode(impl_->cublas, CUBLAS_POINTER_MODE_HOST);
    cusparseCreate(&impl_->cusparse);
#else
    (void) id;
#endif
    return prev;
}

DeviceProperties CudaContext::properties() const
{
    return impl_->props;
}

void CudaContext::refresh_memory_info()
{
#ifdef SIMALL_HAVE_CUDA
    if (impl_->deviceCount == 0)
        return;
    std::size_t free_b = 0, total_b = 0;
    if (cudaMemGetInfo(&free_b, &total_b) == cudaSuccess) {
        impl_->props.free_memory_bytes = free_b;
        impl_->props.total_memory_bytes = total_b;
    }
#endif
}

void CudaContext::synchronize() noexcept
{
#ifdef SIMALL_HAVE_CUDA
    if (impl_->deviceCount > 0)
        cudaDeviceSynchronize();
#endif
}

std::string CudaContext::last_error() const
{
#ifdef SIMALL_HAVE_CUDA
    if (impl_->deviceCount == 0)
        return impl_->lastErrorMessage;
    const cudaError_t e = cudaPeekAtLastError();
    return cuda_err_str(e);
#else
    return impl_->lastErrorMessage;
#endif
}

void* CudaContext::cublas_handle() const noexcept
{
#ifdef SIMALL_HAVE_CUDA
    return static_cast<void*>(impl_->cublas);
#else
    return nullptr;
#endif
}

void* CudaContext::cusparse_handle() const noexcept
{
#ifdef SIMALL_HAVE_CUDA
    return static_cast<void*>(impl_->cusparse);
#else
    return nullptr;
#endif
}

} // namespace simall::gpu
