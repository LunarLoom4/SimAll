// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/Device.hpp
// Phase  : 18 (GPU ARCHITECTURE)
//
// Device-agnostic façade. Currently a CPU fall-back; CUDA kernels register
// themselves through `simall::gpu::register_backend()` when SIMALL_ENABLE_CUDA
// is on. Linear-algebra back-ends (CSR SpMV, GMRES) call `gpu::is_available()`
// and route accordingly.
// =============================================================================
#pragma once

#include <cstddef>
#include <string>

namespace simall::gpu {

struct DeviceInfo {
    std::string name;
    std::size_t total_memory_bytes;
    int         compute_capability_major;
    int         compute_capability_minor;
};

bool        is_available()    noexcept;
int         device_count()    noexcept;
DeviceInfo  query(int device) noexcept;
void        synchronize()     noexcept;

}  // namespace simall::gpu
