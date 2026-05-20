// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/CgGpu.cpp
// Phase  : 18.6 — GPU Conjugate-Gradient solver (W13).
// =============================================================================
#include "gpu/CgGpu.hpp"

#include "core/Logger.hpp"
#include "gpu/CudaContext.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>

namespace simall::gpu
{

namespace
{
/// True when device pointer != host pointer (i.e. real CUDA path).  When
/// the mirror's owns_device() is false, "device" memory aliases host data
/// and we can read it directly without to_host().
inline bool dev(const HostDeviceMirror<double>& m)
{
    return m.owns_device();
}
} // namespace

// -----------------------------------------------------------------------------
CgGpu::CgGpu(GpuSolverConfig cfg) : cfg_(cfg) {}
CgGpu::~CgGpu() = default;

// -----------------------------------------------------------------------------
void CgGpu::set_matrix(const HostCsrView& A)
{
    n_ = A.n;
    nnz_ = A.nnz;

    dRowPtr_.resize(n_ + 1);
    dColIdx_.resize(nnz_);
    dValues_.resize(nnz_);

    std::memcpy(dRowPtr_.host_ptr(), A.rowPtr, (n_ + 1) * sizeof(int));
    std::memcpy(dColIdx_.host_ptr(), A.colIdx, nnz_ * sizeof(int));
    std::memcpy(dValues_.host_ptr(), A.values, nnz_ * sizeof(double));

    // Build preconditioner data on the host where appropriate.
    if (cfg_.preconditioner == GpuPreconditioner::Jacobi) {
        dDiag_.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            double d = 1.0;
            for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k)
                if (A.colIdx[k] == static_cast<int>(i)) {
                    d = A.values[k];
                    break;
                }
            dDiag_.host_ptr()[i] = d;
        }
        dDiag_.to_device();
    } else if (cfg_.preconditioner == GpuPreconditioner::ILU0) {
        dLU_.resize(nnz_);
        std::memcpy(dLU_.host_ptr(), A.values, nnz_ * sizeof(double));
        iluValid_ = ilu0_factor_inplace(n_, A.rowPtr, A.colIdx, dLU_.host_ptr());
        if (!iluValid_) {
            SIMALL_LOG_INFO("Gpu",
                            "ILU(0) factor failed (zero pivot); "
                            "falling back to Jacobi preconditioner.");
            cfg_.preconditioner = GpuPreconditioner::Jacobi;
            // recurse-equivalent: rebuild diagonal
            dDiag_.resize(n_);
            for (std::size_t i = 0; i < n_; ++i) {
                double d = 1.0;
                for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k)
                    if (A.colIdx[k] == static_cast<int>(i)) {
                        d = A.values[k];
                        break;
                    }
                dDiag_.host_ptr()[i] = d;
            }
            dDiag_.to_device();
        } else {
            dLU_.to_device();
        }
    }

    dRowPtr_.to_device();
    dColIdx_.to_device();
    dValues_.to_device();

    // Workspace
    dB_.resize(n_);
    dX_.resize(n_);
    dR_.resize(n_);
    dZ_.resize(n_);
    dP_.resize(n_);
    dQ_.resize(n_);
}

// -----------------------------------------------------------------------------
GpuSolverStats CgGpu::solve(const std::vector<double>& b, std::vector<double>& x)
{
    GpuSolverStats st{};
    const auto t0 = std::chrono::steady_clock::now();
    if (n_ == 0 || b.size() != n_ || x.size() != n_) {
        SIMALL_LOG_INFO("Gpu",
                        "CgGpu::solve called with size mismatch "
                        "(n=",
                        n_,
                        ", b=",
                        b.size(),
                        ", x=",
                        x.size(),
                        ")");
        return st;
    }

    // Upload b and x₀.
    std::memcpy(dB_.host_ptr(), b.data(), n_ * sizeof(double));
    std::memcpy(dX_.host_ptr(), x.data(), n_ * sizeof(double));
    dB_.to_device();
    dX_.to_device();

    auto apply_preconditioner = [&](void* r, void* z) {
        switch (cfg_.preconditioner) {
        case GpuPreconditioner::None:
            copy(r, z, n_);
            break;
        case GpuPreconditioner::Jacobi:
            jacobi_apply(dDiag_.device_ptr(), r, z, n_);
            break;
        case GpuPreconditioner::ILU0:
            ilu0_apply(n_,
                       static_cast<const int*>(dRowPtr_.device_ptr()),
                       static_cast<const int*>(dColIdx_.device_ptr()),
                       static_cast<const double*>(dLU_.device_ptr()),
                       static_cast<const double*>(r),
                       static_cast<double*>(z));
            break;
        }
    };

    // r = b - A x₀
    spmv_csr(n_,
             dRowPtr_.device_ptr(),
             dColIdx_.device_ptr(),
             dValues_.device_ptr(),
             dX_.device_ptr(),
             dR_.device_ptr());
    // r = -r + b  →  axpby(1, b, -1, r) then r ← r       (we want r = b - A x)
    // Simpler: scal(-1, r) then axpy(1, b, r)
    scal(-1.0, dR_.device_ptr(), n_);
    axpy(1.0, dB_.device_ptr(), dR_.device_ptr(), n_);

    const double bNorm = nrm2(dB_.device_ptr(), n_);
    const double bRef = (bNorm > 0.0) ? bNorm : 1.0;
    double rNorm = nrm2(dR_.device_ptr(), n_);
    st.initialResidual = rNorm;

    if (rNorm / bRef <= cfg_.tolerance) {
        st.converged = true;
        st.iterations = 0;
        st.finalResidual = rNorm;
        lastResidual_ = rNorm;
        // x already correct; copy to caller.
        dX_.to_host();
        std::memcpy(x.data(), dX_.host_ptr(), n_ * sizeof(double));
        return st;
    }

    // z = M⁻¹ r
    apply_preconditioner(dR_.device_ptr(), dZ_.device_ptr());
    // p = z
    copy(dZ_.device_ptr(), dP_.device_ptr(), n_);

    double rho_old = dot(dR_.device_ptr(), dZ_.device_ptr(), n_);

    const double dropTol =
        (cfg_.relativeOrders > 0.0) ? rNorm * std::pow(10.0, -cfg_.relativeOrders) : 0.0;

    for (int it = 1; it <= cfg_.maxIterations; ++it) {
        // q = A p
        spmv_csr(n_,
                 dRowPtr_.device_ptr(),
                 dColIdx_.device_ptr(),
                 dValues_.device_ptr(),
                 dP_.device_ptr(),
                 dQ_.device_ptr());

        const double pq = dot(dP_.device_ptr(), dQ_.device_ptr(), n_);
        if (std::abs(pq) < 1e-300) {
            SIMALL_LOG_INFO("Gpu", "CG breakdown: p·q ≈ 0 at it=", it);
            break;
        }
        const double alpha = rho_old / pq;

        axpy(alpha, dP_.device_ptr(), dX_.device_ptr(), n_);
        axpy(-alpha, dQ_.device_ptr(), dR_.device_ptr(), n_);

        rNorm = nrm2(dR_.device_ptr(), n_);
        if (rNorm / bRef <= cfg_.tolerance || (dropTol > 0.0 && rNorm <= dropTol)) {
            st.converged = true;
            st.iterations = it;
            st.finalResidual = rNorm;
            break;
        }

        apply_preconditioner(dR_.device_ptr(), dZ_.device_ptr());
        const double rho_new = dot(dR_.device_ptr(), dZ_.device_ptr(), n_);
        const double beta = rho_new / rho_old;
        xpby(dZ_.device_ptr(), beta, dP_.device_ptr(), n_);
        rho_old = rho_new;

        if (it == cfg_.maxIterations) {
            st.iterations = it;
            st.finalResidual = rNorm;
            st.converged = false;
        }
    }

    lastResidual_ = st.finalResidual;
    dX_.to_host();
    std::memcpy(x.data(), dX_.host_ptr(), n_ * sizeof(double));

    const auto t1 = std::chrono::steady_clock::now();
    st.timeSeconds = std::chrono::duration<double>(t1 - t0).count();
    SIMALL_LOG_INFO("Gpu",
                    "CgGpu: it=",
                    st.iterations,
                    " ||r||=",
                    st.finalResidual,
                    " converged=",
                    st.converged,
                    " time=",
                    st.timeSeconds,
                    "s");
    return st;
}

} // namespace simall::gpu
