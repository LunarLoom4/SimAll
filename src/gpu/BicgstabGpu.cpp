// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/BicgstabGpu.cpp
// Phase  : 18.6 — GPU BiCGSTAB solver (W13).
// =============================================================================
#include "gpu/BicgstabGpu.hpp"
#include "gpu/CudaContext.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>

namespace simall::gpu {

// -----------------------------------------------------------------------------
BicgstabGpu::BicgstabGpu(GpuSolverConfig cfg) : cfg_(cfg) {}
BicgstabGpu::~BicgstabGpu() = default;

// -----------------------------------------------------------------------------
void BicgstabGpu::set_matrix(const HostCsrView& A) {
    n_   = A.n;
    nnz_ = A.nnz;

    dRowPtr_.resize(n_ + 1);
    dColIdx_.resize(nnz_);
    dValues_.resize(nnz_);
    std::memcpy(dRowPtr_.host_ptr(), A.rowPtr, (n_ + 1) * sizeof(int));
    std::memcpy(dColIdx_.host_ptr(), A.colIdx, nnz_     * sizeof(int));
    std::memcpy(dValues_.host_ptr(), A.values, nnz_     * sizeof(double));
    dRowPtr_.to_device();
    dColIdx_.to_device();
    dValues_.to_device();

    if (cfg_.preconditioner == GpuPreconditioner::Jacobi) {
        dDiag_.resize(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            double d = 1.0;
            for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k)
                if (A.colIdx[k] == static_cast<int>(i)) { d = A.values[k]; break; }
            dDiag_.host_ptr()[i] = d;
        }
        dDiag_.to_device();
    } else if (cfg_.preconditioner == GpuPreconditioner::ILU0) {
        dLU_.resize(nnz_);
        std::memcpy(dLU_.host_ptr(), A.values, nnz_ * sizeof(double));
        iluValid_ = ilu0_factor_inplace(n_, A.rowPtr, A.colIdx, dLU_.host_ptr());
        if (!iluValid_) {
            SIMALL_LOG_INFO("Gpu", "BiCGSTAB: ILU(0) zero pivot; fallback Jacobi.");
            cfg_.preconditioner = GpuPreconditioner::Jacobi;
            dDiag_.resize(n_);
            for (std::size_t i = 0; i < n_; ++i) {
                double d = 1.0;
                for (int k = A.rowPtr[i]; k < A.rowPtr[i + 1]; ++k)
                    if (A.colIdx[k] == static_cast<int>(i)) { d = A.values[k]; break; }
                dDiag_.host_ptr()[i] = d;
            }
            dDiag_.to_device();
        } else {
            dLU_.to_device();
        }
    }

    dB_   .resize(n_);
    dX_   .resize(n_);
    dR_   .resize(n_);
    dRhat_.resize(n_);
    dP_   .resize(n_);
    dPhat_.resize(n_);
    dV_   .resize(n_);
    dS_   .resize(n_);
    dShat_.resize(n_);
    dT_   .resize(n_);
}

// -----------------------------------------------------------------------------
GpuSolverStats BicgstabGpu::solve(const std::vector<double>& b,
                                  std::vector<double>&       x)
{
    GpuSolverStats st{};
    const auto t0 = std::chrono::steady_clock::now();
    if (n_ == 0 || b.size() != n_ || x.size() != n_) {
        SIMALL_LOG_INFO("Gpu", "BicgstabGpu size mismatch: n=", n_,
            " b=", b.size(), " x=", x.size());
        return st;
    }
    std::memcpy(dB_.host_ptr(), b.data(), n_ * sizeof(double));
    std::memcpy(dX_.host_ptr(), x.data(), n_ * sizeof(double));
    dB_.to_device();
    dX_.to_device();

    auto apply_M = [&](void* in, void* out) {
        switch (cfg_.preconditioner) {
            case GpuPreconditioner::None:
                copy(in, out, n_);
                break;
            case GpuPreconditioner::Jacobi:
                jacobi_apply(dDiag_.device_ptr(), in, out, n_);
                break;
            case GpuPreconditioner::ILU0:
                ilu0_apply(n_,
                    static_cast<const int*>   (dRowPtr_.device_ptr()),
                    static_cast<const int*>   (dColIdx_.device_ptr()),
                    static_cast<const double*>(dLU_    .device_ptr()),
                    static_cast<const double*>(in),
                    static_cast<      double*>(out));
                break;
        }
    };

    // r = b - A x
    spmv_csr(n_, dRowPtr_.device_ptr(), dColIdx_.device_ptr(),
             dValues_.device_ptr(), dX_.device_ptr(), dR_.device_ptr());
    scal(-1.0, dR_.device_ptr(), n_);
    axpy(1.0, dB_.device_ptr(), dR_.device_ptr(), n_);

    // r̂ = r
    copy(dR_.device_ptr(), dRhat_.device_ptr(), n_);

    const double bNorm = nrm2(dB_.device_ptr(), n_);
    const double bRef  = (bNorm > 0.0) ? bNorm : 1.0;
    double rNorm = nrm2(dR_.device_ptr(), n_);
    st.initialResidual = rNorm;
    if (rNorm / bRef <= cfg_.tolerance) {
        st.converged = true;
        st.finalResidual = rNorm;
        lastResidual_ = rNorm;
        dX_.to_host();
        std::memcpy(x.data(), dX_.host_ptr(), n_ * sizeof(double));
        return st;
    }

    // p = v = 0
    dP_.zero_device();
    dV_.zero_device();

    double rho_old = 1.0, alpha = 1.0, omega = 1.0;

    const double dropTol = (cfg_.relativeOrders > 0.0)
        ? rNorm * std::pow(10.0, -cfg_.relativeOrders)
        : 0.0;

    for (int it = 1; it <= cfg_.maxIterations; ++it) {
        const double rho_new = dot(dRhat_.device_ptr(), dR_.device_ptr(), n_);
        if (std::abs(rho_new) < 1e-300) {
            SIMALL_LOG_INFO("Gpu", "BiCGSTAB breakdown ρ≈0 at it=", it);
            break;
        }
        const double beta = (rho_new / rho_old) * (alpha / omega);

        // p ← r + β (p − ω v)   = β p − β ω v + r
        // Implement as: p ← p − ω v  (axpy −ω)  then  p ← β p + r (axpby).
        axpy(-omega, dV_.device_ptr(), dP_.device_ptr(), n_);
        axpby(1.0, dR_.device_ptr(), beta, dP_.device_ptr(), n_);

        // p̂ = M⁻¹ p
        apply_M(dP_.device_ptr(), dPhat_.device_ptr());

        // v = A p̂
        spmv_csr(n_, dRowPtr_.device_ptr(), dColIdx_.device_ptr(),
                 dValues_.device_ptr(), dPhat_.device_ptr(), dV_.device_ptr());

        const double rhatv = dot(dRhat_.device_ptr(), dV_.device_ptr(), n_);
        if (std::abs(rhatv) < 1e-300) {
            SIMALL_LOG_INFO("Gpu", "BiCGSTAB breakdown (r̂·v ≈ 0) at it=", it);
            break;
        }
        alpha = rho_new / rhatv;

        // s = r − α v
        copy(dR_.device_ptr(), dS_.device_ptr(), n_);
        axpy(-alpha, dV_.device_ptr(), dS_.device_ptr(), n_);

        const double sNorm = nrm2(dS_.device_ptr(), n_);
        if (sNorm / bRef <= cfg_.tolerance ||
            (dropTol > 0.0 && sNorm <= dropTol)) {
            // x ← x + α p̂
            axpy(alpha, dPhat_.device_ptr(), dX_.device_ptr(), n_);
            st.converged     = true;
            st.iterations    = it;
            st.finalResidual = sNorm;
            break;
        }

        // ŝ = M⁻¹ s
        apply_M(dS_.device_ptr(), dShat_.device_ptr());

        // t = A ŝ
        spmv_csr(n_, dRowPtr_.device_ptr(), dColIdx_.device_ptr(),
                 dValues_.device_ptr(), dShat_.device_ptr(), dT_.device_ptr());

        const double tt = dot(dT_.device_ptr(), dT_.device_ptr(), n_);
        if (tt < 1e-300) {
            SIMALL_LOG_INFO("Gpu", "BiCGSTAB breakdown (t·t ≈ 0) at it=", it);
            break;
        }
        omega = dot(dT_.device_ptr(), dS_.device_ptr(), n_) / tt;
        if (std::abs(omega) < 1e-300) {
            SIMALL_LOG_INFO("Gpu", "BiCGSTAB stagnation ω≈0 at it=", it);
            break;
        }

        // x ← x + α p̂ + ω ŝ
        axpy(alpha, dPhat_.device_ptr(), dX_.device_ptr(), n_);
        axpy(omega, dShat_.device_ptr(), dX_.device_ptr(), n_);
        // r ← s − ω t
        copy(dS_.device_ptr(), dR_.device_ptr(), n_);
        axpy(-omega, dT_.device_ptr(), dR_.device_ptr(), n_);

        rNorm = nrm2(dR_.device_ptr(), n_);
        if (rNorm / bRef <= cfg_.tolerance ||
            (dropTol > 0.0 && rNorm <= dropTol)) {
            st.converged     = true;
            st.iterations    = it;
            st.finalResidual = rNorm;
            break;
        }
        rho_old = rho_new;

        if (it == cfg_.maxIterations) {
            st.iterations    = it;
            st.finalResidual = rNorm;
            st.converged     = false;
        }
    }

    lastResidual_ = st.finalResidual;
    dX_.to_host();
    std::memcpy(x.data(), dX_.host_ptr(), n_ * sizeof(double));

    const auto t1 = std::chrono::steady_clock::now();
    st.timeSeconds = std::chrono::duration<double>(t1 - t0).count();
    SIMALL_LOG_INFO("Gpu", "BicgstabGpu: it=", st.iterations,
        " ||r||=", st.finalResidual,
        " converged=", st.converged, " time=", st.timeSeconds, "s");
    return st;
}

}  // namespace simall::gpu
