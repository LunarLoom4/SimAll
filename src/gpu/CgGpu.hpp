// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/CgGpu.hpp
// Phase  : 18.6 — GPU Conjugate-Gradient solver (W13).
//
// Right-preconditioned Conjugate Gradient for SPD systems, mirror of
// solver::CG.  Self-contained: takes a CSR matrix as raw arrays (so the
// gpu/ module need not link solver/), drives the iteration loop on the host,
// and uses GpuKernels for every inner BLAS / sparse / preconditioner call.
//
//   A x = b,   A SPD,   M⁻¹ ≈ A
//
//   r₀ = b - A x₀
//   z₀ = M⁻¹ r₀,   p₀ = z₀
//   for k = 0, 1, …
//       q   = A pₖ
//       α   = (rₖ, zₖ) / (pₖ, q)
//       x   ← x + α pₖ
//       r   ← r - α q
//       check ||r|| ≤ tol
//       z   = M⁻¹ r
//       β   = (rₖ₊₁, zₖ₊₁) / (rₖ, zₖ)
//       p   ← z + β p
//   end
//
// PRECONDITIONERS
// ---------------
//   None     — pure CG
//   Jacobi   — diagonal scaling
//   ILU0     — incomplete LU, host-built then resident on device
//
// FALL-BACK
// ---------
// Identical algorithm runs on the host when SIMALL_HAVE_CUDA is undefined —
// the GpuKernels backend transparently selects the serial implementation.
// =============================================================================
#pragma once

#include "gpu/GpuKernels.hpp"
#include "gpu/HostDeviceMirror.hpp"

#include <cstddef>
#include <memory>
#include <vector>

namespace simall::gpu {

enum class GpuPreconditioner { None, Jacobi, ILU0 };

struct GpuSolverConfig {
    GpuPreconditioner preconditioner = GpuPreconditioner::Jacobi;
    double            tolerance      = 1.0e-8;   ///< ||r|| / ||b||
    int               maxIterations  = 1000;
    /// If positive, stop when ||r||/||b|| drops by this many orders.
    /// Combined with tolerance via "first-trigger-wins" semantics.
    double            relativeOrders = 0.0;
};

struct GpuSolverStats {
    int    iterations    = 0;
    double finalResidual = 0.0;     ///< ||r||₂
    double initialResidual = 0.0;
    bool   converged     = false;
    double timeSeconds   = 0.0;
};

/// Read-only view of a CSR matrix that the GPU solver consumes.  The
/// solver makes its own device copies; the caller retains ownership of
/// the host arrays.
struct HostCsrView {
    std::size_t   n      = 0;
    std::size_t   nnz    = 0;
    const int*    rowPtr = nullptr;     ///< [n+1]
    const int*    colIdx = nullptr;     ///< [nnz]
    const double* values = nullptr;     ///< [nnz]
};

class CgGpu {
public:
    explicit CgGpu(GpuSolverConfig cfg = {});
    ~CgGpu();

    /// Re-build the device-resident matrix copy + preconditioner.  Must be
    /// called whenever ``A`` changes.  ``A.n`` is recorded; subsequent
    /// solve() calls require b/x of matching length.
    void set_matrix(const HostCsrView& A);

    /// Solve A x = b.  ``x`` is both initial guess (in) and solution (out).
    GpuSolverStats solve(const std::vector<double>& b,
                         std::vector<double>&       x);

    const GpuSolverConfig& config() const noexcept { return cfg_; }
    void set_config(GpuSolverConfig c) noexcept    { cfg_ = c; }

    /// Last residual norm returned by solve().
    double last_residual() const noexcept { return lastResidual_; }

private:
    GpuSolverConfig cfg_;
    double          lastResidual_ = 0.0;

    // Device-resident matrix.
    HostDeviceMirror<int>    dRowPtr_;
    HostDeviceMirror<int>    dColIdx_;
    HostDeviceMirror<double> dValues_;
    std::size_t              n_   = 0;
    std::size_t              nnz_ = 0;

    // Preconditioner state.
    HostDeviceMirror<double> dDiag_;        ///< Jacobi
    HostDeviceMirror<double> dLU_;          ///< ILU0 factored values
    bool                     iluValid_ = false;

    // Krylov workspace.
    HostDeviceMirror<double> dB_, dX_, dR_, dZ_, dP_, dQ_;
};

}  // namespace simall::gpu
