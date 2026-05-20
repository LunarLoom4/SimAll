// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/BicgstabGpu.hpp
// Phase  : 18.6 — GPU BiCGSTAB solver (W13).
//
// Bi-Conjugate Gradient Stabilised — handles non-symmetric systems that arise
// from momentum balances with strong convection.  Right-preconditioned form
// (van der Vorst 1992):
//
//   r₀  = b - A x₀
//   r̂₀  = r₀                  (shadow residual; held fixed)
//   ρ₀ = α = ω = 1
//   v₀  = p₀ = 0
//
//   for k = 1, 2, …
//       ρₖ   = (r̂₀, rₖ₋₁)
//       β    = (ρₖ / ρₖ₋₁) (α / ω)
//       p    = r + β (p − ω v)
//       p̂    = M⁻¹ p
//       v    = A p̂
//       α    = ρₖ / (r̂₀, v)
//       s    = r − α v
//       check ||s|| ≤ tol           ;  if so, x ← x + α p̂ and return
//       ŝ    = M⁻¹ s
//       t    = A ŝ
//       ω    = (t, s) / (t, t)
//       x    ← x + α p̂ + ω ŝ
//       r    ← s − ω t
//       check ||r|| ≤ tol
//   end
//
// Same kernel + preconditioner palette as CgGpu.
// =============================================================================
#pragma once

#include "gpu/CgGpu.hpp" // re-uses GpuSolverConfig / GpuSolverStats / HostCsrView
#include "gpu/HostDeviceMirror.hpp"

namespace simall::gpu
{

class BicgstabGpu
{
public:
    explicit BicgstabGpu(GpuSolverConfig cfg = {});
    ~BicgstabGpu();

    void set_matrix(const HostCsrView& A);
    GpuSolverStats solve(const std::vector<double>& b, std::vector<double>& x);

    const GpuSolverConfig& config() const noexcept { return cfg_; }
    void set_config(GpuSolverConfig c) noexcept { cfg_ = c; }
    double last_residual() const noexcept { return lastResidual_; }

private:
    GpuSolverConfig cfg_;
    double lastResidual_ = 0.0;

    HostDeviceMirror<int> dRowPtr_;
    HostDeviceMirror<int> dColIdx_;
    HostDeviceMirror<double> dValues_;
    std::size_t n_ = 0;
    std::size_t nnz_ = 0;

    HostDeviceMirror<double> dDiag_;
    HostDeviceMirror<double> dLU_;
    bool iluValid_ = false;

    // Workspace (Krylov vectors named after van der Vorst's notation).
    HostDeviceMirror<double> dB_, dX_, dR_, dRhat_, dP_, dPhat_, dV_, dS_, dShat_, dT_;
};

} // namespace simall::gpu
