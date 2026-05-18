// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/GpuKernels.hpp
// Phase  : 18.5 — Linear-algebra kernels for GPU Krylov solvers (W13).
//
// Free functions operating on raw device pointers (or host pointers on
// SIMALL_NO_CUDA fall-back).  All sizes are element counts, not bytes.
//
// Each kernel ships in two implementations:
//
//   * SIMALL_HAVE_CUDA defined  →  thin wrapper that launches a cuBLAS /
//                                  cuSPARSE call, or a hand-written kernel
//                                  for the few cases not covered by the
//                                  toolkit (ILU(0) triangular solves use
//                                  cusparseSpSV; Jacobi uses a tiny lambda
//                                  via thrust::transform when present,
//                                  else cuBLAS Dscal+axpy combo).
//
//   * Otherwise                  →  cache-friendly serial reference impl
//                                  in plain C++ over the same pointers.
//                                  Used by tests and by users who deploy on
//                                  CPU-only nodes.
//
// CONTRACT
// --------
//   * Pointers may NOT alias unless explicitly documented.
//   * Result arrays are overwritten (not accumulated) unless the operation
//     is explicitly accumulative (axpy).
//   * Stream argument is opaque (void*); pass StreamScheduler::acquire()
//     output, or nullptr to use the default stream.
//   * Reductions (dot, nrm2) BLOCK on the host until the result is ready.
//
// CSR storage convention matches solver::CSRMatrix:
//   rowPtr[0..n]        — int
//   colIdx[0..nnz)      — int
//   values[0..nnz)      — double
// =============================================================================
#pragma once

#include <cstddef>

namespace simall::gpu {

// ============================================================ vector ops

/// y[i] += alpha * x[i],  i ∈ [0, n)
void axpy(double alpha, const void* x, void* y, std::size_t n,
          void* stream = nullptr);

/// x[i] *= alpha
void scal(double alpha, void* x, std::size_t n, void* stream = nullptr);

/// dst[i] = src[i]
void copy(const void* src, void* dst, std::size_t n, void* stream = nullptr);

/// Returns x · y   (blocks host until done).
double dot(const void* x, const void* y, std::size_t n, void* stream = nullptr);

/// Returns ||x||₂  (blocks host until done).
double nrm2(const void* x, std::size_t n, void* stream = nullptr);

/// y[i] = a*x[i] + b*y[i]  (BLAS-1 fused).
void axpby(double a, const void* x, double b, void* y, std::size_t n,
           void* stream = nullptr);

/// y[i] = x[i] + beta*y[i]      (CG's p ← r + β·p uses this).
void xpby(const void* x, double beta, void* y, std::size_t n,
          void* stream = nullptr);

// ============================================================ sparse ops

/// y = A·x   in CSR format.
///   rowPtr : int[n+1]
///   colIdx : int[nnz]
///   values : double[nnz]
///   x      : double[n]
///   y      : double[n]
void spmv_csr(std::size_t n,
              const void* rowPtr, const void* colIdx, const void* values,
              const void* x, void* y, void* stream = nullptr);

// ============================================================ preconditioners

/// Diagonal (Jacobi) precondition:  z[i] = r[i] / diag[i].
void jacobi_apply(const void* diag, const void* r, void* z,
                  std::size_t n, void* stream = nullptr);

/// In-place ILU(0) factorisation of a CSR matrix.  The CSR pattern of
/// (L+U-I) coincides with A — only the values array is mutated.  Pre-
/// condition: ``A`` is structurally symmetric with non-zero diagonal.
///
/// Returns false if a zero pivot was encountered (skip preconditioning).
bool ilu0_factor_inplace(std::size_t n,
                         const int* rowPtr, const int* colIdx,
                         double* values);

/// Apply M⁻¹ where M = LU was produced by ilu0_factor_inplace.  Performs
/// a forward + backward triangular sweep.  z := L⁻¹·r ; z := U⁻¹·z.
///
/// On CUDA builds this uses cusparseSpSV; on the host fall-back it is a
/// serial sweep — adequate for unit tests, not for production CPU runs
/// (the CPU ILU(0) lives in solver::ILU0Preconditioner).
void ilu0_apply(std::size_t n,
                const int* rowPtr, const int* colIdx, const double* values,
                const double* r, double* z, void* stream = nullptr);

}  // namespace simall::gpu
