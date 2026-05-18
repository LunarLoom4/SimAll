// =============================================================================
// SimAll Beta - GPU Subsystem
// File   : src/gpu/GpuKernels.cpp
// Phase  : 18.5 — Linear-algebra kernels for GPU Krylov solvers (W13).
// =============================================================================
#include "gpu/GpuKernels.hpp"
#include "gpu/CudaContext.hpp"
#include "gpu/StreamScheduler.hpp"
#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#ifdef SIMALL_HAVE_CUDA
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    #include <cusparse.h>
#endif

namespace simall::gpu {

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
namespace {

inline bool use_cuda() {
#ifdef SIMALL_HAVE_CUDA
    return is_cuda_available();
#else
    return false;
#endif
}

#ifdef SIMALL_HAVE_CUDA
inline cudaStream_t to_stream(void* s) {
    return s ? static_cast<cudaStream_t>(s) : cudaStream_t(0);
}
inline cublasHandle_t   blas()   { return static_cast<cublasHandle_t>  (CudaContext::current().cublas_handle());   }
inline cusparseHandle_t sparse() { return static_cast<cusparseHandle_t>(CudaContext::current().cusparse_handle()); }
#endif

// Host fall-back views (interpret void* as double* / int*).
inline       double* hd(void* p)         { return static_cast<double*>(p);       }
inline const double* hd(const void* p)   { return static_cast<const double*>(p); }
inline const int*    hi(const void* p)   { return static_cast<const int*>(p);    }

}  // namespace

// =============================================================================
// vector ops
// =============================================================================

void axpy(double alpha, const void* x, void* y, std::size_t n, void* stream) {
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        auto h = blas();
        cublasSetStream(h, to_stream(stream));
        cublasDaxpy(h, static_cast<int>(n), &alpha, hd(x), 1, hd(y), 1);
        return;
    }
#else
    (void)stream;
#endif
    const double* xp = hd(x); double* yp = hd(y);
    for (std::size_t i = 0; i < n; ++i) yp[i] += alpha * xp[i];
}

void scal(double alpha, void* x, std::size_t n, void* stream) {
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        auto h = blas();
        cublasSetStream(h, to_stream(stream));
        cublasDscal(h, static_cast<int>(n), &alpha, hd(x), 1);
        return;
    }
#else
    (void)stream;
#endif
    double* xp = hd(x);
    for (std::size_t i = 0; i < n; ++i) xp[i] *= alpha;
}

void copy(const void* src, void* dst, std::size_t n, void* stream) {
    const std::size_t bytes = n * sizeof(double);
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        cudaMemcpyAsync(dst, src, bytes, cudaMemcpyDeviceToDevice, to_stream(stream));
        return;
    }
#else
    (void)stream;
#endif
    std::memcpy(dst, src, bytes);
}

double dot(const void* x, const void* y, std::size_t n, void* stream) {
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        auto h = blas();
        cublasSetStream(h, to_stream(stream));
        double r = 0.0;
        cublasDdot(h, static_cast<int>(n), hd(x), 1, hd(y), 1, &r);
        // cublasDdot with host pointer mode blocks until result is ready.
        return r;
    }
#else
    (void)stream;
#endif
    const double* xp = hd(x); const double* yp = hd(y);
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += xp[i] * yp[i];
    return s;
}

double nrm2(const void* x, std::size_t n, void* stream) {
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        auto h = blas();
        cublasSetStream(h, to_stream(stream));
        double r = 0.0;
        cublasDnrm2(h, static_cast<int>(n), hd(x), 1, &r);
        return r;
    }
#else
    (void)stream;
#endif
    const double* xp = hd(x);
    double s = 0.0;
    for (std::size_t i = 0; i < n; ++i) s += xp[i] * xp[i];
    return std::sqrt(s);
}

void axpby(double a, const void* x, double b, void* y, std::size_t n, void* stream) {
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        auto h = blas();
        cublasSetStream(h, to_stream(stream));
        cublasDscal(h, static_cast<int>(n), &b, hd(y), 1);
        cublasDaxpy(h, static_cast<int>(n), &a, hd(x), 1, hd(y), 1);
        return;
    }
#else
    (void)stream;
#endif
    const double* xp = hd(x); double* yp = hd(y);
    for (std::size_t i = 0; i < n; ++i) yp[i] = a * xp[i] + b * yp[i];
}

void xpby(const void* x, double beta, void* y, std::size_t n, void* stream) {
    // y ← x + β·y
    axpby(1.0, x, beta, y, n, stream);
}

// =============================================================================
// sparse ops
// =============================================================================

void spmv_csr(std::size_t n,
              const void* rowPtr, const void* colIdx, const void* values,
              const void* x, void* y, void* stream)
{
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        auto h = sparse();
        cusparseSetStream(h, to_stream(stream));

        // Build descriptors on the fly.  Allocation is cheap (microseconds).
        // For inner-loop SpMV the calling solver should cache a long-lived
        // cusparseSpMatDescr_t; CgGpu / BicgstabGpu do exactly that, so this
        // generic path is the fall-back for one-off launches.
        cusparseSpMatDescr_t mat = nullptr;
        cusparseDnVecDescr_t vx  = nullptr, vy = nullptr;
        const int nnz = static_cast<int>(static_cast<const int*>(rowPtr)[n]);

        cusparseCreateCsr(&mat, n, n, nnz,
            const_cast<void*>(rowPtr), const_cast<void*>(colIdx),
            const_cast<void*>(values),
            CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
            CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F);
        cusparseCreateDnVec(&vx, n, const_cast<void*>(x), CUDA_R_64F);
        cusparseCreateDnVec(&vy, n, y,                    CUDA_R_64F);

        const double alpha = 1.0, beta = 0.0;
        std::size_t bufBytes = 0;
        cusparseSpMV_bufferSize(h, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, mat, vx, &beta, vy, CUDA_R_64F,
            CUSPARSE_SPMV_ALG_DEFAULT, &bufBytes);
        void* buf = nullptr;
        if (bufBytes > 0) cudaMalloc(&buf, bufBytes);

        cusparseSpMV(h, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &alpha, mat, vx, &beta, vy, CUDA_R_64F,
            CUSPARSE_SPMV_ALG_DEFAULT, buf);

        if (buf) cudaFree(buf);
        cusparseDestroyDnVec(vx);
        cusparseDestroyDnVec(vy);
        cusparseDestroySpMat(mat);
        return;
    }
#else
    (void)stream;
#endif
    // Host reference.
    const int*    rp = hi(rowPtr);
    const int*    ci = hi(colIdx);
    const double* vs = hd(values);
    const double* xp = hd(x);
    double*       yp = hd(y);
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (int k = rp[i]; k < rp[i + 1]; ++k)
            s += vs[k] * xp[ci[k]];
        yp[i] = s;
    }
}

// =============================================================================
// preconditioners
// =============================================================================

void jacobi_apply(const void* diag, const void* r, void* z,
                  std::size_t n, void* stream)
{
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        // No native cuBLAS for element-wise division — copy r→z then divide
        // via two cuSPARSE call or a tiny dedicated kernel.  Because the
        // host build keeps this routine portable, we provide a minimal
        // device kernel inline via a fused launch.  For SimAll's typical
        // n=1e6 this single kernel is bandwidth-bound and adequate; the
        // CG/BiCGSTAB drivers also cache an AMG-class preconditioner for
        // larger systems.
        //
        // Implementation strategy: cudaMemcpyAsync r→z (DtoD), then
        // launch a 1D element-wise divide.  We embed the kernel here to
        // avoid splitting the translation unit into a .cu file.
        //
        // NOTE: this requires nvcc / CUDA runtime compilation.  When the
        // build system links GpuKernels.cpp under a host compiler (which
        // is the SimAll convention), we cannot launch raw kernels.  We
        // therefore implement the apply via two cuBLAS calls:
        //   1) cublasDcopy(r, z)
        //   2) build temporary host vector of 1/diag, upload, multiply
        // step (2) costs an HtoD per call — acceptable for unit tests
        // and small problems; production code uses AMG / multi-colour
        // ILU written in a .cu file (W18.5 milestone).
        auto h = blas();
        cublasSetStream(h, to_stream(stream));
        cudaMemcpyAsync(z, r, n * sizeof(double),
            cudaMemcpyDeviceToDevice, to_stream(stream));
        // Per-element multiply: not natively supported by cuBLAS; the
        // serial host fall-back below copies diag down, computes
        // z[i] /= diag[i], and pushes z back.  This is intentionally a
        // correctness-first implementation; performance preconditioners
        // live in a follow-on .cu translation unit.
        std::vector<double> hostDiag(n), hostZ(n);
        cudaMemcpy(hostDiag.data(), diag, n * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(hostZ.data(),    z,    n * sizeof(double), cudaMemcpyDeviceToHost);
        for (std::size_t i = 0; i < n; ++i)
            hostZ[i] /= (hostDiag[i] != 0.0 ? hostDiag[i] : 1.0);
        cudaMemcpy(z, hostZ.data(), n * sizeof(double), cudaMemcpyHostToDevice);
        return;
    }
#else
    (void)stream;
#endif
    const double* dp = hd(diag);
    const double* rp = hd(r);
    double*       zp = hd(z);
    for (std::size_t i = 0; i < n; ++i)
        zp[i] = rp[i] / (dp[i] != 0.0 ? dp[i] : 1.0);
}

// -----------------------------------------------------------------------------
// ILU(0) factorisation — runs on the host even in CUDA builds.  The factored
// matrix is stored in CSR with the same sparsity pattern as A.  We then
// upload values back to the device for the apply step.
//
// Reference: Saad, "Iterative Methods for Sparse Linear Systems", 2nd ed.,
// Algorithm 10.4 ("ILU(0) factorization").
// -----------------------------------------------------------------------------
bool ilu0_factor_inplace(std::size_t n,
                         const int* rowPtr, const int* colIdx, double* values)
{
    // diagPos[i] = position in values[] of A(i,i).  -1 means missing.
    std::vector<int> diagPos(n, -1);
    for (std::size_t i = 0; i < n; ++i) {
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
            if (colIdx[k] == static_cast<int>(i)) { diagPos[i] = k; break; }
        }
        if (diagPos[i] < 0) return false;  // no diagonal entry
    }

    for (std::size_t i = 1; i < n; ++i) {
        // For each k in row i with k < i, scale by U(k,k) and eliminate.
        for (int p = rowPtr[i]; p < rowPtr[i + 1]; ++p) {
            const int k = colIdx[p];
            if (k >= static_cast<int>(i)) break;        // (rowPtr assumed sorted)
            const double Ukk = values[diagPos[k]];
            if (Ukk == 0.0) return false;
            const double lik = values[p] / Ukk;
            values[p] = lik;                            // L(i,k) ← lik
            // Subtract lik * U(k, j) from A(i, j) for j ∈ row k, j > k that
            // also appears in row i (zero-fill preserves pattern of A).
            for (int q = diagPos[k] + 1; q < rowPtr[k + 1]; ++q) {
                const int j = colIdx[q];
                // Linear search inside row i for column j (rows are short).
                for (int r = p + 1; r < rowPtr[i + 1]; ++r) {
                    if (colIdx[r] == j) { values[r] -= lik * values[q]; break; }
                }
            }
        }
        if (values[diagPos[i]] == 0.0) return false;
    }
    return true;
}

// -----------------------------------------------------------------------------
// Apply M^-1 r → z   given L+U-I stored in CSR(values, rowPtr, colIdx).
//
// Forward solve:  L z = r              ( unit diagonal )
//   z[i] = r[i] - Σ_{j<i, A(i,j)≠0} L(i,j) · z[j]
// Backward solve: U z' = z
//   z'[i] = ( z[i] - Σ_{j>i, A(i,j)≠0} U(i,j) · z'[j] ) / U(i,i)
//
// On CUDA builds we currently pull values to host, sweep, push z back — the
// dependency chain of triangular solves needs cusparseSpSV with analysis
// state, which CgGpu/BicgstabGpu set up directly for performance.  This
// generic helper is the correctness reference exercised by the unit tests.
// -----------------------------------------------------------------------------
void ilu0_apply(std::size_t n,
                const int* rowPtr, const int* colIdx, const double* values,
                const double* r, double* z, void* stream)
{
    (void)stream;
#ifdef SIMALL_HAVE_CUDA
    if (use_cuda()) {
        // Host sweep with explicit copies.  See header note.
        std::vector<int>    hRow(n + 1);
        std::vector<int>    hCol;
        std::vector<double> hVal, hR(n), hZ(n);
        cudaMemcpy(hRow.data(), rowPtr, (n + 1) * sizeof(int),    cudaMemcpyDeviceToHost);
        const int nnz = hRow[n];
        hCol.resize(nnz); hVal.resize(nnz);
        cudaMemcpy(hCol.data(), colIdx, nnz * sizeof(int),    cudaMemcpyDeviceToHost);
        cudaMemcpy(hVal.data(), values, nnz * sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(hR.data(),   r,      n   * sizeof(double), cudaMemcpyDeviceToHost);

        // Forward
        for (std::size_t i = 0; i < n; ++i) {
            double s = hR[i];
            for (int k = hRow[i]; k < hRow[i + 1]; ++k) {
                const int j = hCol[k];
                if (j < static_cast<int>(i)) s -= hVal[k] * hZ[j];
                else break;  // remainder belongs to U
            }
            hZ[i] = s;
        }
        // Backward
        for (std::size_t ii = n; ii-- > 0; ) {
            double s = hZ[ii];
            double diag = 1.0;
            for (int k = hRow[ii]; k < hRow[ii + 1]; ++k) {
                const int j = hCol[k];
                if (j > static_cast<int>(ii)) s    -= hVal[k] * hZ[j];
                else if (j == static_cast<int>(ii)) diag = hVal[k];
            }
            hZ[ii] = s / (diag != 0.0 ? diag : 1.0);
        }
        cudaMemcpy(z, hZ.data(), n * sizeof(double), cudaMemcpyHostToDevice);
        return;
    }
#endif
    // Pure host path.
    for (std::size_t i = 0; i < n; ++i) {
        double s = r[i];
        for (int k = rowPtr[i]; k < rowPtr[i + 1]; ++k) {
            const int j = colIdx[k];
            if (j < static_cast<int>(i)) s -= values[k] * z[j];
            else break;
        }
        z[i] = s;
    }
    for (std::size_t ii = n; ii-- > 0; ) {
        double s = z[ii];
        double diag = 1.0;
        for (int k = rowPtr[ii]; k < rowPtr[ii + 1]; ++k) {
            const int j = colIdx[k];
            if (j > static_cast<int>(ii)) s    -= values[k] * z[j];
            else if (j == static_cast<int>(ii)) diag = values[k];
        }
        z[ii] = s / (diag != 0.0 ? diag : 1.0);
    }
}

}  // namespace simall::gpu
