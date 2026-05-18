// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/MpiContext.hpp
// Phase  : 17 — first-class MPI communicator wrapper.
//
// The legacy free-function façade in parallel/Parallel.hpp uses an
// implicit MPI_COMM_WORLD.  Production CFD codes need:
//   * named sub-communicators (split by colour, e.g. one per overset
//     mesh group, or one per multigrid level),
//   * typed point-to-point operations with async request handles,
//   * deterministic tag namespaces so that ghost-exchange tags do not
//     collide with reduction / file-IO tags.
//
// MpiContext owns an MPI_Comm (opaque handle wrapped as `void*` so this
// header does NOT include <mpi.h>) and provides the typed primitives
// consumed by GhostExchange, FieldReducer and LoadBalancer.
// =============================================================================
#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace simall::parallel {

/// Tag-namespace allocator.  Modules pre-allocate a contiguous range so
/// that overlapping non-blocking exchanges never collide on the same tag.
enum class TagNamespace : std::uint16_t {
    GhostExchange   = 1000,
    FieldReducer    = 2000,
    LoadBalancer    = 3000,
    DomainPartition = 4000,
    Adjoint         = 5000,
    User            = 9000
};

/// Opaque request handle.  In serial builds these are dummies; in MPI
/// builds the underlying MPI_Request is stored in `data`.
struct MpiRequest {
    std::uintptr_t data = 0;
    bool           is_active = false;
};

class MpiContext {
public:
    /// Build a context wrapping MPI_COMM_WORLD (or the serial single-rank
    /// fall-back when SIMALL_HAVE_MPI is not defined).
    static MpiContext& world();

    /// Split the current context's communicator by colour/key (analogous
    /// to MPI_Comm_split).  In serial mode returns a copy of `*this`.
    std::unique_ptr<MpiContext> split(int colour, int key) const;

    int  rank() const noexcept { return rank_; }
    int  size() const noexcept { return size_; }
    bool is_master() const noexcept { return rank_ == 0; }
    const std::string& name() const noexcept { return name_; }

    /// Opaque communicator handle (caller must `reinterpret_cast` to
    /// MPI_Comm in modules that include <mpi.h>).
    void* native_handle() const noexcept { return comm_; }

    // ------------- typed point-to-point (blocking) -------------
    void send_doubles  (const double* buf, std::size_t n, int dest, int tag) const;
    void recv_doubles  (double* buf,       std::size_t n, int src,  int tag) const;
    void send_int32s   (const std::int32_t* buf, std::size_t n, int dest, int tag) const;
    void recv_int32s   (std::int32_t* buf,       std::size_t n, int src,  int tag) const;

    // ------------- typed point-to-point (non-blocking) ---------
    MpiRequest isend_doubles(const double* buf, std::size_t n, int dest, int tag) const;
    MpiRequest irecv_doubles(double* buf,       std::size_t n, int src,  int tag) const;

    void wait_all(std::vector<MpiRequest>& reqs) const;

    // ------------- collectives ----------------------------------
    void barrier() const;
    double allreduce_sum (double v) const;
    double allreduce_max (double v) const;
    double allreduce_min (double v) const;
    void   allreduce_sum_array(double* v, std::size_t n) const;
    void   allgather_int(std::int32_t local, std::vector<std::int32_t>& gathered) const;

    ~MpiContext();
    MpiContext(MpiContext&&)            = default;
    MpiContext& operator=(MpiContext&&) = default;
    MpiContext(const MpiContext&)            = delete;
    MpiContext& operator=(const MpiContext&) = delete;

private:
    MpiContext(void* comm, std::string name, bool ownsComm);

    void*       comm_     = nullptr;   // MPI_Comm (or nullptr for serial)
    std::string name_;
    int         rank_     = 0;
    int         size_     = 1;
    bool        ownsComm_ = false;
};

}  // namespace simall::parallel
