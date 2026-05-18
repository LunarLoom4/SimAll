// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/MpiContext.cpp
// =============================================================================
#include "parallel/MpiContext.hpp"
#include "core/Logger.hpp"

#if defined(SIMALL_HAVE_MPI)
#  include <mpi.h>
#endif

#include <algorithm>
#include <cstring>
#include <memory>

namespace simall::parallel {

namespace {
#if defined(SIMALL_HAVE_MPI)
inline MPI_Comm as_comm(void* p) {
    return reinterpret_cast<MPI_Comm>(p);
}
#endif
}  // namespace

MpiContext::MpiContext(void* comm, std::string name, bool ownsComm)
    : comm_(comm), name_(std::move(name)), ownsComm_(ownsComm) {
#if defined(SIMALL_HAVE_MPI)
    if (comm_) {
        MPI_Comm_rank(as_comm(comm_), &rank_);
        MPI_Comm_size(as_comm(comm_), &size_);
    }
#endif
}

MpiContext::~MpiContext() {
#if defined(SIMALL_HAVE_MPI)
    if (ownsComm_ && comm_) {
        MPI_Comm c = as_comm(comm_);
        MPI_Comm_free(&c);
    }
#endif
}

MpiContext& MpiContext::world() {
    static MpiContext w(
#if defined(SIMALL_HAVE_MPI)
        reinterpret_cast<void*>(MPI_COMM_WORLD),
#else
        nullptr,
#endif
        "world", false);
    return w;
}

std::unique_ptr<MpiContext> MpiContext::split(int colour, int key) const {
#if defined(SIMALL_HAVE_MPI)
    MPI_Comm newComm = MPI_COMM_NULL;
    MPI_Comm_split(as_comm(comm_), colour, key, &newComm);
    return std::unique_ptr<MpiContext>(new MpiContext(
        reinterpret_cast<void*>(newComm),
        name_ + ".split", true));
#else
    (void)colour; (void)key;
    return std::unique_ptr<MpiContext>(new MpiContext(nullptr, name_ + ".split", false));
#endif
}

// ----------------------------------------------------------------- blocking p2p
void MpiContext::send_doubles(const double* buf, std::size_t n, int dest, int tag) const {
#if defined(SIMALL_HAVE_MPI)
    MPI_Send(buf, static_cast<int>(n), MPI_DOUBLE, dest, tag, as_comm(comm_));
#else
    (void)buf; (void)n; (void)dest; (void)tag;
#endif
}
void MpiContext::recv_doubles(double* buf, std::size_t n, int src, int tag) const {
#if defined(SIMALL_HAVE_MPI)
    MPI_Recv(buf, static_cast<int>(n), MPI_DOUBLE, src, tag,
             as_comm(comm_), MPI_STATUS_IGNORE);
#else
    (void)buf; (void)n; (void)src; (void)tag;
#endif
}
void MpiContext::send_int32s(const std::int32_t* buf, std::size_t n, int dest, int tag) const {
#if defined(SIMALL_HAVE_MPI)
    MPI_Send(buf, static_cast<int>(n), MPI_INT32_T, dest, tag, as_comm(comm_));
#else
    (void)buf; (void)n; (void)dest; (void)tag;
#endif
}
void MpiContext::recv_int32s(std::int32_t* buf, std::size_t n, int src, int tag) const {
#if defined(SIMALL_HAVE_MPI)
    MPI_Recv(buf, static_cast<int>(n), MPI_INT32_T, src, tag,
             as_comm(comm_), MPI_STATUS_IGNORE);
#else
    (void)buf; (void)n; (void)src; (void)tag;
#endif
}

// --------------------------------------------------------------- non-blocking
MpiRequest MpiContext::isend_doubles(const double* buf, std::size_t n, int dest, int tag) const {
    MpiRequest r{};
#if defined(SIMALL_HAVE_MPI)
    MPI_Request req;
    MPI_Isend(buf, static_cast<int>(n), MPI_DOUBLE, dest, tag, as_comm(comm_), &req);
    static_assert(sizeof(std::uintptr_t) >= sizeof(MPI_Request),
                  "MPI_Request must fit in MpiRequest::data");
    std::memcpy(&r.data, &req, sizeof(req));
    r.is_active = true;
#else
    (void)buf; (void)n; (void)dest; (void)tag;
#endif
    return r;
}
MpiRequest MpiContext::irecv_doubles(double* buf, std::size_t n, int src, int tag) const {
    MpiRequest r{};
#if defined(SIMALL_HAVE_MPI)
    MPI_Request req;
    MPI_Irecv(buf, static_cast<int>(n), MPI_DOUBLE, src, tag, as_comm(comm_), &req);
    std::memcpy(&r.data, &req, sizeof(req));
    r.is_active = true;
#else
    (void)buf; (void)n; (void)src; (void)tag;
#endif
    return r;
}

void MpiContext::wait_all(std::vector<MpiRequest>& reqs) const {
#if defined(SIMALL_HAVE_MPI)
    std::vector<MPI_Request> raw; raw.reserve(reqs.size());
    for (auto& r : reqs) if (r.is_active) {
        MPI_Request req; std::memcpy(&req, &r.data, sizeof(req));
        raw.push_back(req);
    }
    if (!raw.empty())
        MPI_Waitall(static_cast<int>(raw.size()), raw.data(), MPI_STATUSES_IGNORE);
    for (auto& r : reqs) r.is_active = false;
#else
    for (auto& r : reqs) r.is_active = false;
#endif
}

// --------------------------------------------------------------- collectives
void MpiContext::barrier() const {
#if defined(SIMALL_HAVE_MPI)
    MPI_Barrier(as_comm(comm_));
#endif
}

double MpiContext::allreduce_sum(double v) const {
#if defined(SIMALL_HAVE_MPI)
    double g = 0; MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_SUM, as_comm(comm_)); return g;
#else
    return v;
#endif
}
double MpiContext::allreduce_max(double v) const {
#if defined(SIMALL_HAVE_MPI)
    double g = 0; MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_MAX, as_comm(comm_)); return g;
#else
    return v;
#endif
}
double MpiContext::allreduce_min(double v) const {
#if defined(SIMALL_HAVE_MPI)
    double g = 0; MPI_Allreduce(&v, &g, 1, MPI_DOUBLE, MPI_MIN, as_comm(comm_)); return g;
#else
    return v;
#endif
}
void MpiContext::allreduce_sum_array(double* v, std::size_t n) const {
#if defined(SIMALL_HAVE_MPI)
    std::vector<double> g(n);
    MPI_Allreduce(v, g.data(), static_cast<int>(n), MPI_DOUBLE, MPI_SUM, as_comm(comm_));
    std::copy(g.begin(), g.end(), v);
#else
    (void)v; (void)n;
#endif
}
void MpiContext::allgather_int(std::int32_t local, std::vector<std::int32_t>& gathered) const {
    gathered.assign(size_, 0);
#if defined(SIMALL_HAVE_MPI)
    MPI_Allgather(&local, 1, MPI_INT32_T, gathered.data(), 1, MPI_INT32_T, as_comm(comm_));
#else
    gathered[0] = local;
#endif
}

}  // namespace simall::parallel
