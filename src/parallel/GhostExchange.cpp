// =============================================================================
// SimAll Beta - Parallel Subsystem
// File   : src/parallel/GhostExchange.cpp
// =============================================================================
#include "parallel/GhostExchange.hpp"

#include "core/Logger.hpp"

#include <algorithm>

namespace simall::parallel
{

GhostExchange::GhostExchange(MpiContext& ctx) : ctx_(ctx)
{
    SIMALL_LOG_INFO("Parallel", "GhostExchange ctor rank=", ctx_.rank(), "/", ctx_.size());
}

void GhostExchange::set_layout(std::size_t nOwned, std::size_t nGhost)
{
    nOwned_ = nOwned;
    nGhost_ = nGhost;
}

void GhostExchange::add_neighbour(NeighbourComm spec)
{
    nbrs_.push_back(std::move(spec));
}

void GhostExchange::clear_neighbours()
{
    nbrs_.clear();
}

void GhostExchange::sync_scalar(util::aligned_vector<double>& v)
{
    if (v.size() < nOwned_ + nGhost_)
        v.resize(nOwned_ + nGhost_, 0.0);
    if (nbrs_.empty())
        return;

    // Pack send buffers, post Irecv first to avoid deadlock under high
    // message volume (rendezvous protocol corner-cases).
    std::vector<std::vector<double>> sendBufs(nbrs_.size());
    std::vector<std::vector<double>> recvBufs(nbrs_.size());
    std::vector<MpiRequest> reqs;
    const int tagBase = static_cast<int>(TagNamespace::GhostExchange);

    for (std::size_t i = 0; i < nbrs_.size(); ++i) {
        recvBufs[i].assign(nbrs_[i].recvCells.size(), 0.0);
        reqs.push_back(
            ctx_.irecv_doubles(recvBufs[i].data(), recvBufs[i].size(), nbrs_[i].rank, tagBase));
    }
    for (std::size_t i = 0; i < nbrs_.size(); ++i) {
        sendBufs[i].reserve(nbrs_[i].sendCells.size());
        for (auto c : nbrs_[i].sendCells)
            sendBufs[i].push_back(v[c]);
        reqs.push_back(
            ctx_.isend_doubles(sendBufs[i].data(), sendBufs[i].size(), nbrs_[i].rank, tagBase));
    }
    ctx_.wait_all(reqs);

    // Unpack.
    for (std::size_t i = 0; i < nbrs_.size(); ++i) {
        const auto& rc = nbrs_[i].recvCells;
        const auto& rb = recvBufs[i];
        for (std::size_t k = 0; k < rc.size(); ++k)
            v[rc[k]] = rb[k];
    }
}

void GhostExchange::sync_vector(util::aligned_vector<double>& cx,
                                util::aligned_vector<double>& cy,
                                util::aligned_vector<double>& cz)
{
    const std::size_t total = nOwned_ + nGhost_;
    if (cx.size() < total)
        cx.resize(total, 0.0);
    if (cy.size() < total)
        cy.resize(total, 0.0);
    if (cz.size() < total)
        cz.resize(total, 0.0);
    if (nbrs_.empty())
        return;

    std::vector<std::vector<double>> sendBufs(nbrs_.size());
    std::vector<std::vector<double>> recvBufs(nbrs_.size());
    std::vector<MpiRequest> reqs;
    const int tagBase = static_cast<int>(TagNamespace::GhostExchange) + 1;

    for (std::size_t i = 0; i < nbrs_.size(); ++i) {
        recvBufs[i].assign(nbrs_[i].recvCells.size() * 3, 0.0);
        reqs.push_back(
            ctx_.irecv_doubles(recvBufs[i].data(), recvBufs[i].size(), nbrs_[i].rank, tagBase));
    }
    for (std::size_t i = 0; i < nbrs_.size(); ++i) {
        sendBufs[i].reserve(nbrs_[i].sendCells.size() * 3);
        for (auto c : nbrs_[i].sendCells) {
            sendBufs[i].push_back(cx[c]);
            sendBufs[i].push_back(cy[c]);
            sendBufs[i].push_back(cz[c]);
        }
        reqs.push_back(
            ctx_.isend_doubles(sendBufs[i].data(), sendBufs[i].size(), nbrs_[i].rank, tagBase));
    }
    ctx_.wait_all(reqs);

    for (std::size_t i = 0; i < nbrs_.size(); ++i) {
        const auto& rc = nbrs_[i].recvCells;
        const auto& rb = recvBufs[i];
        for (std::size_t k = 0; k < rc.size(); ++k) {
            cx[rc[k]] = rb[3 * k + 0];
            cy[rc[k]] = rb[3 * k + 1];
            cz[rc[k]] = rb[3 * k + 2];
        }
    }
}

} // namespace simall::parallel
