// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/ChangeJournal.cpp
// Phase  : 22 Pass 22.4
// =============================================================================
#include "workbench/ChangeJournal.hpp"

namespace simall::workbench {

void ChangeJournal::record(ChangeDirection dir, std::string description) {
    ChangeEntry e;
    Sink        sink_copy;
    {
        std::lock_guard lk(mtx_);
        e.sequence    = next_seq_++;
        e.direction   = dir;
        e.description = std::move(description);
        ring_.push_back(e);
        while (ring_.size() > capacity_) ring_.pop_front();
        sink_copy = sink_;  // copy under lock; invoke outside
    }
    if (sink_copy) sink_copy(e);
}

void ChangeJournal::set_sink(Sink sink) noexcept {
    std::lock_guard lk(mtx_);
    sink_ = std::move(sink);
}

std::vector<ChangeEntry> ChangeJournal::snapshot() const {
    std::lock_guard lk(mtx_);
    return std::vector<ChangeEntry>{ring_.begin(), ring_.end()};
}

std::size_t ChangeJournal::size() const noexcept {
    std::lock_guard lk(mtx_);
    return ring_.size();
}

void ChangeJournal::clear() noexcept {
    std::lock_guard lk(mtx_);
    ring_.clear();
}

}  // namespace simall::workbench
