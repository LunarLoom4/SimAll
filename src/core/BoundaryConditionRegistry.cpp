// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/BoundaryConditionRegistry.cpp
// Phase  : 22 Pass 2
// =============================================================================
#include "core/BoundaryConditionRegistry.hpp"

#include <algorithm>

namespace simall::core {

BoundaryConditionRegistry::Handle
BoundaryConditionRegistry::add(BoundaryConditionEntry entry) {
    std::lock_guard<std::mutex> lk(mu_);
    const Handle h = next_++;
    entries_.emplace(h, std::move(entry));
    order_.push_back(h);
    return h;
}

bool BoundaryConditionRegistry::remove(Handle h) {
    std::lock_guard<std::mutex> lk(mu_);
    const std::size_t erased = entries_.erase(h);
    if (erased == 0) return false;
    auto it = std::find(order_.begin(), order_.end(), h);
    if (it != order_.end()) order_.erase(it);
    return true;
}

void BoundaryConditionRegistry::clear() noexcept {
    std::lock_guard<std::mutex> lk(mu_);
    entries_.clear();
    order_.clear();
    next_ = 1;
}

const BoundaryConditionEntry*
BoundaryConditionRegistry::find(Handle h) const noexcept {
    std::lock_guard<std::mutex> lk(mu_);
    auto it = entries_.find(h);
    return it == entries_.end() ? nullptr : &it->second;
}

std::vector<BoundaryConditionRegistry::Handle>
BoundaryConditionRegistry::by_zone(BcZoneId z) const {
    std::vector<Handle> out;
    std::lock_guard<std::mutex> lk(mu_);
    out.reserve(order_.size());
    for (auto h : order_) {
        auto it = entries_.find(h);
        if (it != entries_.end() && it->second.zone == z) out.push_back(h);
    }
    return out;
}

std::vector<BoundaryConditionRegistry::Handle>
BoundaryConditionRegistry::by_variable(std::string_view v) const {
    std::vector<Handle> out;
    std::lock_guard<std::mutex> lk(mu_);
    out.reserve(order_.size());
    for (auto h : order_) {
        auto it = entries_.find(h);
        if (it != entries_.end() && it->second.variable == v) out.push_back(h);
    }
    return out;
}

std::vector<BoundaryConditionRegistry::Handle>
BoundaryConditionRegistry::by_subsystem(std::string_view s) const {
    std::vector<Handle> out;
    std::lock_guard<std::mutex> lk(mu_);
    out.reserve(order_.size());
    for (auto h : order_) {
        auto it = entries_.find(h);
        if (it != entries_.end() && it->second.subsystem == s) out.push_back(h);
    }
    return out;
}

BoundaryConditionRegistry::Handle
BoundaryConditionRegistry::first(BcZoneId z, std::string_view variable) const noexcept {
    std::lock_guard<std::mutex> lk(mu_);
    for (auto h : order_) {
        auto it = entries_.find(h);
        if (it != entries_.end() &&
            it->second.zone == z &&
            it->second.variable == variable) {
            return h;
        }
    }
    return kInvalid;
}

std::size_t BoundaryConditionRegistry::size() const noexcept {
    std::lock_guard<std::mutex> lk(mu_);
    return entries_.size();
}

bool BoundaryConditionRegistry::empty() const noexcept {
    std::lock_guard<std::mutex> lk(mu_);
    return entries_.empty();
}

BoundaryConditionRegistry& BoundaryConditionRegistry::instance() {
    static BoundaryConditionRegistry inst;
    return inst;
}

}  // namespace simall::core
