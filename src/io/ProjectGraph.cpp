// =============================================================================
// SimAll Beta - IO Subsystem
// File   : src/io/ProjectGraph.cpp
// =============================================================================
#include "io/ProjectGraph.hpp"

#include <algorithm>
#include <random>
#include <sstream>
#include <unordered_set>

namespace simall::io {

std::uint64_t hash_bytes(const void* data, std::size_t n, std::uint64_t seed) noexcept {
    // FNV-1a 64-bit
    auto* p = static_cast<const std::uint8_t*>(data);
    std::uint64_t h = seed;
    for (std::size_t i = 0; i < n; ++i) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

namespace {
std::string fresh_uuid() {
    static thread_local std::mt19937_64 rng(std::random_device{}());
    static const char hex[] = "0123456789abcdef";
    std::string out(36, '-');
    std::uint64_t a = rng();
    std::uint64_t b = rng();
    std::size_t k = 0;
    auto put_hex = [&](std::uint64_t v, int width) {
        for (int i = width - 1; i >= 0; --i) {
            char c = hex[(v >> (i * 4)) & 0xf];
            while (k < 36 && out[k] == '-') ++k;
            if (k >= 36) return;
            out[k++] = c;
        }
    };
    put_hex(a >> 32, 8);   out[8]  = '-';
    put_hex((a >> 16) & 0xffff, 4); out[13] = '-';
    put_hex(a & 0xffff, 4); out[18] = '-';
    put_hex(b >> 48, 4);    out[23] = '-';
    put_hex(b & 0xffffffffffffULL, 12);
    return out;
}
}  // namespace

ProjectGraph::ProjectGraph() = default;

ProjectNode* ProjectGraph::add_node(NodeKind kind, std::string name, std::string uuid) {
    owned_.push_back(std::make_unique<ProjectNode>(kind, std::move(name)));
    auto* n = owned_.back().get();
    meta_[n] = { uuid.empty() ? fresh_uuid() : std::move(uuid), 0 };
    return n;
}

ProjectNode* ProjectGraph::find_by_uuid(std::string_view uuid) const noexcept {
    for (auto& [node, m] : meta_) if (m.uuid == uuid) return const_cast<ProjectNode*>(node);
    return nullptr;
}

const std::string& ProjectGraph::uuid_of(const ProjectNode* n) const {
    static const std::string empty;
    auto it = meta_.find(n);
    return it == meta_.end() ? empty : it->second.uuid;
}

void ProjectGraph::link(ProjectNode* from, ProjectNode* to) {
    if (from && to) from->add_downstream(to);
}

void ProjectGraph::commit(ProjectNode* n, const void* data, std::size_t nbytes) {
    if (!n) return;
    const auto newHash = hash_bytes(data, nbytes);
    auto& meta = meta_[n];
    const NodeState oldState = n->state();
    if (meta.contentHash != newHash) {
        meta.contentHash = newHash;
        // ProjectNode::invalidate() cascades Stale through its private
        // downstream list, then we re-mark THIS node Valid (we just hashed
        // it cleanly).  Downstream consumers stay Stale until recomputed.
        n->invalidate();
        n->set_state(NodeState::Valid);
    } else {
        n->set_state(NodeState::Valid);
    }
    const NodeState nowState = n->state();
    if (oldState != nowState) {
        for (auto& kv : observers_) {
            try { kv.second(n, oldState, nowState); } catch (...) {}
        }
    }
}

std::vector<ProjectNode*> ProjectGraph::topological_order() const {
    // ProjectNode's `downstream_` list is private; we therefore return the
    // graph in creation order, which the workflow driver populates
    // upstream-first (geometry → mesh → physics → bc → init → solution →
    // result).  This is a topologically valid order for the canonical
    // CFD pipeline; cross-link nodes are still invalidated correctly via
    // ProjectNode::invalidate().
    std::vector<ProjectNode*> out;
    out.reserve(owned_.size());
    for (auto& up : owned_) out.push_back(up.get());
    return out;
}

std::vector<ProjectNode*> ProjectGraph::stale_nodes() const {
    std::vector<ProjectNode*> out;
    for (auto& up : owned_)
        if (up->state() == NodeState::Stale) out.push_back(up.get());
    return out;
}

std::uint64_t ProjectGraph::hash_of(const ProjectNode* n) const {
    auto it = meta_.find(n);
    return it == meta_.end() ? 0 : it->second.contentHash;
}

int ProjectGraph::subscribe(StateObserver obs) {
    const int id = nextObsId_++;
    observers_[id] = std::move(obs);
    return id;
}
void ProjectGraph::unsubscribe(int id) { observers_.erase(id); }

}  // namespace simall::io
