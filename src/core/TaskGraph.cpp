// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/TaskGraph.cpp
// =============================================================================
#include "core/TaskGraph.hpp"

#include <stack>

namespace simall::core
{

void TaskGraph::dependOn(TaskHandle child, TaskHandle dep)
{
    if (started_)
        throw std::runtime_error("TaskGraph: dependOn() after run()");
    if (!child.valid() || !dep.valid() || child.graph_ != this || dep.graph_ != this)
        throw std::runtime_error("TaskGraph: dependOn handle from another graph");
    if (child.id() == dep.id())
        throw std::runtime_error("TaskGraph: self-dependency");
    nodes_[dep.id()]->successors.push_back(child.id());
    ++nodes_[child.id()]->predecessorCount;
}

void TaskGraph::detectCycles() const
{
    enum
    {
        kUnseen,
        kOpen,
        kDone
    };
    std::vector<int> color(nodes_.size(), kUnseen);
    std::vector<std::size_t> stack;
    for (std::size_t root = 0; root < nodes_.size(); ++root) {
        if (color[root] != kUnseen)
            continue;
        stack.push_back(root);
        while (!stack.empty()) {
            std::size_t u = stack.back();
            if (color[u] == kUnseen) {
                color[u] = kOpen;
                for (std::size_t v : nodes_[u]->successors) {
                    if (color[v] == kOpen)
                        throw std::runtime_error("TaskGraph: cycle detected");
                    if (color[v] == kUnseen)
                        stack.push_back(v);
                }
            } else {
                color[u] = kDone;
                stack.pop_back();
            }
        }
    }
}

void TaskGraph::schedule(std::size_t id)
{
    pool_->submit([this, id]() {
        auto& node = *nodes_[id];
        node.fn();
        for (std::size_t succ : node.successors) {
            auto& s = *nodes_[succ];
            if (s.remaining.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                schedule(succ);
            }
        }
        if (outstanding_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            complete_.set_value();
        }
    });
}

void TaskGraph::run()
{
    if (started_)
        throw std::runtime_error("TaskGraph: run() twice");
    started_ = true;
    detectCycles();

    outstanding_.store(nodes_.size(), std::memory_order_release);
    if (nodes_.empty()) {
        complete_.set_value();
        return;
    }
    // Initialise remaining counters and seed roots.
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        nodes_[i]->remaining.store(nodes_[i]->predecessorCount, std::memory_order_release);
    }
    for (std::size_t i = 0; i < nodes_.size(); ++i) {
        if (nodes_[i]->predecessorCount == 0)
            schedule(i);
    }
}

void TaskGraph::waitAll()
{
    if (!started_)
        throw std::runtime_error("TaskGraph: waitAll() before run()");
    completeFut_.wait();
    // Propagate first exception, if any.
    for (auto& node : nodes_) {
        auto fut = node->done->get_future();
        if (fut.valid())
            fut.get(); // rethrows on stored exception
    }
}

} // namespace simall::core
