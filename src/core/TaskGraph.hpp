// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/TaskGraph.hpp
// Phase  : 1.3 (APPLICATION CORE → Task DAG)
//
// Directed acyclic graph of tasks executed on a ThreadPool. Each node returns
// a value of type T (or void); dependencies are declared with then()/after().
// The graph is one-shot: build it, call run(), then collect results from the
// returned futures. Reusable graphs can be obtained by `clone()`.
//
// Used by:
//   - meshing pipelines (parse → tessellate → improve → partition)
//   - solver per-step assembly (gradient → convection → diffusion → source)
//   - postprocessing filter graphs.
// =============================================================================
#pragma once

#include "core/ThreadPool.hpp"

#include <atomic>
#include <functional>
#include <future>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace simall::core
{

class TaskGraph;

class TaskHandle
{
public:
    TaskHandle() = default;

    /// Convenience: schedule a continuation that runs after this task.
    template <class F> TaskHandle then(F&& f) const;

    std::size_t id() const noexcept { return id_; }
    bool valid() const noexcept { return graph_ != nullptr; }

private:
    friend class TaskGraph;
    TaskHandle(TaskGraph* g, std::size_t i) : graph_(g), id_(i) {}
    TaskGraph* graph_ = nullptr;
    std::size_t id_ = 0;
};

class TaskGraph
{
public:
    explicit TaskGraph(ThreadPool& pool = ThreadPool::global()) : pool_(&pool) {}
    ~TaskGraph() = default;

    TaskGraph(const TaskGraph&) = delete;
    TaskGraph& operator=(const TaskGraph&) = delete;
    TaskGraph(TaskGraph&&) = default;
    TaskGraph& operator=(TaskGraph&&) = default;

    /// Add a task with no dependencies.
    template <class F> TaskHandle add(F&& f, std::string label = {});

    /// Declare `dep` must complete before `child`.
    void dependOn(TaskHandle child, TaskHandle dep);

    /// Submit the graph; returns once all roots are scheduled. Use waitAll()
    /// for completion. Throws std::runtime_error on cycle detection.
    void run();

    /// Block until every task has completed.
    void waitAll();

    /// Total number of tasks added.
    std::size_t size() const noexcept { return nodes_.size(); }

private:
    struct Node
    {
        std::function<void()> fn;
        std::vector<std::size_t> successors;
        std::atomic<std::size_t> remaining{0}; // unmet predecessors
        std::size_t predecessorCount = 0;
        std::string label;
        std::shared_ptr<std::promise<void>> done = std::make_shared<std::promise<void>>();
    };

    void schedule(std::size_t id);
    void detectCycles() const;

    ThreadPool* pool_;
    std::vector<std::unique_ptr<Node>> nodes_;
    std::atomic<std::size_t> outstanding_{0};
    std::promise<void> complete_;
    std::shared_future<void> completeFut_ = complete_.get_future().share();
    bool started_ = false;
};

template <class F> TaskHandle TaskGraph::add(F&& f, std::string label)
{
    if (started_)
        throw std::runtime_error("TaskGraph: add() after run()");
    auto node = std::make_unique<Node>();
    node->fn = [fn = std::forward<F>(f), p = node->done]() mutable {
        try {
            fn();
            p->set_value();
        } catch (...) {
            p->set_exception(std::current_exception());
        }
    };
    node->label = std::move(label);
    std::size_t id = nodes_.size();
    nodes_.emplace_back(std::move(node));
    return TaskHandle{this, id};
}

template <class F> TaskHandle TaskHandle::then(F&& f) const
{
    if (!graph_)
        throw std::runtime_error("TaskHandle::then on invalid handle");
    TaskHandle next = graph_->add(std::forward<F>(f));
    graph_->dependOn(next, *this);
    return next;
}

} // namespace simall::core
