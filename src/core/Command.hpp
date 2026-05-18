// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/Command.hpp
// Phase  : 1.3 (APPLICATION CORE → Undo/Redo System)
//
// Command pattern with bounded history. Every user-facing mutation
// (geometry edit, BC change, solver setting, named-selection assignment)
// MUST be performed through a Command so that undo/redo is universal.
// Spec: "Every operation must support undo/redo." (Master Blueprint §1.3)
// =============================================================================
#pragma once

#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simall::core {

class ICommand {
public:
    virtual ~ICommand()                = default;
    virtual void        execute()      = 0;
    virtual void        undo()         = 0;
    virtual std::string description() const = 0;
};

class CommandHistory {
public:
    explicit CommandHistory(std::size_t limit = 256) : limit_(limit) {}

    void execute(std::unique_ptr<ICommand> cmd) {
        cmd->execute();
        std::lock_guard lk(mtx_);
        redo_.clear();
        undo_.emplace_back(std::move(cmd));
        if (undo_.size() > limit_) undo_.pop_front();
    }

    bool can_undo() const { std::lock_guard lk(mtx_); return !undo_.empty(); }
    bool can_redo() const { std::lock_guard lk(mtx_); return !redo_.empty(); }

    void undo() {
        std::unique_ptr<ICommand> cmd;
        {
            std::lock_guard lk(mtx_);
            if (undo_.empty()) return;
            cmd = std::move(undo_.back());
            undo_.pop_back();
        }
        cmd->undo();
        std::lock_guard lk(mtx_);
        redo_.emplace_back(std::move(cmd));
    }

    void redo() {
        std::unique_ptr<ICommand> cmd;
        {
            std::lock_guard lk(mtx_);
            if (redo_.empty()) return;
            cmd = std::move(redo_.back());
            redo_.pop_back();
        }
        cmd->execute();
        std::lock_guard lk(mtx_);
        undo_.emplace_back(std::move(cmd));
    }

    void clear() {
        std::lock_guard lk(mtx_);
        undo_.clear(); redo_.clear();
    }

private:
    mutable std::mutex                       mtx_;
    std::size_t                              limit_;
    std::deque<std::unique_ptr<ICommand>>    undo_;
    std::deque<std::unique_ptr<ICommand>>    redo_;
};

// ----------------------------------------------------------------------------
// CompositeCommand — atomically executes / undoes a sequence of child
// commands. Used to record GUI macros and bundle multi-step user edits
// (e.g. "Apply BC group" = create BC + assign faces + reset solver state).
//
// Executes children in registration order; undoes in reverse order.
// If any child throws during execute(), already-executed children are
// rolled back so the system remains consistent.
// ----------------------------------------------------------------------------
class CompositeCommand : public ICommand {
public:
    explicit CompositeCommand(std::string desc = "Composite")
        : desc_(std::move(desc)) {}

    void add(std::unique_ptr<ICommand> child) {
        if (!child) throw std::invalid_argument("CompositeCommand: null child");
        children_.push_back(std::move(child));
    }

    /// Build a composite from any number of commands in one expression:
    ///   auto c = CompositeCommand::of("Edit BC", std::move(a), std::move(b));
    template <class... Cmds>
    static std::unique_ptr<CompositeCommand> of(std::string desc, Cmds&&... cs) {
        auto c = std::make_unique<CompositeCommand>(std::move(desc));
        (c->add(std::forward<Cmds>(cs)), ...);
        return c;
    }

    std::size_t size() const noexcept { return children_.size(); }
    bool        empty() const noexcept { return children_.empty(); }

    void execute() override {
        std::size_t i = 0;
        try {
            for (; i < children_.size(); ++i) children_[i]->execute();
        } catch (...) {
            // Roll back already-executed children in reverse order.
            for (std::size_t k = i; k-- > 0;) {
                try { children_[k]->undo(); } catch (...) { /* swallow */ }
            }
            throw;
        }
    }

    void undo() override {
        for (std::size_t k = children_.size(); k-- > 0;) {
            children_[k]->undo();
        }
    }

    std::string description() const override { return desc_; }

private:
    std::string                                  desc_;
    std::vector<std::unique_ptr<ICommand>>       children_;
};

// ----------------------------------------------------------------------------
// LambdaCommand — convenience adapter for one-off scripted commands:
//   history.execute(std::make_unique<LambdaCommand>(
//       "Set viscosity",
//       [&]{ field.set("mu", v_new); },
//       [&]{ field.set("mu", v_old); }));
// ----------------------------------------------------------------------------
class LambdaCommand : public ICommand {
public:
    LambdaCommand(std::string desc,
                  std::function<void()> doFn,
                  std::function<void()> undoFn)
        : desc_(std::move(desc)),
          do_(std::move(doFn)),
          undo_(std::move(undoFn)) {
        if (!do_ || !undo_) throw std::invalid_argument("LambdaCommand: null functor");
    }
    void execute() override { do_(); }
    void undo()    override { undo_(); }
    std::string description() const override { return desc_; }
private:
    std::string desc_;
    std::function<void()> do_;
    std::function<void()> undo_;
};

}  // namespace simall::core
