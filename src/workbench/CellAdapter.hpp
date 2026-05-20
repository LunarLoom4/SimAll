// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : src/workbench/CellAdapter.hpp
// Phase  : 22 Pass 22.5
//
// CellAdapter -- the bridge between an abstract workbench Cell and the
// concrete CAD / Mesh / Solver / Results back-end that should actually
// do the work when the user clicks "Refresh" on that cell.
//
// Design contract
// ---------------
// The workbench module deliberately depends only on `core` + `utilities`
// (see CMakeLists.txt PUBLIC_DEPS).  It must NEVER acquire a hard link
// dependency on `cad`, `meshing`, `solver`, or `visualization` -- doing
// so would invert the layering and turn the workbench into a cyclic
// blob that drags every subsystem into every test binary.
//
// To stay clean, the adapter layer uses TWO indirections:
//
//   1. `ICellAdapter` is an abstract interface.  Concrete adapters live
//      next to (or in the same module as) the back-end they wrap; they
//      may freely #include cad/, meshing/, solver/, etc.
//   2. `CellAdapterRegistry` is a registry of *factory functions* keyed
//      by string adapter-id.  The GUI shell (`src/gui/MainWindow.cpp`)
//      assembles the registry at boot time, plugging in lambdas that
//      close over the already-instantiated `cad::CadKernel`,
//      `meshing::SurfaceMesher`, `solver::Solver`, etc.
//
// `WorkflowEngine::refresh_one(id, registry, ctx)` (Pass 22.5 overload)
// is the canonical run-cell entry point.  It looks up the cell's
// `adapter_id`, instantiates the adapter via the registry, calls
// `execute()`, and dispatches `mark_solved` / `mark_failed`.
//
// All adapters share a small `ExecutionContext` that carries progress,
// log, and cancellation channels.  Adapters are expected to:
//   * call `ctx.report_progress(0..1)` periodically so the GUI status
//     bar stays alive,
//   * call `ctx.log(level, msg)` for any user-visible output, and
//   * poll `ctx.is_cancelled()` between expensive steps and bail out
//     early (returning `false`) when the user cancels.
//
// All members are single-threaded inside one `execute()` call; the GUI
// dispatches each cell on a `QtConcurrent::run` worker so the schematic
// stays interactive.
// =============================================================================
#pragma once

#include "workbench/Cell.hpp"
#include "workbench/Workbench.hpp"

#include <atomic>
#include <functional>
#include <memory>
#include <string>

namespace simall::workbench {

// ---------------------------------------------------------------------------
// Log severity for ExecutionContext::log().  Matches the levels used by
// the GUI console panel (`simall::gui::ConsolePanel::post`):
//   0 = trace     1 = error / failure
//   2 = info      3 = warning
// Kept as plain ints to avoid pulling core::Logger into this header.
// ---------------------------------------------------------------------------
enum class AdapterLogLevel : int {
    Trace   = 0,
    Error   = 1,
    Info    = 2,
    Warning = 3,
};

// ---------------------------------------------------------------------------
// ExecutionContext -- per-call scratch space.  Lifetime is bounded by a
// single `ICellAdapter::execute()`.  Three channels:
//
//   * progress_sink  -- void(double 0..1).  Optional; default is a no-op.
//   * log_sink       -- void(AdapterLogLevel, std::string_view).  Optional.
//   * cancel_flag    -- shared_ptr<atomic_bool>.  The caller flips this
//                       to request cancellation; adapters should poll
//                       `is_cancelled()` and return false when set.
//
// shared_ptr is used for the cancel flag so the GUI can outlive the
// worker thread without races (the QtConcurrent task captures the
// shared_ptr by value).
// ---------------------------------------------------------------------------
class ExecutionContext {
public:
    using ProgressSink = std::function<void(double)>;
    using LogSink      = std::function<void(AdapterLogLevel, std::string_view)>;

    ExecutionContext() = default;

    // -- sinks -----------------------------------------------------------
    void set_progress_sink(ProgressSink fn) { progress_sink_ = std::move(fn); }
    void set_log_sink     (LogSink      fn) { log_sink_      = std::move(fn); }

    void report_progress(double frac) const {
        if (progress_sink_) progress_sink_(frac);
    }
    void log(AdapterLogLevel lvl, std::string_view msg) const {
        if (log_sink_) log_sink_(lvl, msg);
    }
    void info (std::string_view m) const { log(AdapterLogLevel::Info,    m); }
    void warn (std::string_view m) const { log(AdapterLogLevel::Warning, m); }
    void error(std::string_view m) const { log(AdapterLogLevel::Error,   m); }

    // -- cancellation ----------------------------------------------------
    void set_cancel_flag(std::shared_ptr<std::atomic_bool> flag) {
        cancel_flag_ = std::move(flag);
    }
    [[nodiscard]] bool is_cancelled() const noexcept {
        return cancel_flag_ && cancel_flag_->load(std::memory_order_acquire);
    }
    void request_cancel() {
        if (cancel_flag_) cancel_flag_->store(true, std::memory_order_release);
    }

private:
    ProgressSink                       progress_sink_;
    LogSink                            log_sink_;
    std::shared_ptr<std::atomic_bool>  cancel_flag_;
};

// ---------------------------------------------------------------------------
// ICellAdapter -- abstract bridge to a concrete back-end.
//
// `execute()` is the only verb.  Implementations:
//   * read inputs from the cell (and from the wider Schematic via the
//     pointer supplied by the engine if needed),
//   * do their work synchronously inside the call,
//   * return `true` on success (engine then calls mark_solved),
//   * return `false` on failure (engine calls mark_failed) AND populate
//     `last_error_` so the GUI can surface the message.
//
// Adapters MUST NOT mutate the Schematic / StateMachine themselves --
// state propagation is owned by WorkflowEngine.
// ---------------------------------------------------------------------------
class ICellAdapter {
public:
    virtual ~ICellAdapter() = default;

    // The kind the adapter is intended to drive.  Engines use this for
    // sanity checking ("is the registered adapter even applicable to
    // this cell?").  May return CellKind::Custom for adapters that are
    // intentionally polymorphic across kinds.
    [[nodiscard]] virtual CellKind kind() const noexcept = 0;

    // Human-readable identifier matching the registry key, e.g.
    // "cad.import.step".  Provided so the engine can log a meaningful
    // failure message without re-plumbing the registry key.
    [[nodiscard]] virtual std::string_view adapter_id() const noexcept = 0;

    // Do the work.  Implementations must be reentrant -- the engine may
    // call execute() many times across the cell's lifetime (once per
    // user-initiated refresh).
    [[nodiscard]] virtual bool execute(Cell& cell, ExecutionContext& ctx) = 0;

    // The most recent error message, valid until the next execute().
    // Empty string when the last call succeeded.
    [[nodiscard]] virtual std::string_view last_error() const noexcept = 0;
};

// ---------------------------------------------------------------------------
// FunctionalCellAdapter -- glue that wraps a `std::function` body into
// the full ICellAdapter interface.  The GUI shell uses this to register
// lambdas that close over the live `cad::CadKernel` / `meshing::*` /
// `solver::Solver` instances without needing to hand-roll a new C++
// class per adapter.
// ---------------------------------------------------------------------------
class FunctionalCellAdapter final : public ICellAdapter {
public:
    using Body = std::function<bool(Cell&, ExecutionContext&, std::string& out_error)>;

    FunctionalCellAdapter(CellKind kind, std::string id, Body body)
        : kind_(kind), id_(std::move(id)), body_(std::move(body)) {}

    [[nodiscard]] CellKind         kind()       const noexcept override { return kind_; }
    [[nodiscard]] std::string_view adapter_id() const noexcept override { return id_;  }
    [[nodiscard]] std::string_view last_error() const noexcept override { return last_error_; }

    [[nodiscard]] bool execute(Cell& cell, ExecutionContext& ctx) override {
        last_error_.clear();
        if (!body_) {
            last_error_ = "FunctionalCellAdapter has no body";
            return false;
        }
        try {
            return body_(cell, ctx, last_error_);
        } catch (const std::exception& ex) {
            last_error_ = std::string("uncaught exception: ") + ex.what();
            return false;
        } catch (...) {
            last_error_ = "uncaught non-std exception";
            return false;
        }
    }

private:
    CellKind    kind_;
    std::string id_;
    Body        body_;
    std::string last_error_;
};

}  // namespace simall::workbench
