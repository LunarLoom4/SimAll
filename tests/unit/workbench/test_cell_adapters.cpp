// =============================================================================
// SimAll Beta — Workbench Subsystem
// File   : tests/unit/workbench/test_cell_adapters.cpp
// Phase  : 22 Pass 22.5
//
// Covers the CellAdapter / CellAdapterRegistry layer and the new
// registry-aware WorkflowEngine::refresh_one() overload that turns an
// in-graph cell into a real back-end invocation.
// =============================================================================
#include <catch2/catch_test_macros.hpp>

#include "workbench/Cell.hpp"
#include "workbench/CellAdapter.hpp"
#include "workbench/CellAdapterRegistry.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/WorkflowEngine.hpp"

#include <atomic>
#include <memory>
#include <string>
#include <vector>

using namespace simall::workbench;

namespace {

// Minimal hand-rolled adapter that records every call so the tests can
// inspect exactly what the engine handed it.  Uses a shared_ptr to a
// recorder so multiple factory instantiations can converge on the same
// log -- mimicking how a real adapter might capture a back-end singleton.
struct Recorder {
    int  calls = 0;
    bool last_cancelled = false;
    Cell* last_cell     = nullptr;
};

class TestAdapter final : public ICellAdapter {
public:
    TestAdapter(std::shared_ptr<Recorder> rec, bool succeed, std::string err)
        : rec_(std::move(rec)), succeed_(succeed), err_(std::move(err)) {}

    CellKind         kind()       const noexcept override { return CellKind::Custom; }
    std::string_view adapter_id() const noexcept override { return "test.adapter"; }
    std::string_view last_error() const noexcept override { return err_;          }

    bool execute(Cell& cell, ExecutionContext& ctx) override {
        rec_->calls++;
        rec_->last_cell      = &cell;
        rec_->last_cancelled = ctx.is_cancelled();
        ctx.report_progress(0.5);
        ctx.info("test adapter running");
        if (ctx.is_cancelled()) {
            err_ = "cancelled";
            return false;
        }
        return succeed_;
    }

private:
    std::shared_ptr<Recorder> rec_;
    bool                      succeed_;
    std::string               err_;
};

}  // namespace


TEST_CASE("CellAdapterRegistry: register/has/make/unregister/keys",
          "[workbench][adapter][registry]") {
    CellAdapterRegistry reg;
    REQUIRE(reg.empty());

    auto factory = [] {
        auto rec = std::make_shared<Recorder>();
        return std::make_unique<TestAdapter>(rec, true, "");
    };

    SECTION("register new id returns true") {
        REQUIRE(reg.register_factory("alpha", factory));
        REQUIRE(reg.has("alpha"));
        REQUIRE(reg.size() == 1);
    }
    SECTION("re-register returns false but overwrites") {
        reg.register_factory("alpha", factory);
        REQUIRE_FALSE(reg.register_factory("alpha", factory));
    }
    SECTION("make() returns null for unknown id") {
        REQUIRE(reg.make("nope") == nullptr);
    }
    SECTION("make() returns a fresh adapter") {
        reg.register_factory("alpha", factory);
        auto a = reg.make("alpha");
        REQUIRE(a != nullptr);
        REQUIRE(a->adapter_id() == "test.adapter");
    }
    SECTION("unregister is idempotent + reports removal") {
        reg.register_factory("alpha", factory);
        REQUIRE(reg.unregister("alpha"));
        REQUIRE_FALSE(reg.unregister("alpha"));
    }
    SECTION("keys() are sorted") {
        reg.register_factory("z", factory);
        reg.register_factory("a", factory);
        reg.register_factory("m", factory);
        const auto k = reg.keys();
        REQUIRE(k.size() == 3);
        REQUIRE(k[0] == "a");
        REQUIRE(k[1] == "m");
        REQUIRE(k[2] == "z");
    }
}


TEST_CASE("FunctionalCellAdapter: body invoked, exceptions captured",
          "[workbench][adapter][functional]") {
    SECTION("body returning true succeeds") {
        FunctionalCellAdapter a{
            CellKind::Geometry, "g.test",
            [](Cell&, ExecutionContext& ctx, std::string&) {
                ctx.info("ok");
                return true;
            }};
        Schematic s;
        s.add_cell(CellKind::Geometry, "G");
        ExecutionContext ctx;
        REQUIRE(a.execute(*s.cell(1), ctx));
        REQUIRE(a.last_error().empty());
    }
    SECTION("body returning false propagates error string") {
        FunctionalCellAdapter a{
            CellKind::Mesh, "m.test",
            [](Cell&, ExecutionContext&, std::string& err) {
                err = "bad mesh";
                return false;
            }};
        Schematic s;
        s.add_cell(CellKind::Mesh, "M");
        ExecutionContext ctx;
        REQUIRE_FALSE(a.execute(*s.cell(1), ctx));
        REQUIRE(a.last_error() == "bad mesh");
    }
    SECTION("thrown std::exception is caught and reported") {
        FunctionalCellAdapter a{
            CellKind::Custom, "boom",
            [](Cell&, ExecutionContext&, std::string&) -> bool {
                throw std::runtime_error("kaboom");
            }};
        Schematic s;
        s.add_cell(CellKind::Custom, "X");
        ExecutionContext ctx;
        REQUIRE_FALSE(a.execute(*s.cell(1), ctx));
        REQUIRE(a.last_error().find("kaboom") != std::string_view::npos);
    }
    SECTION("missing body fails safely") {
        FunctionalCellAdapter a{CellKind::Custom, "empty", {}};
        Schematic s;
        s.add_cell(CellKind::Custom, "X");
        ExecutionContext ctx;
        REQUIRE_FALSE(a.execute(*s.cell(1), ctx));
        REQUIRE(a.last_error().find("no body") != std::string_view::npos);
    }
}


TEST_CASE("ExecutionContext: cancel flag + sinks",
          "[workbench][adapter][context]") {
    ExecutionContext ctx;
    auto flag = std::make_shared<std::atomic_bool>(false);
    ctx.set_cancel_flag(flag);
    REQUIRE_FALSE(ctx.is_cancelled());
    flag->store(true);
    REQUIRE(ctx.is_cancelled());

    std::vector<double> prog;
    ctx.set_progress_sink([&](double f){ prog.push_back(f); });
    ctx.report_progress(0.25);
    ctx.report_progress(0.75);
    REQUIRE(prog.size() == 2);
    REQUIRE(prog[1] == 0.75);

    std::vector<std::pair<AdapterLogLevel, std::string>> logs;
    ctx.set_log_sink([&](AdapterLogLevel l, std::string_view m){
        logs.emplace_back(l, std::string(m));
    });
    ctx.info("hello");
    ctx.warn("careful");
    ctx.error("boom");
    REQUIRE(logs.size() == 3);
    REQUIRE(logs[0].first  == AdapterLogLevel::Info);
    REQUIRE(logs[2].second == "boom");
}


TEST_CASE("WorkflowEngine::refresh_one(registry) success + state propagation",
          "[workbench][adapter][engine]") {
    Schematic s;
    StateMachine sm{s};
    WorkflowEngine eng{s, sm};
    CellAdapterRegistry reg;

    const CellId g = s.add_cell(CellKind::Geometry, "G");
    s.cell(g)->add_output("brep", "shape");
    s.cell(g)->set_adapter_id("cad.import.fake");

    auto rec = std::make_shared<Recorder>();
    reg.register_factory("cad.import.fake",
        [rec]{ return std::make_unique<TestAdapter>(rec, true, ""); });

    sm.recompute_all();
    REQUIRE(s.cell(g)->state() == CellState::UpToDate);  // no required inputs
    sm.mark_modified(g);
    REQUIRE(s.cell(g)->state() == CellState::RefreshRequired);

    ExecutionContext ctx;
    REQUIRE(eng.refresh_one(g, reg, ctx));
    REQUIRE(rec->calls == 1);
    REQUIRE(s.cell(g)->state() == CellState::UpToDate);
}


TEST_CASE("WorkflowEngine::refresh_one(registry) failure paths mark_failed",
          "[workbench][adapter][engine]") {
    Schematic s;
    StateMachine sm{s};
    WorkflowEngine eng{s, sm};
    CellAdapterRegistry reg;

    const CellId c = s.add_cell(CellKind::Custom, "X");

    SECTION("empty adapter_id -> mark_failed + error logged") {
        std::vector<std::string> errs;
        ExecutionContext ctx;
        ctx.set_log_sink([&](AdapterLogLevel l, std::string_view m){
            if (l == AdapterLogLevel::Error) errs.emplace_back(m);
        });
        REQUIRE_FALSE(eng.refresh_one(c, reg, ctx));
        REQUIRE(s.cell(c)->state() == CellState::Failed);
        REQUIRE_FALSE(errs.empty());
    }
    SECTION("unknown adapter_id -> mark_failed") {
        s.cell(c)->set_adapter_id("never.registered");
        ExecutionContext ctx;
        REQUIRE_FALSE(eng.refresh_one(c, reg, ctx));
        REQUIRE(s.cell(c)->state() == CellState::Failed);
    }
    SECTION("adapter returns false -> mark_failed + last_error propagated") {
        s.cell(c)->set_adapter_id("fail");
        auto rec = std::make_shared<Recorder>();
        reg.register_factory("fail",
            [rec]{ return std::make_unique<TestAdapter>(rec, false, "boom"); });
        std::vector<std::string> errs;
        ExecutionContext ctx;
        ctx.set_log_sink([&](AdapterLogLevel l, std::string_view m){
            if (l == AdapterLogLevel::Error) errs.emplace_back(m);
        });
        REQUIRE_FALSE(eng.refresh_one(c, reg, ctx));
        REQUIRE(s.cell(c)->state() == CellState::Failed);
        REQUIRE_FALSE(errs.empty());
        REQUIRE(errs.front() == "boom");
    }
    SECTION("unknown cell id -> false, no state change") {
        ExecutionContext ctx;
        REQUIRE_FALSE(eng.refresh_one(static_cast<CellId>(999), reg, ctx));
    }
}


TEST_CASE("Cell adapter_id round-trips through Set/SetCellAdapter command",
          "[workbench][adapter][command]") {
    Schematic s;
    const CellId c = s.add_cell(CellKind::Geometry, "G");
    REQUIRE(s.cell(c)->adapter_id().empty());

    s.cell(c)->set_adapter_id("cad.import.step");
    REQUIRE(s.cell(c)->adapter_id() == "cad.import.step");

    s.cell(c)->set_adapter_id("");
    REQUIRE(s.cell(c)->adapter_id().empty());
}
