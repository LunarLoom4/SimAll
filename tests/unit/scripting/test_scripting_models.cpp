// =============================================================================
// SimAll Beta — Scripting unit tests (Week 16)
//
// Exercises every Python-free component of the scripting subsystem:
//
//   * Repl::feed  — bracket / triple-string / block-statement detection.
//   * ScriptContext — get/set/erase/keys/subscriber.
//   * CommandBus — handler routing + listener fan-out + Python repr.
//   * MacroRecorder — listener attach, ordered capture, script emission.
//   * UdfHost — scalar(t), scalar(t,x,y,z), vector(t,x,y,z) registration
//                and signature-aware dispatch fall-back.
//
// Python interpreter is NOT required.  PyBindings.cpp links the stub
// branch when SIMALL_HAVE_PYTHON is undefined and `is_built_in()`
// returns false; tests assert exactly that.
// =============================================================================
#include "scripting/CommandBus.hpp"
#include "scripting/MacroRecorder.hpp"
#include "scripting/PyBindings.hpp"
#include "scripting/Repl.hpp"
#include "scripting/ScriptContext.hpp"
#include "scripting/UdfHost.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fstream>
#include <thread>

using namespace simall::scripting;
using Catch::Matchers::WithinAbs;

// ----------------------------------------------------------------------------
// Repl
// ----------------------------------------------------------------------------
TEST_CASE("Repl: single expression returns Ready", "[scripting][repl]")
{
    Repl r;
    auto fr = r.feed("1 + 2");
    REQUIRE(fr.state == ReplState::Ready);
    REQUIRE(r.take_buffer() == "1 + 2");
}

TEST_CASE("Repl: open bracket requires continuation", "[scripting][repl]")
{
    Repl r;
    REQUIRE(r.feed("x = (1 +").state == ReplState::NeedsMore);
    REQUIRE(r.feed("     2)").state == ReplState::Ready);
}

TEST_CASE("Repl: triple quoted string spans lines", "[scripting][repl]")
{
    Repl r;
    REQUIRE(r.feed("s = \"\"\"hello").state == ReplState::NeedsMore);
    REQUIRE(r.feed("world\"\"\"").state == ReplState::Ready);
}

TEST_CASE("Repl: def block terminated by blank line", "[scripting][repl]")
{
    Repl r;
    REQUIRE(r.feed("def f(x):").state == ReplState::NeedsMore);
    REQUIRE(r.feed("    return x*x").state == ReplState::NeedsMore);
    REQUIRE(r.feed("").state == ReplState::Ready);
}

TEST_CASE("Repl: too many closing brackets is a syntax error", "[scripting][repl]")
{
    Repl r;
    auto fr = r.feed("foo(1, 2))");
    REQUIRE(fr.state == ReplState::SyntaxError);
    REQUIRE_FALSE(fr.message.empty());
    REQUIRE(r.take_buffer().empty());
}

TEST_CASE("Repl: history navigation up/down", "[scripting][repl]")
{
    Repl r;
    r.push_history("one");
    r.push_history("two");
    r.push_history("three");
    REQUIRE(r.history_prev() == "three");
    REQUIRE(r.history_prev() == "two");
    REQUIRE(r.history_prev() == "one");
    REQUIRE(r.history_prev() == "one"); // clamp at top
    REQUIRE(r.history_next() == "two");
    REQUIRE(r.history_next() == "three");
    REQUIRE(r.history_next().empty()); // bottom edge
}

TEST_CASE("Repl: comments and strings do not affect bracket count", "[scripting][repl]")
{
    Repl r;
    auto fr = r.feed("y = '((((' + \")))))\"  # )))");
    REQUIRE(fr.state == ReplState::Ready);
}

// ----------------------------------------------------------------------------
// ScriptContext
// ----------------------------------------------------------------------------
TEST_CASE("ScriptContext: get_or returns fallback on miss", "[scripting][context]")
{
    auto& ctx = ScriptContext::instance();
    ctx.clear();
    REQUIRE(ctx.get_or("missing", "default") == "default");
    ctx.set("project.name", "cavity");
    REQUIRE(ctx.get_or("project.name", "x") == "cavity");
    REQUIRE(ctx.erase("project.name"));
    REQUIRE_FALSE(ctx.erase("project.name"));
}

TEST_CASE("ScriptContext: subscriber receives change events", "[scripting][context]")
{
    auto& ctx = ScriptContext::instance();
    ctx.clear();
    int hits = 0;
    std::string lastKey, lastOld, lastNew;
    int id = ctx.subscribe([&](const std::string& k, const std::string& o, const std::string& n) {
        ++hits;
        lastKey = k;
        lastOld = o;
        lastNew = n;
    });
    ctx.set("a", "1");
    ctx.set("a", "2");
    ctx.unsubscribe(id);
    ctx.set("a", "3"); // no longer counted
    REQUIRE(hits == 2);
    REQUIRE(lastKey == "a");
    REQUIRE(lastOld == "1");
    REQUIRE(lastNew == "2");
}

// ----------------------------------------------------------------------------
// CommandBus
// ----------------------------------------------------------------------------
TEST_CASE("CommandBus: handler runs and listeners fan-out", "[scripting][bus]")
{
    auto& bus = CommandBus::instance();
    bus.clear();

    int handlerCount = 0;
    int listenerCount = 0;
    bus.register_handler("Mesh.Surface", [&](const ScriptableCommand&) { ++handlerCount; });
    int lid = bus.add_listener([&](const ScriptableCommand&) { ++listenerCount; });

    ScriptableCommand c;
    c.name = "Mesh.Surface";
    c.args = {{"size", "0.01"}, {"label", "fine"}};
    REQUIRE(bus.dispatch(c));
    REQUIRE(handlerCount == 1);
    REQUIRE(listenerCount == 1);

    // Unknown command — listeners still fire, handler does not.
    ScriptableCommand u{"Mesh.NotAHandler", {{"x", "1"}}};
    REQUIRE_FALSE(bus.dispatch(u));
    REQUIRE(listenerCount == 2);
    REQUIRE(handlerCount == 1);

    bus.remove_listener(lid);
    bus.clear();
}

TEST_CASE("CommandBus: ScriptableCommand renders Python repr", "[scripting][bus]")
{
    ScriptableCommand c;
    c.name = "BC.Set";
    c.args = {{"zone", "inlet"}, {"u", "2.5"}, {"on", "true"}};
    REQUIRE(c.to_python_call("simall") == "simall.BC.Set(zone=\"inlet\", u=2.5, on=true)");
}

// ----------------------------------------------------------------------------
// MacroRecorder
// ----------------------------------------------------------------------------
TEST_CASE("MacroRecorder: captures only between start and stop", "[scripting][macro]")
{
    auto& bus = CommandBus::instance();
    bus.clear();
    MacroRecorder rec;
    bus.dispatch({"Before.Start", {}}); // not recorded
    rec.start();
    bus.dispatch({"Mesh.Surface", {{"size", "0.01"}}});
    bus.dispatch({"Mesh.Volume", {{"algorithm", "tetra"}}});
    rec.stop();
    bus.dispatch({"After.Stop", {}}); // not recorded

    REQUIRE(rec.size() == 2);
    auto script = rec.to_python_script("regression test");
    REQUIRE(script.find("simall.Mesh.Surface(size=0.01)") != std::string::npos);
    REQUIRE(script.find("simall.Mesh.Volume(algorithm=\"tetra\")") != std::string::npos);
    bus.clear();
}

// ----------------------------------------------------------------------------
// UdfHost
// ----------------------------------------------------------------------------
TEST_CASE("UdfHost: scalar(t) registration and eval", "[scripting][udf]")
{
    auto& h = UdfHost::instance();
    h.clear();
    h.register_udf("ramp", std::make_shared<ScalarTUdf>([](double t) { return 2.0 * t; }));
    REQUIRE(h.has("ramp"));
    REQUIRE_THAT(h.eval_scalar("ramp", 3.0), WithinAbs(6.0, 1e-12));
    REQUIRE(h.unregister_udf("ramp"));
    REQUIRE_FALSE(h.has("ramp"));
}

TEST_CASE("UdfHost: vector(t,x,y,z) signature-aware dispatch", "[scripting][udf]")
{
    auto& h = UdfHost::instance();
    h.clear();
    h.register_udf("inlet", std::make_shared<VectorTXYZUdf>([](double, double, double y, double) {
                       std::array<double, 3> v{1.0 - y * y, 0.0, 0.0};
                       return v;
                   }));
    auto v = h.eval_vector("inlet", 0.0, 0.0, 0.5, 0.0);
    REQUIRE_THAT(v[0], WithinAbs(0.75, 1e-12));
    REQUIRE_THAT(v[1], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(v[2], WithinAbs(0.0, 1e-12));

    // scalar(t,x,y,z) call against vector UDF falls back to 0 (mismatch).
    REQUIRE_THAT(h.eval_scalar("inlet", 0.0, 0.0, 0.5, 0.0), WithinAbs(0.0, 1e-12));
    h.clear();
}

TEST_CASE("UdfHost: missing name returns neutral element", "[scripting][udf]")
{
    auto& h = UdfHost::instance();
    h.clear();
    REQUIRE_THAT(h.eval_scalar("does_not_exist", 1.0), WithinAbs(0.0, 1e-12));
    auto v = h.eval_vector("does_not_exist", 0, 0, 0, 0);
    REQUIRE(v[0] == 0.0);
    REQUIRE(v[1] == 0.0);
    REQUIRE(v[2] == 0.0);
}

// ----------------------------------------------------------------------------
// PyBindings stub mode
// ----------------------------------------------------------------------------
TEST_CASE("PyBindings: stub mode advertises absence of Python", "[scripting][python]")
{
#ifndef SIMALL_HAVE_PYTHON
    REQUIRE_FALSE(python::is_built_in());
    REQUIRE_FALSE(python::is_running());
    auto err = python::exec("x = 1");
    REQUIRE(err.find("not built") != std::string::npos);
#else
    REQUIRE(python::is_built_in());
#endif
}
