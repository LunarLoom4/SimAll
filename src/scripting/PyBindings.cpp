// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/PyBindings.cpp
//
// When SIMALL_HAVE_PYTHON is defined (CMake found pybind11 + Python), this
// TU embeds CPython via pybind11::scoped_interpreter, registers the
// `simall` module with bindings for ScriptContext, CommandBus, UdfHost,
// and routes stdout/stderr through a Python `_StdoutBridge` class so the
// GUI's PythonConsolePanel sees printf-style output.
//
// When Python is not available, every function is a no-op stub — callers
// see `is_built_in() == false` and bypass the engine entirely.
// =============================================================================
#include "scripting/PyBindings.hpp"
#include "scripting/CommandBus.hpp"
#include "scripting/ScriptContext.hpp"
#include "scripting/UdfHost.hpp"

#include <mutex>
#include <sstream>

namespace simall::scripting::python {

namespace {
std::mutex                              g_outMu;
OutputSink                              g_sink;
std::string                             g_buffer;
}  // namespace

#ifndef SIMALL_HAVE_PYTHON
// -----------------------------------------------------------------------------
// Stub build — no Python interpreter linked.
// -----------------------------------------------------------------------------
bool is_built_in() noexcept { return false; }
bool is_running()  noexcept { return false; }
void start() {}
void stop()  {}
std::string exec(const std::string&) {
    return "Python scripting is not built into this binary "
           "(rebuild with -DSIMALL_ENABLE_PYTHON=ON).";
}
std::string eval(const std::string&) { return exec({}); }
void install_output(OutputSink sink) {
    std::lock_guard lk(g_outMu);
    g_sink = std::move(sink);
}
std::string drain_captured_output() {
    std::lock_guard lk(g_outMu);
    std::string out;
    out.swap(g_buffer);
    return out;
}

#else  // SIMALL_HAVE_PYTHON
// -----------------------------------------------------------------------------
// Real pybind11 implementation.
// -----------------------------------------------------------------------------
#include <pybind11/embed.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
namespace py = pybind11;

namespace {
std::unique_ptr<py::scoped_interpreter> g_interp;
py::object                              g_mainModule;
py::object                              g_mainDict;

void publish(const std::string& s) {
    std::lock_guard lk(g_outMu);
    g_buffer += s;
    if (g_sink) {
        OutputSink local = g_sink;
        lk.~lock_guard();                       // release before user callback
        try { local(s); } catch (...) {}
        return;
    }
}

void register_simall_module(py::module_& m) {
    // ScriptContext bindings.
    py::class_<ScriptContext, std::unique_ptr<ScriptContext, py::nodelete>>(m, "Context")
        .def_static("instance", &ScriptContext::instance,
                    py::return_value_policy::reference)
        .def("get",  [](ScriptContext& s, const std::string& k) {
                          return s.get_or(k, "");
                      })
        .def("set",  &ScriptContext::set)
        .def("erase",&ScriptContext::erase)
        .def("keys", &ScriptContext::keys)
        .def("clear",&ScriptContext::clear);

    // CommandBus bindings.
    py::class_<ScriptableCommand>(m, "Command")
        .def(py::init<>())
        .def_readwrite("name", &ScriptableCommand::name)
        .def_readwrite("args", &ScriptableCommand::args)
        .def("to_python_call", &ScriptableCommand::to_python_call,
             py::arg("module") = "simall");
    py::class_<CommandBus, std::unique_ptr<CommandBus, py::nodelete>>(m, "CommandBus")
        .def_static("instance", &CommandBus::instance,
                    py::return_value_policy::reference)
        .def("dispatch",          &CommandBus::dispatch)
        .def("has_handler",       [](CommandBus& b, const std::string& n) {
                                       return b.has_handler(n);
                                   })
        .def("registered_commands", &CommandBus::registered_commands);

    // UdfHost — accept any Python callable, classify by return type.
    py::class_<UdfHost, std::unique_ptr<UdfHost, py::nodelete>>(m, "UdfHost")
        .def_static("instance", &UdfHost::instance,
                    py::return_value_policy::reference)
        .def("register_scalar_t",
             [](UdfHost& h, const std::string& n, py::function fn) {
                 h.register_udf(n, std::make_shared<ScalarTUdf>(
                     [fn](double t) {
                         py::gil_scoped_acquire g;
                         return fn(t).cast<double>();
                     }));
             })
        .def("register_scalar_txyz",
             [](UdfHost& h, const std::string& n, py::function fn) {
                 h.register_udf(n, std::make_shared<ScalarTXYZUdf>(
                     [fn](double t, double x, double y, double z) {
                         py::gil_scoped_acquire g;
                         return fn(t, x, y, z).cast<double>();
                     }));
             })
        .def("register_vector_txyz",
             [](UdfHost& h, const std::string& n, py::function fn) {
                 h.register_udf(n, std::make_shared<VectorTXYZUdf>(
                     [fn](double t, double x, double y, double z) {
                         py::gil_scoped_acquire g;
                         auto seq = fn(t, x, y, z).cast<std::vector<double>>();
                         std::array<double,3> r{0,0,0};
                         for (size_t i = 0; i < seq.size() && i < 3; ++i) r[i] = seq[i];
                         return r;
                     }));
             })
        .def("has",    [](UdfHost& h, const std::string& n) { return h.has(n); })
        .def("names",  &UdfHost::names)
        .def("clear",  &UdfHost::clear)
        .def("unregister", [](UdfHost& h, const std::string& n) {
                               return h.unregister_udf(n);
                           });

    // Convenience top-level helpers.
    m.def("log", [](const std::string& msg) { publish(msg + "\n"); });
    m.def("ctx", []() -> ScriptContext& { return ScriptContext::instance(); },
          py::return_value_policy::reference);
}

PYBIND11_EMBEDDED_MODULE(simall, m) {
    m.doc() = "SimAll Beta scripting bridge";
    register_simall_module(m);
}

// Replace sys.stdout/sys.stderr with an object that calls back into publish().
void install_stdio_bridge() {
    py::exec(R"PY(
import sys
class _SimAllStdio:
    def __init__(self):  self.buf = ""
    def write(self, s):
        import _simall_stdio_sink
        _simall_stdio_sink.emit(s)
    def flush(self):     pass
sys.stdout = _SimAllStdio()
sys.stderr = _SimAllStdio()
)PY");
}
}  // namespace

bool is_built_in() noexcept { return true; }
bool is_running()  noexcept { return static_cast<bool>(g_interp); }

void start() {
    if (g_interp) return;
    g_interp = std::make_unique<py::scoped_interpreter>();
    // Side-module that exposes a C++ sink to the Python-side stdio bridge.
    auto sinkMod = py::module_::create_extension_module(
        "_simall_stdio_sink", "internal", new py::module_::module_def{});
    sinkMod.def("emit", [](const std::string& s) { publish(s); });
    py::module_::import("sys").attr("modules")["_simall_stdio_sink"] = sinkMod;
    install_stdio_bridge();
    g_mainModule = py::module_::import("__main__");
    g_mainDict   = g_mainModule.attr("__dict__");
    py::exec("import simall", g_mainDict);
}

void stop() {
    if (!g_interp) return;
    g_mainDict = py::object();
    g_mainModule = py::object();
    g_interp.reset();
}

std::string exec(const std::string& code) {
    if (!g_interp) return "Python interpreter not started";
    py::gil_scoped_acquire g;
    try {
        py::exec(code, g_mainDict);
        return {};
    } catch (const py::error_already_set& e) {
        return e.what();
    } catch (const std::exception& e) {
        return e.what();
    }
}

std::string eval(const std::string& expression) {
    if (!g_interp) return "Python interpreter not started";
    py::gil_scoped_acquire g;
    try {
        py::object r = py::eval(expression, g_mainDict);
        return py::repr(r).cast<std::string>();
    } catch (const py::error_already_set& e) {
        publish(std::string(e.what()) + "\n");
        return {};
    } catch (const std::exception& e) {
        publish(std::string(e.what()) + "\n");
        return {};
    }
}

void install_output(OutputSink sink) {
    std::lock_guard lk(g_outMu);
    g_sink = std::move(sink);
}

std::string drain_captured_output() {
    std::lock_guard lk(g_outMu);
    std::string out;
    out.swap(g_buffer);
    return out;
}

#endif  // SIMALL_HAVE_PYTHON

}  // namespace simall::scripting::python
