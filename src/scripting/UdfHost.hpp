// =============================================================================
// SimAll Beta — Scripting subsystem
// File   : src/scripting/UdfHost.hpp
//
// User-Defined Function host.  Solvers query UDFs by name + signature; the
// host dispatches to either a native C++ callable or (when Python is built
// in) a Python callback wrapped in a `PythonUdf`.
//
// Three signatures cover ~all CFD UDF use cases:
//
//   Scalar(t)               — time-varying source / inlet ramp.
//   Scalar(t, x, y, z)      — spatially varying inlet profile, source term.
//   Vector(t, x, y, z)→Vec3 — vector inlet velocity, body force.
//
// PythonUdf lives in PyBindings.cpp (compiled only when Python is enabled).
// =============================================================================
#pragma once

#include <array>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace simall::scripting {

enum class UdfSignature { ScalarT, ScalarTXYZ, VectorTXYZ };

class IUdf {
public:
    virtual ~IUdf() = default;
    [[nodiscard]] virtual UdfSignature signature() const = 0;
    virtual double                  eval_scalar(double t) {
        (void)t; return 0.0;
    }
    virtual double                  eval_scalar(double t, double x, double y, double z) {
        (void)t; (void)x; (void)y; (void)z; return 0.0;
    }
    virtual std::array<double,3>    eval_vector(double t, double x, double y, double z) {
        (void)t; (void)x; (void)y; (void)z; return {0,0,0};
    }
};

// -- adapters for native std::function-based UDFs ----------------------------
class ScalarTUdf : public IUdf {
public:
    using Fn = std::function<double(double)>;
    explicit ScalarTUdf(Fn f) : fn_(std::move(f)) {}
    UdfSignature signature() const override { return UdfSignature::ScalarT; }
    double eval_scalar(double t) override { return fn_ ? fn_(t) : 0.0; }
private:
    Fn fn_;
};
class ScalarTXYZUdf : public IUdf {
public:
    using Fn = std::function<double(double,double,double,double)>;
    explicit ScalarTXYZUdf(Fn f) : fn_(std::move(f)) {}
    UdfSignature signature() const override { return UdfSignature::ScalarTXYZ; }
    double eval_scalar(double t, double x, double y, double z) override {
        return fn_ ? fn_(t, x, y, z) : 0.0;
    }
private:
    Fn fn_;
};
class VectorTXYZUdf : public IUdf {
public:
    using Fn = std::function<std::array<double,3>(double,double,double,double)>;
    explicit VectorTXYZUdf(Fn f) : fn_(std::move(f)) {}
    UdfSignature signature() const override { return UdfSignature::VectorTXYZ; }
    std::array<double,3> eval_vector(double t, double x, double y, double z) override {
        return fn_ ? fn_(t, x, y, z) : std::array<double,3>{0,0,0};
    }
private:
    Fn fn_;
};

class UdfHost {
public:
    static UdfHost& instance();

    void                  register_udf(std::string name, std::shared_ptr<IUdf> udf);
    bool                  unregister_udf(std::string_view name);
    [[nodiscard]] bool    has(std::string_view name) const;
    [[nodiscard]] std::shared_ptr<IUdf> get(std::string_view name) const;
    [[nodiscard]] std::vector<std::string> names() const;
    void                  clear();

    // Convenience eval — returns 0/zero-vector on missing / signature mismatch.
    double                eval_scalar(std::string_view name, double t) const;
    double                eval_scalar(std::string_view name,
                                       double t, double x, double y, double z) const;
    std::array<double,3>  eval_vector(std::string_view name,
                                       double t, double x, double y, double z) const;

private:
    UdfHost() = default;
    mutable std::mutex                                     mu_;
    std::unordered_map<std::string, std::shared_ptr<IUdf>> udfs_;
};

}  // namespace simall::scripting
