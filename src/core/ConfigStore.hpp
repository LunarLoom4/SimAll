// =============================================================================
// SimAll Beta - Core Subsystem
// File   : src/core/ConfigStore.hpp
// Phase  : 1.3 (APPLICATION CORE → Configuration)
//
// Layered configuration store. Four tiers, last-wins precedence:
//   Defaults  → built-in defaults loaded by the application at boot
//   Site      → /etc/simall/site.json (or %PROGRAMDATA%\SimAll\site.json)
//   User      → ~/.config/simall/user.json (or %APPDATA%\SimAll\user.json)
//   CLI       → ad-hoc overrides from argv: --key=value or --key.sub=value
//
// Keys use dotted paths ("solver.linear.tolerance"). Values are typed via
// the small `ConfigValue` variant — bool, int64, double, string, array,
// object — kept independent of nlohmann/json so the public API is stable
// regardless of which backend parses the source files.
// =============================================================================
#pragma once

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>
#include <filesystem>

namespace simall::core {

class ConfigValue;
using ConfigArray  = std::vector<ConfigValue>;
using ConfigObject = std::map<std::string, ConfigValue>;

class ConfigValue {
public:
    using Storage = std::variant<std::monostate,
                                 bool,
                                 std::int64_t,
                                 double,
                                 std::string,
                                 ConfigArray,
                                 ConfigObject>;

    ConfigValue() = default;
    template <class T> ConfigValue(T v) : v_(std::move(v)) {}

    bool isNull()   const noexcept { return std::holds_alternative<std::monostate>(v_); }
    bool isBool()   const noexcept { return std::holds_alternative<bool>(v_); }
    bool isInt()    const noexcept { return std::holds_alternative<std::int64_t>(v_); }
    bool isDouble() const noexcept { return std::holds_alternative<double>(v_); }
    bool isString() const noexcept { return std::holds_alternative<std::string>(v_); }
    bool isArray()  const noexcept { return std::holds_alternative<ConfigArray>(v_); }
    bool isObject() const noexcept { return std::holds_alternative<ConfigObject>(v_); }

    bool                asBool(bool def = false) const;
    std::int64_t        asInt(std::int64_t def = 0) const;
    double              asDouble(double def = 0.0) const;
    std::string         asString(std::string_view def = {}) const;
    const ConfigArray&  asArray() const;
    const ConfigObject& asObject() const;

    const Storage& raw() const noexcept { return v_; }
          Storage& raw()       noexcept { return v_; }

private:
    Storage v_;
};

enum class ConfigLayer : int { Defaults = 0, Site = 1, User = 2, Cli = 3, _Count };

class ConfigStore {
public:
    ConfigStore();
    ~ConfigStore();

    /// Replace one layer wholesale. Pass an object value.
    void setLayer(ConfigLayer layer, ConfigValue root);

    /// Merge an overlay on top of the existing layer (object merge, recursive).
    void mergeLayer(ConfigLayer layer, const ConfigValue& overlay);

    /// Load a JSON file into the given layer. Returns false on parse/IO error.
    bool loadJsonFile(ConfigLayer layer, const std::filesystem::path& file);

    /// Parse CLI overrides of the form `--key.subkey=value`. Unknown flags
    /// (those not starting with `--`) are ignored. Values are coerced to
    /// bool/int/double/string heuristically.
    void applyCli(int argc, const char* const* argv);

    /// Query a dotted-path key. Returns nullopt if absent in every layer.
    std::optional<ConfigValue> get(std::string_view dottedPath) const;

    /// Typed convenience wrappers (return `def` on missing or wrong type).
    bool          getBool  (std::string_view path, bool          def = false) const;
    std::int64_t  getInt   (std::string_view path, std::int64_t  def = 0)     const;
    double        getDouble(std::string_view path, double        def = 0.0)   const;
    std::string   getString(std::string_view path, std::string_view def = {}) const;

    /// Final view: defaults < site < user < cli merged into one object.
    ConfigValue effective() const;

private:
    struct Impl;
    std::unique_ptr<Impl> p_;
};

}  // namespace simall::core
