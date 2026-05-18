// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/core/test_config_store.cpp
// =============================================================================
#include "core/ConfigStore.hpp"

#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <filesystem>
#include <fstream>

using namespace simall::core;

namespace {
std::filesystem::path writeTempJson(const std::string& body) {
    auto p = std::filesystem::temp_directory_path() /
             ("simall_cfg_" + std::to_string(std::rand()) + ".json");
    std::ofstream(p) << body;
    return p;
}
}

TEST_CASE("ConfigStore layered precedence", "[core][config]") {
    ConfigStore cs;

    ConfigObject defaults;
    defaults["solver"] = ConfigValue{ConfigObject{
        {"tol",   ConfigValue{1e-6}},
        {"iter",  ConfigValue{std::int64_t{100}}}}};
    cs.setLayer(ConfigLayer::Defaults, ConfigValue{std::move(defaults)});

    ConfigObject user;
    user["solver"] = ConfigValue{ConfigObject{{"tol", ConfigValue{1e-8}}}};
    cs.setLayer(ConfigLayer::User, ConfigValue{std::move(user)});

    REQUIRE(cs.getDouble("solver.tol") == 1e-8);
    REQUIRE(cs.getInt   ("solver.iter") == 100);
    REQUIRE(cs.getString("missing.key", "fallback") == "fallback");
}

TEST_CASE("ConfigStore loads JSON", "[core][config]") {
    auto p = writeTempJson(R"({"a":{"b":42,"c":[1,2,3]},"flag":true})");
    ConfigStore cs;
    REQUIRE(cs.loadJsonFile(ConfigLayer::Site, p));
    REQUIRE(cs.getInt("a.b") == 42);
    REQUIRE(cs.getBool("flag") == true);
    auto arr = cs.get("a.c");
    REQUIRE(arr.has_value());
    REQUIRE(arr->isArray());
    REQUIRE(arr->asArray().size() == 3);
    std::filesystem::remove(p);
}

TEST_CASE("ConfigStore CLI overrides win", "[core][config]") {
    ConfigStore cs;
    cs.setLayer(ConfigLayer::Defaults,
        ConfigValue{ConfigObject{{"n", ConfigValue{std::int64_t{1}}}}});

    const char* argv[] = {"prog", "--n=99", "--solver.tol=1e-12", "--enabled"};
    cs.applyCli(4, argv);

    REQUIRE(cs.getInt   ("n") == 99);
    REQUIRE(cs.getDouble("solver.tol") == 1e-12);
    REQUIRE(cs.getBool  ("enabled"));
}
