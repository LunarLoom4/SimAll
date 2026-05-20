// =============================================================================
// SimAll Beta — Unit tests for GuiCore (Week 15)
// File   : tests/unit/gui_core/test_gui_core_models.cpp
// =============================================================================
#include "gui_core/DockPerspective.hpp"
#include "gui_core/PanelRegistry.hpp"
#include "gui_core/PropertyDescriptor.hpp"
#include "gui_core/ThemePaths.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using namespace simall::gui_core;

TEST_CASE("PropertyDescriptor validates Double bounds", "[gui_core][property]")
{
    PropertyDescriptor d;
    d.propertyName = "p";
    d.type = PropertyType::Double;
    d.minimum = 0.0;
    d.maximum = 10.0;
    d.defaultValue = 1.0;
    d.currentValue = 1.0;

    REQUIRE(d.validate(Variant{5.0}).ok);
    REQUIRE_FALSE(d.validate(Variant{-1.0}).ok);
    REQUIRE_FALSE(d.validate(Variant{11.0}).ok);

    auto r = d.commit(Variant{7.5});
    REQUIRE(r.ok);
    REQUIRE(std::get<double>(d.currentValue) == 7.5);
}

TEST_CASE("PropertyDescriptor enum index range", "[gui_core][property]")
{
    PropertyDescriptor d;
    d.propertyName = "algo";
    d.type = PropertyType::Enum;
    d.enumOptions = {"SIMPLE", "PISO", "Coupled"};
    d.defaultValue = 0;
    d.currentValue = 0;
    REQUIRE(d.validate(Variant{2}).ok);
    REQUIRE_FALSE(d.validate(Variant{3}).ok);
    REQUIRE_FALSE(d.validate(Variant{-1}).ok);
}

TEST_CASE("PropertyDescriptor Range enforces min<=max and Color stays in [0,1]",
          "[gui_core][property]")
{
    PropertyDescriptor rng;
    rng.propertyName = "iterRange";
    rng.type = PropertyType::Range;
    rng.minimum = 0.0;
    rng.maximum = 1e9;
    rng.defaultValue = std::array<double, 3>{0, 1, 0};
    rng.currentValue = std::array<double, 3>{0, 1, 0};
    REQUIRE(rng.validate(Variant{std::array<double, 3>{10, 20, 0}}).ok);
    REQUIRE_FALSE(rng.validate(Variant{std::array<double, 3>{20, 10, 0}}).ok);

    PropertyDescriptor col;
    col.propertyName = "tint";
    col.type = PropertyType::Color;
    col.minimum = 0.0;
    col.maximum = 1.0;
    col.defaultValue = std::array<double, 3>{1, 1, 1};
    col.currentValue = std::array<double, 3>{1, 1, 1};
    REQUIRE(col.validate(Variant{std::array<double, 3>{0.5, 0.5, 0.5}}).ok);
    REQUIRE_FALSE(col.validate(Variant{std::array<double, 3>{1.5, 0, 0}}).ok);
}

TEST_CASE("PropertyDescriptor custom validator runs after built-ins", "[gui_core][property]")
{
    PropertyDescriptor d;
    d.propertyName = "even";
    d.type = PropertyType::Int;
    d.minimum = 0;
    d.maximum = 100;
    d.currentValue = 0;
    d.customValidator = [](const Variant& v) -> ValidationResult {
        const int x = std::get<int>(v);
        return (x % 2 == 0) ? ValidationResult::success()
                            : ValidationResult::failure("must be even");
    };
    REQUIRE(d.validate(Variant{4}).ok);
    auto r = d.validate(Variant{5});
    REQUIRE_FALSE(r.ok);
    REQUIRE_THAT(r.message, Catch::Matchers::ContainsSubstring("even"));
}

TEST_CASE("PropertyDescriptor commit fires onChanged once and updates currentValue",
          "[gui_core][property]")
{
    PropertyDescriptor d;
    d.propertyName = "x";
    d.type = PropertyType::Double;
    d.currentValue = 0.0;
    int hits = 0;
    d.onChanged = [&](const Variant&) { ++hits; };
    REQUIRE(d.commit(Variant{2.5}).ok);
    REQUIRE(hits == 1);
    REQUIRE(std::get<double>(d.currentValue) == 2.5);
    REQUIRE_FALSE(d.commit(Variant{double(NAN)}).ok);
    REQUIRE(hits == 1); // unchanged on failure
}

TEST_CASE("PropertyDescriptor bag find + textual snapshot", "[gui_core][property]")
{
    std::vector<PropertyDescriptor> bag;
    PropertyDescriptor a;
    a.propertyName = "a";
    a.type = PropertyType::Int;
    a.currentValue = 3;
    PropertyDescriptor b;
    b.propertyName = "b";
    b.type = PropertyType::Bool;
    b.currentValue = true;
    bag.push_back(a);
    bag.push_back(b);
    REQUIRE(find(bag, "a") != nullptr);
    REQUIRE(find(bag, "missing") == nullptr);
    const std::string s = to_textual_snapshot(bag);
    REQUIRE_THAT(s, Catch::Matchers::ContainsSubstring("a = 3"));
    REQUIRE_THAT(s, Catch::Matchers::ContainsSubstring("b = true"));
}

TEST_CASE("DockPerspective text round-trip preserves placements", "[gui_core][dock]")
{
    DockPerspective p;
    p.name = "Default";
    p.qtBlob = "Zm9vYmFy"; // arbitrary base64 payload
    DockPlacement a{"Console", DockArea::Bottom, true, false, 0, 0, 0, 0, 1, 0};
    DockPlacement b{"Properties", DockArea::Right, true, false, 0, 0, 0, 0, -1, 0};
    DockPlacement c{"PythonConsole", DockArea::Floating, false, true, 100, 120, 640, 480, -1, 0};
    p.placements = {a, b, c};

    const std::string txt = serialize(p);
    auto p2 = deserialize(txt);
    REQUIRE(p2.has_value());
    REQUIRE(p2->name == "Default");
    REQUIRE(p2->qtBlob == "Zm9vYmFy");
    REQUIRE(p2->placements.size() == 3);
    REQUIRE(p2->find("Console") != nullptr);
    REQUIRE(p2->find("Console")->area == DockArea::Bottom);
    auto* py = p2->find("PythonConsole");
    REQUIRE(py != nullptr);
    REQUIRE(py->floating);
    REQUIRE(py->w == 640);
    REQUIRE(py->h == 480);
}

TEST_CASE("DockPerspective::upsert is idempotent", "[gui_core][dock]")
{
    DockPerspective p;
    auto& first = p.upsert("X");
    first.area = DockArea::Top;
    auto& again = p.upsert("X");
    again.visible = false;
    REQUIRE(p.placements.size() == 1);
    REQUIRE(p.placements[0].area == DockArea::Top);
    REQUIRE_FALSE(p.placements[0].visible);
}

TEST_CASE("PanelRegistry registers and dispatches factories", "[gui_core][panels]")
{
    auto& reg = PanelRegistry::instance();
    reg.clear();
    int sentinel = 0;
    reg.register_panel(PanelInfo{"Spam", "Spam", "Diagnostics", "", true}, [&sentinel]() -> void* {
        ++sentinel;
        return &sentinel;
    });
    REQUIRE(reg.has("Spam"));
    REQUIRE(reg.info("Spam")->displayName == "Spam");
    void* p = reg.create("Spam");
    REQUIRE(p == &sentinel);
    REQUIRE(sentinel == 1);
    REQUIRE(reg.list().size() == 1);
    REQUIRE(reg.create("nope") == nullptr);
    reg.clear();
    REQUIRE_FALSE(reg.has("Spam"));
}

TEST_CASE("ThemePaths gracefully reports empty result for unknown roots", "[gui_core][theme]")
{
    // Pass two nonexistent roots; the loader should return without crashing
    // and the searchedDirectories list should be populated.
    auto loc = locate_theme_assets("/no/such/exe/dir", "/no/such/source/root", "dark");
    REQUIRE_FALSE(loc.searchedDirectories.empty());
    // stylesheet may or may not be found depending on installed prefix.
}
