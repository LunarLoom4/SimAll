// =============================================================================
// SimAll Beta - Tests
// File   : tests/unit/core/test_composite_command.cpp
// =============================================================================
#include "core/Command.hpp"

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>
#include <vector>

using namespace simall::core;

namespace {

class IncrementCmd : public ICommand {
public:
    IncrementCmd(int& slot, int delta, bool failOnExecute = false)
        : slot_(slot), delta_(delta), fail_(failOnExecute) {}
    void execute() override {
        if (fail_) throw std::runtime_error("inject");
        slot_ += delta_;
    }
    void undo() override { slot_ -= delta_; }
    std::string description() const override { return "Increment"; }
private:
    int& slot_;
    int  delta_;
    bool fail_;
};

}  // namespace

TEST_CASE("CompositeCommand executes children in order and undoes in reverse",
          "[core][command]") {
    int x = 0;
    auto comp = CompositeCommand::of("Adds",
        std::make_unique<IncrementCmd>(x, 1),
        std::make_unique<IncrementCmd>(x, 2),
        std::make_unique<IncrementCmd>(x, 3));
    REQUIRE(comp->size() == 3);
    comp->execute();
    REQUIRE(x == 6);
    comp->undo();
    REQUIRE(x == 0);
}

TEST_CASE("CompositeCommand rolls back on partial failure", "[core][command]") {
    int x = 0;
    auto comp = std::make_unique<CompositeCommand>("WithFailure");
    comp->add(std::make_unique<IncrementCmd>(x, 5));
    comp->add(std::make_unique<IncrementCmd>(x, 7));
    comp->add(std::make_unique<IncrementCmd>(x, 0, /*failOnExecute=*/true));
    REQUIRE_THROWS(comp->execute());
    REQUIRE(x == 0);   // both successes were rolled back
}

TEST_CASE("LambdaCommand integrates with history", "[core][command]") {
    int x = 10;
    CommandHistory h;
    h.execute(std::make_unique<LambdaCommand>(
        "set", [&]{ x = 42; }, [&]{ x = 10; }));
    REQUIRE(x == 42);
    h.undo();
    REQUIRE(x == 10);
    h.redo();
    REQUIRE(x == 42);
}
