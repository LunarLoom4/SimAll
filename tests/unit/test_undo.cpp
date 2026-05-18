#include <catch2/catch_test_macros.hpp>
#include "core/Command.hpp"
#include <memory>

using namespace simall::core;
namespace {
struct SetCmd : ICommand {
    int& target; int newV, oldV;
    SetCmd(int& t, int v) : target(t), newV(v), oldV(t) {}
    void execute() override { target = newV; }
    void undo()    override { target = oldV; }
    std::string description() const override { return "Set"; }
};
}

TEST_CASE("CommandHistory undo/redo", "[core]") {
    CommandHistory h;
    int v = 0;
    h.execute(std::make_unique<SetCmd>(v, 10));
    h.execute(std::make_unique<SetCmd>(v, 20));
    REQUIRE(v == 20);
    h.undo(); REQUIRE(v == 10);
    h.undo(); REQUIRE(v == 0);
    h.redo(); REQUIRE(v == 10);
}
