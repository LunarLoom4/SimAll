// =============================================================================
// SimAll Beta — Workbench Unit Tests
// File   : tests/unit/workbench/test_workbench_commands.cpp
// Phase  : 22 Pass 22.4
//
// Reversibility / journaling tests for the six workbench command
// factories.  Each command is exercised through core::CommandHistory so
// the do/undo/redo cycle exactly mirrors how the UI will drive them.
// =============================================================================
#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/ChangeJournal.hpp"
#include "workbench/Commands.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"
#include "workbench/WorkflowEngine.hpp"

#include "core/Command.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <memory>
#include <string>

using namespace simall::workbench;
using simall::core::CommandHistory;

namespace {

// Minimal "Geometry -> Mesh" fixture.  Two cells, no links by default; the
// individual tests wire links as needed.
struct Fixture {
    Schematic     s;
    StateMachine  sm{s};
    WorkflowEngine eng{s, sm};
    ChangeJournal journal;
    CommandHistory hist;

    CellId g{}, m{};
    PortId g_out{}, m_in{};

    Fixture() {
        g = s.add_cell(CellKind::Geometry, "G");
        m = s.add_cell(CellKind::Mesh,     "M");
        g_out = s.cell(g)->add_output("x", "x");
        m_in  = s.cell(m)->add_input ("x", "x");
    }
};

}  // namespace


// ---------------------------------------------------------------------------
// ChangeJournal
// ---------------------------------------------------------------------------
TEST_CASE("ChangeJournal records ordered Do/Undo entries and respects capacity",
          "[workbench][journal]") {
    ChangeJournal j(3);
    j.record(ChangeDirection::Do,   "a");
    j.record(ChangeDirection::Do,   "b");
    j.record(ChangeDirection::Undo, "a");
    j.record(ChangeDirection::Do,   "c");
    const auto entries = j.snapshot();
    REQUIRE(entries.size() == 3);
    REQUIRE(entries.front().description == "b");
    REQUIRE(entries.back ().description == "c");
    REQUIRE(entries.front().sequence    == 1);
    REQUIRE(entries.back ().sequence    == 3);
    REQUIRE(entries[1].direction        == ChangeDirection::Undo);
}

TEST_CASE("ChangeJournal sink is invoked synchronously on every record",
          "[workbench][journal][sink]") {
    ChangeJournal j;
    std::vector<std::string> seen;
    j.set_sink([&seen](const ChangeEntry& e) { seen.push_back(e.description); });
    j.record(ChangeDirection::Do, "alpha");
    j.record(ChangeDirection::Do, "beta");
    REQUIRE(seen == std::vector<std::string>{"alpha", "beta"});
}


// ---------------------------------------------------------------------------
// AddCellCommand
// ---------------------------------------------------------------------------
TEST_CASE("AddCellCommand do/undo/redo round-trips a new cell with the same id",
          "[workbench][commands][add-cell]") {
    Fixture f;
    CellId minted = kInvalidCellId;
    f.hist.execute(make_add_cell_command(f.s, CellKind::Setup, "Su",
                                          &minted, &f.journal));
    REQUIRE(minted != kInvalidCellId);
    REQUIRE(f.s.cell(minted) != nullptr);
    REQUIRE(f.s.size() == 3);

    f.hist.undo();
    REQUIRE(f.s.cell(minted) == nullptr);
    REQUIRE(f.s.size() == 2);

    f.hist.redo();
    REQUIRE(f.s.cell(minted) != nullptr);
    REQUIRE(f.s.cell(minted)->label() == "Su");
    REQUIRE(f.s.size() == 3);

    const auto entries = f.journal.snapshot();
    REQUIRE(entries.size() == 3);
    REQUIRE(entries[0].direction == ChangeDirection::Do);
    REQUIRE(entries[1].direction == ChangeDirection::Undo);
    REQUIRE(entries[2].direction == ChangeDirection::Do);
    REQUIRE(entries[0].description.find("Add Setup cell 'Su'") != std::string::npos);
}


// ---------------------------------------------------------------------------
// RemoveCellCommand: cell + every incident link survives undo verbatim.
// ---------------------------------------------------------------------------
TEST_CASE("RemoveCellCommand restores cell ports + every incident link on undo",
          "[workbench][commands][remove-cell]") {
    Fixture f;
    // Wire G.out -> M.in so removal must drop one link.
    f.hist.execute(make_connect_command(f.eng, CellLink{f.g, f.g_out, f.m, f.m_in},
                                        &f.journal));
    REQUIRE(f.s.links().size() == 1);

    f.hist.execute(make_remove_cell_command(f.s, f.m, &f.journal));
    REQUIRE(f.s.cell(f.m) == nullptr);
    REQUIRE(f.s.links().empty());

    f.hist.undo();    // re-add M and the link
    REQUIRE(f.s.cell(f.m) != nullptr);
    REQUIRE(f.s.cell(f.m)->label() == "M");
    REQUIRE(f.s.cell(f.m)->ports().size() == 1);
    REQUIRE(f.s.cell(f.m)->ports().front().name == "x");
    REQUIRE(f.s.links().size() == 1);

    // Redo (remove again) then undo once more -- snapshot must survive.
    f.hist.redo();
    REQUIRE(f.s.cell(f.m) == nullptr);
    REQUIRE(f.s.links().empty());
    f.hist.undo();
    REQUIRE(f.s.cell(f.m) != nullptr);
    REQUIRE(f.s.links().size() == 1);
}


// ---------------------------------------------------------------------------
// ConnectCommand drives WorkflowEngine -> state propagation fires.
// ---------------------------------------------------------------------------
TEST_CASE("ConnectCommand fires downstream invalidation; undo restores",
          "[workbench][commands][connect]") {
    Fixture f;
    f.sm.mark_solved(f.g);
    f.sm.mark_solved(f.m);   // M is unwired so will be Unfulfilled afterwards
    // Manually move M into UpToDate via state-set command so we can prove
    // connect() demotes it on execute.
    f.hist.execute(make_set_cell_state_command(f.s, f.m, CellState::UpToDate,
                                               &f.journal));
    REQUIRE(f.s.cell(f.m)->state() == CellState::UpToDate);

    f.hist.execute(make_connect_command(f.eng,
                                        CellLink{f.g, f.g_out, f.m, f.m_in},
                                        &f.journal));
    REQUIRE(f.s.links().size() == 1);
    // The sink (M) was UpToDate; connect() invalidates it.
    REQUIRE(f.s.cell(f.m)->state() != CellState::UpToDate);

    f.hist.undo();   // disconnect
    REQUIRE(f.s.links().empty());

    f.hist.redo();
    REQUIRE(f.s.links().size() == 1);
}


// ---------------------------------------------------------------------------
// DisconnectCommand on an established link.
// ---------------------------------------------------------------------------
TEST_CASE("DisconnectCommand undo re-adds the same link",
          "[workbench][commands][disconnect]") {
    Fixture f;
    REQUIRE(f.eng.connect({f.g, f.g_out, f.m, f.m_in}));
    REQUIRE(f.s.links().size() == 1);

    f.hist.execute(make_disconnect_command(f.eng,
                                           CellLink{f.g, f.g_out, f.m, f.m_in},
                                           &f.journal));
    REQUIRE(f.s.links().empty());
    f.hist.undo();
    REQUIRE(f.s.links().size() == 1);
}


// ---------------------------------------------------------------------------
// SetCellLabelCommand
// ---------------------------------------------------------------------------
TEST_CASE("SetCellLabelCommand is symmetric",
          "[workbench][commands][label]") {
    Fixture f;
    f.hist.execute(make_set_cell_label_command(f.s, f.g, "Imported geometry",
                                               &f.journal));
    REQUIRE(f.s.cell(f.g)->label() == "Imported geometry");
    f.hist.undo();
    REQUIRE(f.s.cell(f.g)->label() == "G");
    f.hist.redo();
    REQUIRE(f.s.cell(f.g)->label() == "Imported geometry");
}


// ---------------------------------------------------------------------------
// SetCellStateCommand
// ---------------------------------------------------------------------------
TEST_CASE("SetCellStateCommand records prior state for undo",
          "[workbench][commands][state]") {
    Fixture f;
    REQUIRE(f.s.cell(f.g)->state() == CellState::Unfulfilled);
    f.hist.execute(make_set_cell_state_command(f.s, f.g, CellState::UpToDate,
                                               &f.journal));
    REQUIRE(f.s.cell(f.g)->state() == CellState::UpToDate);
    f.hist.undo();
    REQUIRE(f.s.cell(f.g)->state() == CellState::Unfulfilled);
}


// ---------------------------------------------------------------------------
// CompositeCommand: bundle add-cell + add-port-less-skip + connect into one
// undoable user action.  Validates that core::CompositeCommand cooperates
// with workbench commands (rollback on mid-sequence failure).
// ---------------------------------------------------------------------------
TEST_CASE("CompositeCommand rolls back partial workbench mutations on failure",
          "[workbench][commands][composite]") {
    Fixture f;

    auto composite = std::make_unique<simall::core::CompositeCommand>("Add+Connect");
    CellId minted = kInvalidCellId;
    composite->add(make_add_cell_command(f.s, CellKind::Setup, "Su",
                                         &minted, &f.journal));
    // The follow-up Connect references a port id that does NOT exist on
    // the freshly-minted Setup cell.  ConnectCommand::execute() throws,
    // which triggers Composite's rollback: minted cell must be removed.
    composite->add(make_connect_command(f.eng,
        CellLink{f.g, f.g_out, /*to_cell*/ 0xDEADBEEF, /*to_port*/ 0},
        &f.journal));

    REQUIRE_THROWS(f.hist.execute(std::move(composite)));
    // After rollback the schematic must look untouched.
    REQUIRE(f.s.size() == 2);
    REQUIRE(f.s.cell(minted) == nullptr);
    REQUIRE(f.s.links().empty());
}


// ---------------------------------------------------------------------------
// Schematic::restore_cell rejects id collisions.
// ---------------------------------------------------------------------------
TEST_CASE("Schematic::restore_cell refuses to clobber an existing id",
          "[workbench][schematic][restore]") {
    Schematic s;
    const CellId id = s.add_cell(CellKind::Geometry, "G");
    Cell duplicate(id, CellKind::Mesh, "Duplicate");
    REQUIRE_FALSE(s.restore_cell(std::move(duplicate)));
    REQUIRE(s.size() == 1);
    REQUIRE(s.cell(id)->label() == "G");   // original untouched
}
