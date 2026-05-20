// =============================================================================
// SimAll Beta — Workbench Unit Tests
// File   : tests/unit/workbench/test_workbench_models.cpp
// Phase  : 22 Pass 22.1
//
// Exercises the pure-data workbench layer (no UI):
//   * Cell + CellPort: id minting, find_port, label + state mutators.
//   * Schematic: cycle prevention, type checking, single fan-in on inputs,
//     duplicate rejection, topological order (deterministic ties),
//     upstream / downstream queries, inputs_satisfied predicate.
//   * StateMachine: mark_modified / mark_solved / mark_failed propagation,
//     recompute_all, Failed-is-sticky semantics.
// =============================================================================
#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/CellPort.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>

using namespace simall::workbench;

namespace {

// Build the canonical Fluent-style Geometry -> Mesh -> Setup -> Solution
// -> Results pipeline.  Returns the five cell ids in pipeline order.
struct Pipeline {
    Schematic s;
    CellId geometry{}, mesh{}, setup{}, solution{}, results{};
};

Pipeline make_fluent_pipeline() {
    Pipeline p;

    p.geometry = p.s.add_cell(CellKind::Geometry, "Geometry");
    p.mesh     = p.s.add_cell(CellKind::Mesh,     "Mesh");
    p.setup    = p.s.add_cell(CellKind::Setup,    "Setup");
    p.solution = p.s.add_cell(CellKind::Solution, "Solution");
    p.results  = p.s.add_cell(CellKind::Results,  "Results");

    Cell* g = p.s.cell(p.geometry);
    Cell* m = p.s.cell(p.mesh);
    Cell* su= p.s.cell(p.setup);
    Cell* so= p.s.cell(p.solution);
    Cell* r = p.s.cell(p.results);

    const PortId g_out  = g ->add_output("geometry", "geometry");
    const PortId m_in   = m ->add_input ("geometry", "geometry");
    const PortId m_out  = m ->add_output("mesh",     "mesh");
    const PortId su_in  = su->add_input ("mesh",     "mesh");
    const PortId su_out = su->add_output("case",     "case");
    const PortId so_in  = so->add_input ("case",     "case");
    const PortId so_out = so->add_output("fields",   "fields");
    const PortId r_in   = r ->add_input ("fields",   "fields");

    REQUIRE(p.s.add_link({p.geometry, g_out,  p.mesh,     m_in }));
    REQUIRE(p.s.add_link({p.mesh,     m_out,  p.setup,    su_in}));
    REQUIRE(p.s.add_link({p.setup,    su_out, p.solution, so_in}));
    REQUIRE(p.s.add_link({p.solution, so_out, p.results,  r_in }));
    return p;
}

}  // namespace


// ---------------------------------------------------------------------------
// Enum stringification
// ---------------------------------------------------------------------------
TEST_CASE("workbench enums stringify", "[workbench][enums]") {
    REQUIRE(to_string(CellKind::Geometry)         == "Geometry");
    REQUIRE(to_string(CellKind::Custom)           == "Custom");
    REQUIRE(to_string(CellState::Unfulfilled)     == "Unfulfilled");
    REQUIRE(to_string(CellState::RefreshRequired) == "RefreshRequired");
    REQUIRE(to_string(CellState::UpToDate)        == "UpToDate");
    REQUIRE(to_string(CellState::Failed)          == "Failed");
    REQUIRE(to_string(PortDirection::Input)       == "Input");
    REQUIRE(to_string(PortDirection::Output)      == "Output");
}


// ---------------------------------------------------------------------------
// Cell basics
// ---------------------------------------------------------------------------
TEST_CASE("Cell mints distinct port ids and find_port returns them",
          "[workbench][cell]") {
    Cell c(/*id=*/7, CellKind::Mesh, "Mesh");

    REQUIRE(c.id()    == 7);
    REQUIRE(c.kind()  == CellKind::Mesh);
    REQUIRE(c.label() == "Mesh");
    REQUIRE(c.state() == CellState::Unfulfilled);

    const PortId in_a  = c.add_input ("geo",  "geometry");
    const PortId in_b  = c.add_input ("opt",  "options", /*required=*/false);
    const PortId out_a = c.add_output("mesh", "mesh");

    REQUIRE(in_a  != in_b);
    REQUIRE(in_b  != out_a);

    const CellPort* p_in_a = c.find_port(in_a);
    REQUIRE(p_in_a);
    REQUIRE(p_in_a->name      == "geo");
    REQUIRE(p_in_a->direction == PortDirection::Input);
    REQUIRE(p_in_a->required  == true);

    const CellPort* p_in_b = c.find_port(in_b);
    REQUIRE(p_in_b);
    REQUIRE(p_in_b->required == false);

    const CellPort* p_out_a = c.find_port(out_a);
    REQUIRE(p_out_a);
    REQUIRE(p_out_a->direction == PortDirection::Output);

    REQUIRE(c.find_port(kInvalidPortId) == nullptr);
}


// ---------------------------------------------------------------------------
// Schematic link validation
// ---------------------------------------------------------------------------
TEST_CASE("Schematic rejects ill-formed links", "[workbench][schematic][links]") {
    Schematic s;
    const CellId a = s.add_cell(CellKind::Geometry, "A");
    const CellId b = s.add_cell(CellKind::Mesh,     "B");

    Cell* ca = s.cell(a);
    Cell* cb = s.cell(b);
    const PortId a_out = ca->add_output("geo", "geometry");
    const PortId a_in  = ca->add_input ("scratch", "scratch");
    const PortId b_in  = cb->add_input ("geo", "geometry");
    const PortId b_in2 = cb->add_input ("alt", "geometry");
    const PortId b_out = cb->add_output("mesh", "mesh");

    SECTION("unknown cell or port") {
        REQUIRE_FALSE(s.add_link({99,  a_out,           b,   b_in}));
        REQUIRE_FALSE(s.add_link({a,   kInvalidPortId,  b,   b_in}));
        REQUIRE_FALSE(s.add_link({a,   a_out,           99,  b_in}));
        REQUIRE_FALSE(s.add_link({a,   a_out,           b,   kInvalidPortId}));
    }
    SECTION("direction errors: input as source / output as sink / self") {
        REQUIRE_FALSE(s.add_link({a, a_in,  b, b_in }));   // source is Input
        REQUIRE_FALSE(s.add_link({a, a_out, b, b_out}));   // sink is Output
    }
    SECTION("data_type mismatch") {
        // a_out is "geometry"; build a "mesh" sink port and try to link it.
        const PortId b_mesh_in = cb->add_input("bogus", "mesh");
        REQUIRE_FALSE(s.add_link({a, a_out, b, b_mesh_in}));
    }
    SECTION("duplicate edge rejected on second add") {
        REQUIRE(s.add_link({a, a_out, b, b_in}));
        REQUIRE_FALSE(s.add_link({a, a_out, b, b_in}));
    }
    SECTION("input ports have fan-in of 1") {
        const CellId c = s.add_cell(CellKind::Geometry, "C");
        const PortId c_out = s.cell(c)->add_output("geo", "geometry");
        REQUIRE(s.add_link({a, a_out, b, b_in}));
        REQUIRE_FALSE(s.add_link({c, c_out, b, b_in}));
        // ...but a *different* input port on b is fine.
        REQUIRE(s.add_link({c, c_out, b, b_in2}));
    }
}


TEST_CASE("Schematic rejects cycles", "[workbench][schematic][cycle]") {
    Schematic s;
    const CellId a = s.add_cell(CellKind::Custom, "A");
    const CellId b = s.add_cell(CellKind::Custom, "B");
    const CellId c = s.add_cell(CellKind::Custom, "C");

    const PortId a_o = s.cell(a)->add_output("x", "x");
    const PortId a_i = s.cell(a)->add_input ("x", "x");
    const PortId b_o = s.cell(b)->add_output("x", "x");
    const PortId b_i = s.cell(b)->add_input ("x", "x");
    const PortId c_o = s.cell(c)->add_output("x", "x");
    const PortId c_i = s.cell(c)->add_input ("x", "x");

    REQUIRE      (s.add_link({a, a_o, b, b_i}));   // A -> B
    REQUIRE      (s.add_link({b, b_o, c, c_i}));   // B -> C
    REQUIRE_FALSE(s.add_link({c, c_o, a, a_i}));   // C -> A would cycle
    REQUIRE_FALSE(s.add_link({a, a_o, a, a_i}));   // self-loop
}


TEST_CASE("Schematic topological_order is deterministic",
          "[workbench][schematic][topo]") {
    auto p = make_fluent_pipeline();
    const auto order = p.s.topological_order();
    REQUIRE(order.size() == 5);
    REQUIRE(order[0] == p.geometry);
    REQUIRE(order[1] == p.mesh);
    REQUIRE(order[2] == p.setup);
    REQUIRE(order[3] == p.solution);
    REQUIRE(order[4] == p.results);
}


TEST_CASE("Schematic upstream / downstream / inputs_satisfied",
          "[workbench][schematic][queries]") {
    auto p = make_fluent_pipeline();

    REQUIRE(p.s.upstream(p.geometry).empty());
    REQUIRE(p.s.upstream(p.setup)   == std::vector<CellId>{p.mesh});
    REQUIRE(p.s.downstream(p.mesh)  == std::vector<CellId>{p.setup});
    REQUIRE(p.s.downstream(p.results).empty());

    REQUIRE(p.s.inputs_satisfied(p.geometry));  // no required inputs
    REQUIRE(p.s.inputs_satisfied(p.mesh));
    REQUIRE(p.s.inputs_satisfied(p.solution));

    // Drop the geometry -> mesh link: mesh's required input is now unwired.
    REQUIRE(p.s.remove_link({p.geometry,
                              p.s.cell(p.geometry)->ports()[0].id,
                              p.mesh,
                              p.s.cell(p.mesh)->ports()[0].id}));
    REQUIRE_FALSE(p.s.inputs_satisfied(p.mesh));
}


// ---------------------------------------------------------------------------
// StateMachine propagation
// ---------------------------------------------------------------------------
TEST_CASE("StateMachine recompute_all derives initial states",
          "[workbench][state][recompute]") {
    auto         p  = make_fluent_pipeline();
    StateMachine sm(p.s);
    sm.recompute_all();

    // Nobody has been "solved" yet, so every cell with satisfied inputs
    // sits at RefreshRequired (== "ready, but stale").
    for (const CellId id : {p.geometry, p.mesh, p.setup, p.solution, p.results}) {
        REQUIRE(p.s.cell(id)->state() == CellState::RefreshRequired);
    }
}


TEST_CASE("StateMachine mark_solved walks downstream to RefreshRequired",
          "[workbench][state][solve]") {
    auto         p  = make_fluent_pipeline();
    StateMachine sm(p.s);
    sm.recompute_all();

    sm.mark_solved(p.geometry);
    REQUIRE(p.s.cell(p.geometry)->state() == CellState::UpToDate);
    // Everything downstream stays RefreshRequired.
    REQUIRE(p.s.cell(p.mesh)    ->state() == CellState::RefreshRequired);
    REQUIRE(p.s.cell(p.setup)   ->state() == CellState::RefreshRequired);

    sm.mark_solved(p.mesh);
    sm.mark_solved(p.setup);
    sm.mark_solved(p.solution);
    sm.mark_solved(p.results);
    for (const CellId id : {p.geometry, p.mesh, p.setup, p.solution, p.results}) {
        REQUIRE(p.s.cell(id)->state() == CellState::UpToDate);
    }
}


TEST_CASE("StateMachine mark_modified invalidates the downstream subgraph",
          "[workbench][state][modify]") {
    auto         p  = make_fluent_pipeline();
    StateMachine sm(p.s);
    sm.recompute_all();
    for (const CellId id : {p.geometry, p.mesh, p.setup, p.solution, p.results}) {
        sm.mark_solved(id);
    }
    REQUIRE(p.s.cell(p.results)->state() == CellState::UpToDate);

    // Edit the geometry -- everything downstream goes stale.
    sm.mark_modified(p.geometry);
    REQUIRE(p.s.cell(p.geometry)->state() == CellState::RefreshRequired);
    REQUIRE(p.s.cell(p.mesh)    ->state() == CellState::RefreshRequired);
    REQUIRE(p.s.cell(p.setup)   ->state() == CellState::RefreshRequired);
    REQUIRE(p.s.cell(p.solution)->state() == CellState::RefreshRequired);
    REQUIRE(p.s.cell(p.results) ->state() == CellState::RefreshRequired);
}


TEST_CASE("StateMachine mark_failed poisons downstream to Unfulfilled",
          "[workbench][state][fail]") {
    auto         p  = make_fluent_pipeline();
    StateMachine sm(p.s);
    sm.recompute_all();
    for (const CellId id : {p.geometry, p.mesh, p.setup, p.solution, p.results}) {
        sm.mark_solved(id);
    }

    sm.mark_failed(p.setup);
    REQUIRE(p.s.cell(p.geometry)->state() == CellState::UpToDate);    // upstream unaffected
    REQUIRE(p.s.cell(p.mesh)    ->state() == CellState::UpToDate);
    REQUIRE(p.s.cell(p.setup)   ->state() == CellState::Failed);
    REQUIRE(p.s.cell(p.solution)->state() == CellState::Unfulfilled);
    REQUIRE(p.s.cell(p.results) ->state() == CellState::Unfulfilled);

    // Failed is sticky: re-solving an ancestor does NOT clear it.
    sm.mark_solved(p.mesh);
    REQUIRE(p.s.cell(p.setup)->state() == CellState::Failed);

    // Only the user re-asserting the cell clears Failed.
    sm.mark_modified(p.setup);
    REQUIRE(p.s.cell(p.setup)->state() == CellState::RefreshRequired);
}


TEST_CASE("Unwired required input keeps a cell Unfulfilled after recompute",
          "[workbench][state][unfulfilled]") {
    Schematic s;
    const CellId a = s.add_cell(CellKind::Geometry, "A");
    const CellId b = s.add_cell(CellKind::Mesh,     "B");

    s.cell(a)->add_output("geo", "geometry");
    s.cell(b)->add_input ("geo", "geometry");   // required, unwired

    StateMachine sm(s);
    sm.recompute_all();
    REQUIRE(s.cell(a)->state() == CellState::RefreshRequired);
    REQUIRE(s.cell(b)->state() == CellState::Unfulfilled);
}
