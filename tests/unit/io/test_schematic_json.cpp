// =============================================================================
// SimAll Beta — IO Unit Tests
// File   : tests/unit/io/test_schematic_json.cpp
// Phase  : 22 Pass 22.3
//
// Round-trip persistence for workbench::Schematic via io::SchematicJson.
// =============================================================================
#include "io/SchematicJson.hpp"

#include "workbench/Cell.hpp"
#include "workbench/CellLink.hpp"
#include "workbench/CellPort.hpp"
#include "workbench/Schematic.hpp"
#include "workbench/StateMachine.hpp"
#include "workbench/Workbench.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace simall::workbench;
using simall::io::load_schematic_json;
using simall::io::save_schematic_json;
using simall::io::schematic_from_json;
using simall::io::schematic_to_json;

namespace {

// Build the canonical 5-cell Workbench pipeline with one optional input on
// the Setup cell so we can prove the `required` flag round-trips.
struct Pipe {
    Schematic s;
    StateMachine sm{s};
    CellId g{}, m{}, su{}, so{}, r{};
    PortId g_o{}, m_i{}, m_o{}, su_i{}, su_opt{}, su_o{}, so_i{}, so_o{}, r_i{};

    void build() {
        g  = s.add_cell(CellKind::Geometry, "G");
        m  = s.add_cell(CellKind::Mesh,     "M");
        su = s.add_cell(CellKind::Setup,    "Su");
        so = s.add_cell(CellKind::Solution, "So");
        r  = s.add_cell(CellKind::Results,  "R");
        g_o    = s.cell(g) ->add_output("x",   "x");
        m_i    = s.cell(m) ->add_input ("x",   "x");
        m_o    = s.cell(m) ->add_output("x",   "x");
        su_i   = s.cell(su)->add_input ("x",   "x");
        su_opt = s.cell(su)->add_input ("opt", "k", /*required*/false);
        su_o   = s.cell(su)->add_output("x",   "x");
        so_i   = s.cell(so)->add_input ("x",   "x");
        so_o   = s.cell(so)->add_output("x",   "x");
        r_i    = s.cell(r) ->add_input ("x",   "x");
        REQUIRE(s.add_link({g,  g_o,  m,  m_i }));
        REQUIRE(s.add_link({m,  m_o,  su, su_i}));
        REQUIRE(s.add_link({su, su_o, so, so_i}));
        REQUIRE(s.add_link({so, so_o, r,  r_i }));
    }
};

// Stable per-cell structural fingerprint (independent of minted ids).
struct CellFP {
    CellKind    kind{};
    std::string label;
    CellState   state{};
    struct PortFP {
        PortDirection dir{};
        std::string   name;
        std::string   type;
        bool          required{};
        bool operator==(const PortFP&) const = default;
    };
    std::vector<PortFP> ports;
    bool operator==(const CellFP&) const = default;
};

std::vector<CellFP> fingerprint_cells(const Schematic& s) {
    std::vector<CellFP> out;
    for (const Cell& c : s.cells()) {
        CellFP cf{c.kind(), c.label(), c.state(), {}};
        for (const CellPort& p : c.ports()) {
            cf.ports.push_back({p.direction, p.name, p.data_type, p.required});
        }
        out.push_back(std::move(cf));
    }
    return out;
}

// Topology fingerprint by *cell-label* pairs (id-independent).
struct LinkFP {
    std::string from_label, to_label;
    std::string from_port,  to_port;
    bool operator==(const LinkFP&) const = default;
};
std::vector<LinkFP> fingerprint_links(const Schematic& s) {
    std::vector<LinkFP> out;
    for (const CellLink& l : s.links()) {
        const Cell* fc = s.cell(l.from_cell);
        const Cell* tc = s.cell(l.to_cell);
        REQUIRE(fc); REQUIRE(tc);
        const CellPort* fp = fc->find_port(l.from_port);
        const CellPort* tp = tc->find_port(l.to_port);
        REQUIRE(fp); REQUIRE(tp);
        out.push_back({fc->label(), tc->label(), fp->name, tp->name});
    }
    std::sort(out.begin(), out.end(),
              [](const LinkFP& a, const LinkFP& b) {
                  return std::tie(a.from_label, a.from_port, a.to_label, a.to_port)
                       < std::tie(b.from_label, b.from_port, b.to_label, b.to_port);
              });
    return out;
}

}  // namespace


// ---------------------------------------------------------------------------
// to_json() emits stable content the parser can handle.
// ---------------------------------------------------------------------------
TEST_CASE("schematic_to_json emits the v1 schema marker and lists every cell + link",
          "[io][schematic-json][emit]") {
    Pipe p; p.build();
    const std::string j = schematic_to_json(p.s);
    REQUIRE(j.find(R"("schema":"simall.workbench.schematic/v1")") != std::string::npos);
    REQUIRE(j.find(R"("label":"G")") != std::string::npos);
    REQUIRE(j.find(R"("kind":"Setup")") != std::string::npos);
    REQUIRE(j.find(R"("dir":"Input")")  != std::string::npos);
    REQUIRE(j.find(R"("required":false)") != std::string::npos);   // su_opt
    REQUIRE(j.find(R"("from_cell":)") != std::string::npos);
}


// ---------------------------------------------------------------------------
// In-memory round trip preserves structure + states.
// ---------------------------------------------------------------------------
TEST_CASE("schematic_from_json reproduces structure, ports, and states",
          "[io][schematic-json][round-trip]") {
    Pipe p; p.build();
    p.sm.mark_solved(p.g);
    p.sm.mark_solved(p.m);
    p.sm.mark_failed(p.su);

    const std::string j = schematic_to_json(p.s);

    Schematic loaded;
    StateMachine loaded_sm{loaded};
    REQUIRE(schematic_from_json(j, loaded, &loaded_sm));

    REQUIRE(fingerprint_cells(p.s) == fingerprint_cells(loaded));
    REQUIRE(fingerprint_links(p.s) == fingerprint_links(loaded));

    // Re-encode and re-decode -> must be a fixed point.
    const std::string j2 = schematic_to_json(loaded);
    Schematic loaded2;
    REQUIRE(schematic_from_json(j2, loaded2));
    REQUIRE(fingerprint_cells(loaded)  == fingerprint_cells(loaded2));
    REQUIRE(fingerprint_links(loaded)  == fingerprint_links(loaded2));
    REQUIRE(j == j2);
}


// ---------------------------------------------------------------------------
// File round trip via temp path.
// ---------------------------------------------------------------------------
TEST_CASE("save_schematic_json + load_schematic_json round-trip on disk",
          "[io][schematic-json][file]") {
    Pipe p; p.build();
    p.sm.mark_solved(p.g);

    const auto path = (std::filesystem::temp_directory_path()
                       / "simall_schematic_roundtrip.json").string();
    REQUIRE(save_schematic_json(path, p.s));

    Schematic loaded;
    REQUIRE(load_schematic_json(path, loaded));
    REQUIRE(fingerprint_cells(p.s) == fingerprint_cells(loaded));
    REQUIRE(fingerprint_links(p.s) == fingerprint_links(loaded));

    std::remove(path.c_str());
}


// ---------------------------------------------------------------------------
// Malformed input rejected.
// ---------------------------------------------------------------------------
TEST_CASE("schematic_from_json rejects malformed payloads",
          "[io][schematic-json][errors]") {
    Schematic out;
    REQUIRE_FALSE(schematic_from_json("", out));
    REQUIRE_FALSE(schematic_from_json("not json at all", out));
    REQUIRE_FALSE(schematic_from_json(R"({"cells":[{)", out));
}


// ---------------------------------------------------------------------------
// Unknown JSON keys are ignored (forward compat).
// ---------------------------------------------------------------------------
TEST_CASE("schematic_from_json ignores unknown top-level and per-cell keys",
          "[io][schematic-json][forward-compat]") {
    const std::string j = R"({
        "schema":"simall.workbench.schematic/v1",
        "metadata":{"author":"qa","when":"2026-05-19"},
        "cells":[
            {"id":1,"kind":"Geometry","label":"G","state":"UpToDate",
             "annotation":"future field",
             "ports":[{"id":0,"dir":"Output","name":"x","type":"x","required":false,"hint":"foo"}]}
        ],
        "links":[]
    })";
    Schematic out;
    REQUIRE(schematic_from_json(j, out));
    REQUIRE(out.size() == 1);
    REQUIRE(out.cells().front().label() == "G");
    REQUIRE(out.cells().front().state() == CellState::UpToDate);
}


// ---------------------------------------------------------------------------
// Link validation failures during load propagate.
// ---------------------------------------------------------------------------
TEST_CASE("schematic_from_json fails when a link references unknown ids",
          "[io][schematic-json][errors]") {
    const std::string j = R"({
        "cells":[
            {"id":1,"kind":"Geometry","label":"G","state":"UpToDate",
             "ports":[{"id":0,"dir":"Output","name":"x","type":"x","required":false}]}
        ],
        "links":[
            {"from_cell":1,"from_port":0,"to_cell":99,"to_port":0}
        ]
    })";
    Schematic out;
    REQUIRE_FALSE(schematic_from_json(j, out));
}
