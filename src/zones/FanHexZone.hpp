// =============================================================================
// SimAll Beta - Lumped Zone Models Subsystem
// File   : src/zones/FanHexZone.hpp
// Phase  : 22 — Lumped 0-D / 1-D zone models embedded in the 3-D mesh.
//
// Fan zone: applies a discrete pressure rise Δp(Q) across a cell zone
// described by a fan curve Δp = c0 + c1 Q + c2 Q². The implementation
// computes the volumetric flow Q through the fan's defining face zone,
// looks up Δp, and inserts it as a body force per unit volume oriented
// along a user-specified axis.
//
// Heat-exchanger zone: applies a heat source per cell q'''(T) = h_eff
// (T_amb - T) using a fitted heat-transfer-versus-velocity correlation
// h_eff(|u|) = a + b |u|^c, plus an optional Darcy-Forchheimer momentum
// resistance for the friction drop across the core.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <vector>

namespace simall::zones {

struct FanCurve {
    // Δp [Pa] = c0 + c1 Q + c2 Q²,  Q [m³/s]
    double c0 = 200.0;
    double c1 = 0.0;
    double c2 = -50.0;
};

struct FanZoneSpec {
    meshing::ZoneId cellZone   = 0;
    meshing::ZoneId throughFaceZone = 0;
    util::Vec3d     axis{0, 0, 1};   // unit vector of flow direction
    FanCurve        curve{};
};

struct HexZoneSpec {
    meshing::ZoneId cellZone = 0;
    double          T_ambient = 300.0;
    double          h_a = 100.0;     // W/(m³·K) — coefficient
    double          h_b = 800.0;
    double          h_c = 0.7;
    // Optional Darcy-Forchheimer (set to 0 to disable).
    double          darcy      = 0.0;   // 1/K_perm  [1/m²]
    double          forchheimer= 0.0;   // C_2       [1/m]
};

class FanZone {
public:
    void initialize(const meshing::Mesh& mesh, std::vector<FanZoneSpec> specs);
    /// Updates "S_FanMom" vector source per cell from current "U".
    void apply(solver::FieldRegistry& fields);

private:
    const meshing::Mesh* mesh_ = nullptr;
    std::vector<FanZoneSpec> specs_;
    /// Pre-computed cell-zone membership per cell.
    std::vector<std::int32_t> cellZoneId_;
    /// Boundary face IDs belonging to each fan's throughFaceZone, in order
    /// of `specs_`.
    std::vector<std::vector<meshing::FaceId>> throughFaces_;
};

class HexZone {
public:
    void initialize(const meshing::Mesh& mesh, std::vector<HexZoneSpec> specs);
    /// Writes per-cell "S_HexEn" (scalar, W/m³) and "S_HexMom" (vector, N/m³).
    void apply(solver::FieldRegistry& fields);

private:
    const meshing::Mesh* mesh_ = nullptr;
    std::vector<HexZoneSpec> specs_;
    std::vector<std::int32_t> cellZoneId_;
};

}  // namespace simall::zones
