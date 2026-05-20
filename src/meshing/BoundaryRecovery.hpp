// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/BoundaryRecovery.hpp
// Phase  : 6.11 — Boundary edge / face recovery for constrained Delaunay
//          tetrahedrisation (Si 2010 / Shewchuk 1998).
//
// After an unconstrained 3-D Delaunay triangulation, some boundary edges
// and faces may not appear in the tetrahedral mesh (the so-called
// "missing-segment" or "missing-face" problem).  Boundary recovery
// enforces the input PSLG by either:
//
//   * Edge swapping (2D analogue: flip diagonals)
//   * Steiner-point insertion: midpoint subdivision of the missing
//     segment / face, repeated until recovery succeeds
//
// API operates on an existing tetra mesh + a list of required edges (pairs
// of node IDs) and required faces (triples).  Returns the number of
// missing entities recovered and the number of Steiner points added.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"

#include <cstdint>
#include <utility>
#include <vector>

namespace simall::meshing
{

struct RequiredEdge
{
    NodeId a, b;
};
struct RequiredFace
{
    NodeId a, b, c;
};

struct BoundaryRecoveryProps
{
    std::size_t maxSteinerInsertions = 1024;
    double coordTol = 1e-9;
};

struct RecoveryReport
{
    std::size_t edgesMissing = 0;
    std::size_t edgesRecovered = 0;
    std::size_t facesMissing = 0;
    std::size_t facesRecovered = 0;
    std::size_t steinerAdded = 0;
};

class BoundaryRecovery
{
public:
    void initialize(BoundaryRecoveryProps props);

    void add_required_edge(RequiredEdge e) { edges_.push_back(e); }
    void add_required_face(RequiredFace f) { faces_.push_back(f); }

    RecoveryReport recover(Mesh& mesh);

    const BoundaryRecoveryProps& props() const noexcept { return p_; }

private:
    bool edge_present(const Mesh& mesh, const RequiredEdge& e) const;
    bool face_present(const Mesh& mesh, const RequiredFace& f) const;
    NodeId insert_steiner(Mesh& mesh, const RequiredEdge& e);

    BoundaryRecoveryProps p_{};
    std::vector<RequiredEdge> edges_;
    std::vector<RequiredFace> faces_;
};

} // namespace simall::meshing
