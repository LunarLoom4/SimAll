// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/bc/Bc.hpp
// Phase  : 16 — Boundary-condition framework.
//
// IBoundaryCondition is the common interface for every boundary condition in
// SimAll Beta.  A BC is bound to a single mesh zone (`meshing::ZoneId`) and a
// single field variable name (e.g. "U.x", "p", "T", "k", "omega", "Y_CH4").
//
// During matrix assembly the solver calls `apply(ctx)` once per (zone,
// variable) pair.  The BC modifies the diagonal / off-diagonal entries of
// the CSR matrix and the RHS vector for every face owned by the bound zone.
//
//   * Dirichlet  : a_pp += large, b_p += large·value           (penalty)
//   * Neumann    : b_p += flux·area
//   * Robin      : a_pp += h·area, b_p += h·area·T_inf
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/CSRMatrix.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <memory>
#include <string>

namespace simall::solver::bc
{

/// Identifier for the variable whose equation is currently being assembled.
/// String-based to avoid coupling the BC library to every physics module.
using VariableName = std::string;

/// Per-assemble context handed to a boundary condition.  All non-owning.
struct BcContext
{
    const meshing::Mesh* mesh = nullptr;
    FieldRegistry* fields = nullptr;
    CSRMatrix* matrix = nullptr;                 // may be null for explicit BCs
    util::aligned_vector<double>* rhs = nullptr; // may be null
    VariableName variable;
    double dt = 0.0;
    double timeNow = 0.0;
};

/// Canonical BC family — used by BcFactory & UI inspector.
enum class BcKind
{
    Wall,
    Inlet,
    Outlet,
    Symmetry,
    Axisymmetric,
    Periodic,
    PorousJump,
    Fan,
    Overset,
    Interface
};

class IBoundaryCondition
{
public:
    virtual ~IBoundaryCondition() = default;

    /// Mesh zone this BC is bound to.
    meshing::ZoneId zone() const noexcept { return zone_; }
    void setZone(meshing::ZoneId z) noexcept { zone_ = z; }

    /// Canonical kind. UI / serialization use this.
    virtual BcKind kind() const noexcept = 0;
    virtual const char* name() const noexcept = 0;

    /// Apply BC contribution for the current (mesh, variable) context.
    /// Returns the number of faces touched.
    virtual std::size_t apply(BcContext& ctx) = 0;

    /// Optional initialization (e.g. compute precomputed face indices).
    virtual void initialize(const meshing::Mesh& /*mesh*/) {}

    /// Deep-copy.
    virtual std::unique_ptr<IBoundaryCondition> clone() const = 0;

protected:
    /// Penalty applied to enforce Dirichlet conditions (Neumann boundary
    /// conditions skip the penalty).  Magnitudes around 1e30 are standard.
    static constexpr double kPenalty = 1.0e30;

    /// Walk the FaceStorage and invoke `fn(faceId)` for every face whose
    /// boundaryZone matches this BC's zone.  Convenience helper for
    /// concrete BCs.
    template <class Fn> std::size_t for_each_face(const meshing::Mesh& m, Fn&& fn) const
    {
        const auto& F = m.faces();
        const std::size_t n = F.size();
        std::size_t hits = 0;
        for (std::size_t i = 0; i < n; ++i) {
            if (F.boundaryZone[i] == zone_) {
                fn(static_cast<meshing::FaceId>(i));
                ++hits;
            }
        }
        return hits;
    }

    /// Locate (or skip) the diagonal entry of `row` in CSRMatrix and add δ.
    /// Diagonal is assumed to exist (true after solver assembly).
    static void addToDiagonal(CSRMatrix& A, std::size_t row, double delta)
    {
        const int r = static_cast<int>(row);
        for (int k = A.rowPtr[r]; k < A.rowPtr[r + 1]; ++k) {
            if (A.colIdx[k] == r) {
                A.values[k] += delta;
                return;
            }
        }
    }

    /// Penalty-Dirichlet on owner cell of face f. Safe to call repeatedly.
    static void applyDirichlet(CSRMatrix& A,
                               util::aligned_vector<double>& b,
                               meshing::CellId cell,
                               double value)
    {
        addToDiagonal(A, cell, kPenalty);
        b[cell] += kPenalty * value;
    }

    /// Neumann (flux) contribution.
    static void applyNeumann(util::aligned_vector<double>& b,
                             meshing::CellId cell,
                             double flux,
                             double area)
    {
        b[cell] += flux * area;
    }

    /// Robin contribution: ∂φ/∂n + h(φ - φ∞) = 0.
    static void applyRobin(CSRMatrix& A,
                           util::aligned_vector<double>& b,
                           meshing::CellId cell,
                           double h,
                           double area,
                           double phi_inf)
    {
        addToDiagonal(A, cell, h * area);
        b[cell] += h * area * phi_inf;
    }

    meshing::ZoneId zone_ = 0;
};

} // namespace simall::solver::bc
