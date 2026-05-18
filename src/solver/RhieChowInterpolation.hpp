// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/RhieChowInterpolation.hpp
// Phase  : 6.5 — Rhie & Chow (1983) face-velocity interpolation
//
// Collocated finite-volume schemes suffer from pressure-velocity decoupling
// (the "checkerboard" mode) because a centred face velocity:
//
//      U_f = ½(U_O + U_N)
//
// is insensitive to a saw-tooth pressure field whose gradient at the face is
// in fact non-zero.  Rhie & Chow's remedy adds a pressure-stencil correction:
//
//      U_f = ½(U_O + U_N)
//          - (V/aP)_f · [ (p_N - p_O)/|d_ON|·(d̂_ON)  -  ½(∇p_O + ∇p_N) ]
//
// where (V/aP)_f is the harmonic-averaged momentum-diagonal volume coefficient
// and ∇p_O, ∇p_N are the cell-centred pressure gradients (least-squares).
//
// References:
//   Rhie & Chow, AIAA J 21 (11), 1525-1532 (1983).
//   Ferziger & Perić, Computational Methods for Fluid Dynamics, §7.5.3.
//   Jasak, PhD thesis, Imperial College (1996), §3.10.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "FieldRegistry.hpp"
#include "solver/Solver.hpp"

namespace simall::solver {

class RhieChowInterpolation {
public:
    /// rho is constant density (incompressible). For compressible/variable-
    /// density flows pass a face-averaged value through compute_face_mass_flux_var.
    RhieChowInterpolation(const meshing::Mesh& mesh,
                          const std::vector<BoundarySpec>& bcs,
                          double rho)
        : mesh_(mesh), bcs_(bcs), rho_(rho) {}

    /// Compute Rhie-Chow stabilised mass flux at every face:
    ///   φ_f = ρ · (A · U_f^RC)
    /// where U_f^RC is the corrected face velocity. Writes into outFlux
    /// (resized to faces().size()).
    ///
    /// Inputs:
    ///   U      cell-centred velocity field (existing in registry)
    ///   p      cell-centred pressure field
    ///   aP     per-cell momentum diagonal coefficient
    ///   gradP  optional cell-centred ∇p (if null, computed internally)
    void compute_face_mass_flux(const VectorField& U,
                                const ScalarField& p,
                                const util::aligned_vector<double>& aP,
                                util::aligned_vector<double>& outFlux,
                                const VectorField* gradP = nullptr) const;

    /// Variable-density variant: ρ supplied per cell. Used by VOF and
    /// compressible-incompressible bridges.
    void compute_face_mass_flux_var(const ScalarField& rho,
                                    const VectorField& U,
                                    const ScalarField& p,
                                    const util::aligned_vector<double>& aP,
                                    util::aligned_vector<double>& outFlux,
                                    const VectorField* gradP = nullptr) const;

private:
    const meshing::Mesh&              mesh_;
    const std::vector<BoundarySpec>&  bcs_;
    double                            rho_;
};

}  // namespace simall::solver
