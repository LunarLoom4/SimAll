// =============================================================================
// SimAll Beta - Solver Subsystem
// File   : src/solver/RhieChowInterpolation.cpp
// =============================================================================
#include "solver/RhieChowInterpolation.hpp"

#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::solver
{

namespace
{
const BoundarySpec* find_bc(const std::vector<BoundarySpec>& bcs, meshing::ZoneId z)
{
    for (const auto& b : bcs)
        if (b.zone == z)
            return &b;
    return nullptr;
}

inline double harmonic(double a, double b)
{
    const double s = a + b;
    return (s > 1.0e-30) ? (2.0 * a * b / s) : 0.0;
}
} // namespace

void RhieChowInterpolation::compute_face_mass_flux(const VectorField& U,
                                                   const ScalarField& p,
                                                   const util::aligned_vector<double>& aP,
                                                   util::aligned_vector<double>& outFlux,
                                                   const VectorField* gradPIn) const
{
    const auto& F = mesh_.faces();
    const auto& C = mesh_.cells();
    const std::size_t nF = F.size();
    const std::size_t nC = C.size();
    outFlux.assign(nF, 0.0);

    // Cell-centred ∇p (least-squares) — required for the second piece of
    // the Rhie-Chow correction.
    VectorField gradPLocal;
    const VectorField* gradP = gradPIn;
    if (!gradP) {
        LeastSquaresGradient G(mesh_);
        G.evaluate(p, gradPLocal);
        gradP = &gradPLocal;
    }

    for (std::size_t f = 0; f < nF; ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];

        if (n != meshing::kBoundaryCell) {
            // Average (linear) face velocity
            const double Uxf = 0.5 * (U.x[o] + U.x[n]);
            const double Uyf = 0.5 * (U.y[o] + U.y[n]);
            const double Uzf = 0.5 * (U.z[o] + U.z[n]);

            // d_ON = (centroid_n - centroid_o)
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double d2 = dx * dx + dy * dy + dz * dz;
            const double dlen = std::sqrt(std::max(d2, 1.0e-60));

            // Cell-volume ratio (V/aP)_f — harmonic average suppresses
            // contamination from poorly-conditioned cells.
            const double VaP_o = C.volume[o] / std::max(aP[o], 1.0e-30);
            const double VaP_n = C.volume[n] / std::max(aP[n], 1.0e-30);
            const double VaPf = harmonic(VaP_o, VaP_n);

            // Compact pressure gradient along d̂_ON
            const double dpdL_compact = (p[n] - p[o]) / dlen;
            // Mean of cell-centred gradients projected onto d̂_ON
            const double meanGx = 0.5 * (gradP->x[o] + gradP->x[n]);
            const double meanGy = 0.5 * (gradP->y[o] + gradP->y[n]);
            const double meanGz = 0.5 * (gradP->z[o] + gradP->z[n]);
            const double dpdL_mean = (meanGx * dx + meanGy * dy + meanGz * dz) / dlen;

            // Rhie-Chow correction projected along d̂_ON:
            // U_f_RC = U_f_avg - (V/aP)_f * (dpdL_compact - dpdL_mean) * d̂_ON
            const double corr = VaPf * (dpdL_compact - dpdL_mean) / dlen;
            const double UxRC = Uxf - corr * dx;
            const double UyRC = Uyf - corr * dy;
            const double UzRC = Uzf - corr * dz;

            outFlux[f] = rho_ * (UxRC * Ax + UyRC * Ay + UzRC * Az);
        } else {
            // Boundary: use BC-specific value (no Rhie-Chow on boundary;
            // pressure-velocity decoupling is interior phenomenon).
            const auto* bc = find_bc(bcs_, F.boundaryZone[f]);
            double Uxf, Uyf, Uzf;
            if (bc && (bc->type == BCType::VelocityInlet || bc->type == BCType::MovingWall)) {
                Uxf = bc->vectorValue[0];
                Uyf = bc->vectorValue[1];
                Uzf = bc->vectorValue[2];
            } else if (bc
                       && (bc->type == BCType::NoSlipWall || bc->type == BCType::Wall
                           || bc->type == BCType::Symmetry)) {
                Uxf = Uyf = Uzf = 0.0;
            } else {
                Uxf = U.x[o];
                Uyf = U.y[o];
                Uzf = U.z[o];
            }
            outFlux[f] = rho_ * (Uxf * Ax + Uyf * Ay + Uzf * Az);
        }
    }
    (void) nC;
}

void RhieChowInterpolation::compute_face_mass_flux_var(const ScalarField& rho,
                                                       const VectorField& U,
                                                       const ScalarField& p,
                                                       const util::aligned_vector<double>& aP,
                                                       util::aligned_vector<double>& outFlux,
                                                       const VectorField* gradPIn) const
{
    const auto& F = mesh_.faces();
    const auto& C = mesh_.cells();
    outFlux.assign(F.size(), 0.0);

    VectorField gradPLocal;
    const VectorField* gradP = gradPIn;
    if (!gradP) {
        LeastSquaresGradient G(mesh_);
        G.evaluate(p, gradPLocal);
        gradP = &gradPLocal;
    }

    for (std::size_t f = 0; f < F.size(); ++f) {
        const meshing::CellId o = F.owner[f];
        const meshing::CellId n = F.neighbor[f];
        const double Ax = F.areaX[f], Ay = F.areaY[f], Az = F.areaZ[f];
        if (n != meshing::kBoundaryCell) {
            const double Uxf = 0.5 * (U.x[o] + U.x[n]);
            const double Uyf = 0.5 * (U.y[o] + U.y[n]);
            const double Uzf = 0.5 * (U.z[o] + U.z[n]);
            const double rhoF = 0.5 * (rho[o] + rho[n]);
            const double dx = C.centroidX[n] - C.centroidX[o];
            const double dy = C.centroidY[n] - C.centroidY[o];
            const double dz = C.centroidZ[n] - C.centroidZ[o];
            const double dlen = std::sqrt(std::max(dx * dx + dy * dy + dz * dz, 1.0e-60));
            const double VaP_o = C.volume[o] / std::max(aP[o], 1.0e-30);
            const double VaP_n = C.volume[n] / std::max(aP[n], 1.0e-30);
            const double VaPf = harmonic(VaP_o, VaP_n);
            const double dpdL_compact = (p[n] - p[o]) / dlen;
            const double mGx = 0.5 * (gradP->x[o] + gradP->x[n]);
            const double mGy = 0.5 * (gradP->y[o] + gradP->y[n]);
            const double mGz = 0.5 * (gradP->z[o] + gradP->z[n]);
            const double dpdL_mean = (mGx * dx + mGy * dy + mGz * dz) / dlen;
            const double corr = VaPf * (dpdL_compact - dpdL_mean) / dlen;
            outFlux[f] =
                rhoF * ((Uxf - corr * dx) * Ax + (Uyf - corr * dy) * Ay + (Uzf - corr * dz) * Az);
        } else {
            const auto* bc = find_bc(bcs_, F.boundaryZone[f]);
            double Uxf = U.x[o], Uyf = U.y[o], Uzf = U.z[o];
            if (bc && (bc->type == BCType::VelocityInlet || bc->type == BCType::MovingWall)) {
                Uxf = bc->vectorValue[0];
                Uyf = bc->vectorValue[1];
                Uzf = bc->vectorValue[2];
            } else if (bc
                       && (bc->type == BCType::NoSlipWall || bc->type == BCType::Wall
                           || bc->type == BCType::Symmetry)) {
                Uxf = Uyf = Uzf = 0.0;
            }
            outFlux[f] = rho[o] * (Uxf * Ax + Uyf * Ay + Uzf * Az);
        }
    }
}

} // namespace simall::solver
