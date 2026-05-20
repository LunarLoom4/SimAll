// =============================================================================
// SimAll Beta - Multiphase Subsystem
// File   : src/multiphase/IsoAdvector.cpp
// =============================================================================
#include "multiphase/IsoAdvector.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::multiphase
{

void IsoAdvector::initialize(const meshing::Mesh& m, IsoAdvectorParams params)
{
    mesh_ = &m;
    p_ = params;
    SIMALL_LOG_INFO("Multiphase",
                    "IsoAdvector initialised: α_iso=",
                    p_.alpha_iso,
                    " subSteps=",
                    p_.subSteps,
                    " (",
                    m.cells().size(),
                    " cells)");
}

double IsoAdvector::face_flux_alpha(meshing::FaceId f,
                                    const solver::ScalarField& alpha,
                                    const PlicPatch& patch,
                                    double dt,
                                    double massFlux) const
{
    if (!mesh_)
        return 0.0;
    const auto& Fc = mesh_->faces();
    const auto o = Fc.owner[f];
    const auto n = Fc.neighbor[f];
    const double areaMag = std::sqrt(Fc.areaX[f] * Fc.areaX[f] + Fc.areaY[f] * Fc.areaY[f]
                                     + Fc.areaZ[f] * Fc.areaZ[f]);
    if (areaMag < 1e-30)
        return 0.0;

    // Upwind donor cell for α-flux.
    const auto donor = (massFlux >= 0.0) ? o : n;
    const double aUp = (donor != meshing::kBoundaryCell) ? alpha[donor] : 0.0;

    // Volumetric face flux Q_f = massFlux / ρ (here we treat massFlux as a
    // volumetric flux; calling code must pre-normalise).  Swept volume:
    // V_sw = |Q_f| · dt.
    const double Vsw = std::abs(massFlux) * dt;

    // If donor cell has a valid PLIC interface, integrate the fraction of
    // V_sw lying in the {n·(x-x_c) < δ} half-space.  Otherwise α_f ≈ α_up.
    double alphaFace = aUp;
    if (patch.valid && donor != meshing::kBoundaryCell) {
        // Approximate the swept volume's α as the cell's interface
        // distance to the face: project face centre onto the cell PLIC plane.
        const double cx = Fc.centroidX[f];
        const double cy = Fc.centroidY[f];
        const double cz = Fc.centroidZ[f];
        const auto& C = mesh_->cells();
        const double dx = cx - C.centroidX[donor];
        const double dy = cy - C.centroidY[donor];
        const double dz = cz - C.centroidZ[donor];
        const double s = patch.normal.x * dx + patch.normal.y * dy + patch.normal.z * dz;
        // Linear blend across interface plane within ±h (h = cube-root V):
        const double h = std::cbrt(std::max(C.volume[donor], 1e-30));
        const double xi = std::clamp((patch.offset - s) / std::max(h, 1e-30), -1.0, 1.0);
        alphaFace = std::clamp(0.5 + 0.5 * xi, 0.0, 1.0);
    }

    // α-volume crossing the face during this sub-step (signed by owner→neigh).
    const double sign = (massFlux >= 0.0) ? +1.0 : -1.0;
    return sign * Vsw * alphaFace;
}

double IsoAdvector::step(double dt, solver::FieldRegistry& F)
{
    if (!mesh_)
        return 0.0;
    auto* aField = F.find_scalar("alpha");
    auto* mfField = F.find_scalar("massFlux");
    if (!aField || !mfField)
        return 0.0;

    const auto& Cc = mesh_->cells();
    const auto& Fc = mesh_->faces();
    const std::size_t nC = Cc.size();
    const std::size_t nF = Fc.size();

    // Reconstruct PLIC patches for the current α field.
    plic_.reconstruct(*mesh_, *aField);
    const auto& patches = plic_.patches();

    const int nSub = std::max(1, p_.subSteps);
    const double dts = dt / nSub;

    double l1Excess = 0.0;
    for (int s = 0; s < nSub; ++s) {
        util::aligned_vector<double> dAlpha(nC, 0.0);
        // Loop faces, accumulate α-flux on owner (+) and neighbour (−).
        for (std::size_t f = 0; f < nF; ++f) {
            const auto o = Fc.owner[f];
            const auto n = Fc.neighbor[f];
            const auto donor = ((*mfField)[f] >= 0.0) ? o : n;
            const PlicPatch& patch = (donor != meshing::kBoundaryCell && donor < patches.size())
                                         ? patches[donor]
                                         : PlicPatch{};
            const double dAlphaFlux = face_flux_alpha(
                static_cast<meshing::FaceId>(f), *aField, patch, dts, (*mfField)[f]);
            if (o != meshing::kBoundaryCell)
                dAlpha[o] -= dAlphaFlux / std::max(Cc.volume[o], 1e-30);
            if (n != meshing::kBoundaryCell)
                dAlpha[n] += dAlphaFlux / std::max(Cc.volume[n], 1e-30);
        }
        // Cell update.
        for (std::size_t c = 0; c < nC; ++c) {
            const double aNew = (*aField)[c] + dAlpha[c];
            if (p_.boundedness) {
                if (aNew < 0.0)
                    l1Excess += -aNew * Cc.volume[c];
                if (aNew > 1.0)
                    l1Excess += (aNew - 1.0) * Cc.volume[c];
                (*aField)[c] = std::clamp(aNew, 0.0, 1.0);
            } else {
                (*aField)[c] = aNew;
            }
        }
    }
    return l1Excess;
}

} // namespace simall::multiphase
