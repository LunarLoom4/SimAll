// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/MonteCarloRte.cpp
// =============================================================================
#include "radiation/MonteCarloRte.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::radiation
{

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kSigmaSB = 5.670374419e-8;

struct Vec3
{
    double x, y, z;
};

inline double dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline double mag(const Vec3& a)
{
    return std::sqrt(dot(a, a));
}
} // namespace

bool MonteCarloRte::initialize(const meshing::Mesh& mesh,
                               solver::FieldRegistry& F,
                               const MonteCarloProps& props)
{
    mesh_ = &mesh;
    F_ = &F;
    p_ = props;
    rng_.seed(p_.rngSeed);
    const std::size_t nC = mesh.cells().size();
    F.scalar("T", nC);
    F.scalar("S_rad", nC);
    q_wall_.assign(mesh.faces().size(), 0.0);
    SIMALL_LOG_INFO("Radiation",
                    "MonteCarloRte init: N_p/cell=",
                    p_.nPhotonsPerCell,
                    " κ=",
                    p_.absorption,
                    " σ_s=",
                    p_.scattering);
    return true;
}

const McWallSpec* MonteCarloRte::find_wall(meshing::ZoneId z) const
{
    for (const auto& w : walls_)
        if (w.zone == z)
            return &w;
    return nullptr;
}

double MonteCarloRte::wall_face_flux(std::size_t f) const
{
    return (f < q_wall_.size()) ? q_wall_[f] : 0.0;
}

double MonteCarloRte::run()
{
    if (!mesh_ || !F_)
        return 0.0;
    const auto& C = mesh_->cells();
    const auto& Ff = mesh_->faces();
    const std::size_t nC = C.size();
    const auto* T = F_->find_scalar("T");
    auto* S = F_->find_scalar("S_rad");
    if (!T || !S)
        return 0.0;
    std::fill(q_wall_.begin(), q_wall_.end(), 0.0);

    const double beta = std::max(p_.absorption + p_.scattering, 1e-30);
    const double albedo = p_.scattering / beta;

    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::vector<double> Eabs(nC, 0.0);
    double Etotal = 0.0;

    for (std::size_t c = 0; c < nC; ++c) {
        const double Tc = std::max((*T)[c], 1.0);
        const double E_cell = 4.0 * p_.absorption * kSigmaSB * std::pow(Tc, 4) * C.volume[c];
        if (E_cell <= 0.0 || p_.nPhotonsPerCell == 0)
            continue;
        const double Ep = E_cell / static_cast<double>(p_.nPhotonsPerCell);
        Etotal += E_cell;

        for (std::size_t ip = 0; ip < p_.nPhotonsPerCell; ++ip) {
            // Isotropic direction.
            const double xi1 = U(rng_), xi2 = U(rng_);
            const double cosTh = 1.0 - 2.0 * xi1;
            const double sinTh = std::sqrt(std::max(1.0 - cosTh * cosTh, 0.0));
            const double phi = 2.0 * kPi * xi2;
            Vec3 dir{sinTh * std::cos(phi), sinTh * std::sin(phi), cosTh};
            Vec3 pos{C.centroidX[c], C.centroidY[c], C.centroidZ[c]};
            meshing::CellId cur = static_cast<meshing::CellId>(c);
            double E = Ep;

            for (std::size_t seg = 0; seg < p_.maxSegments && E > 1e-20; ++seg) {
                // Sample free path.
                const double u = U(rng_);
                const double Lcoll = -std::log(std::max(u, 1e-300)) / beta;

                // Find next-face exit distance from cell `cur`.
                const std::size_t off = C.faceOffsets[cur];
                const std::size_t end = C.faceOffsets[cur + 1];
                double tExit = std::numeric_limits<double>::infinity();
                std::size_t fHit = static_cast<std::size_t>(-1);
                for (std::size_t k = off; k < end; ++k) {
                    const std::size_t fi = C.faceIndices[k];
                    Vec3 nrm{Ff.areaX[fi], Ff.areaY[fi], Ff.areaZ[fi]};
                    const double am = mag(nrm);
                    if (am < 1e-30)
                        continue;
                    nrm.x /= am;
                    nrm.y /= am;
                    nrm.z /= am;
                    // Owner-side outward normal; flip if `cur` is the neighbour.
                    if (Ff.neighbor[fi] != meshing::kBoundaryCell && Ff.owner[fi] != cur) {
                        nrm.x = -nrm.x;
                        nrm.y = -nrm.y;
                        nrm.z = -nrm.z;
                    }
                    const double sn = dot(dir, nrm);
                    if (sn <= 1e-30)
                        continue;
                    Vec3 fc{Ff.centroidX[fi], Ff.centroidY[fi], Ff.centroidZ[fi]};
                    const double t =
                        ((fc.x - pos.x) * nrm.x + (fc.y - pos.y) * nrm.y + (fc.z - pos.z) * nrm.z)
                        / sn;
                    if (t > 1e-12 && t < tExit) {
                        tExit = t;
                        fHit = fi;
                    }
                }
                if (fHit == static_cast<std::size_t>(-1))
                    break;

                if (Lcoll < tExit) {
                    // Volumetric event at distance Lcoll inside `cur`.
                    pos.x += Lcoll * dir.x;
                    pos.y += Lcoll * dir.y;
                    pos.z += Lcoll * dir.z;
                    if (U(rng_) < albedo) {
                        // Scatter (isotropic).
                        const double a1 = U(rng_), a2 = U(rng_);
                        const double cs = 1 - 2 * a1;
                        const double ss = std::sqrt(std::max(1 - cs * cs, 0.0));
                        const double ph = 2 * kPi * a2;
                        dir = {ss * std::cos(ph), ss * std::sin(ph), cs};
                        continue;
                    } else {
                        Eabs[cur] += E;
                        E = 0.0;
                        break;
                    }
                }

                // Step to the face and either cross or terminate at wall.
                pos.x += tExit * dir.x;
                pos.y += tExit * dir.y;
                pos.z += tExit * dir.z;
                if (Ff.neighbor[fHit] == meshing::kBoundaryCell) {
                    const auto z = Ff.boundaryZone[fHit];
                    const auto* wp = find_wall(z);
                    if (!wp) { // open / black sink
                        q_wall_[fHit] += E;
                        E = 0.0;
                        break;
                    }
                    const double eps = std::clamp(wp->emissivity, 0.0, 1.0);
                    q_wall_[fHit] += eps * E;
                    E *= (1.0 - eps);
                    if (E < 1e-20)
                        break;
                    // Diffuse reflection: cosine-weighted hemisphere about
                    // the inward face normal.
                    Vec3 nrm{-Ff.areaX[fHit], -Ff.areaY[fHit], -Ff.areaZ[fHit]};
                    const double am = mag(nrm);
                    nrm.x /= am;
                    nrm.y /= am;
                    nrm.z /= am;
                    const double r1 = U(rng_), r2 = U(rng_);
                    const double sR = std::sqrt(r1);
                    const double th = 2 * kPi * r2;
                    Vec3 ref{
                        sR * std::cos(th), sR * std::sin(th), std::sqrt(std::max(1.0 - r1, 0.0))};
                    // Build local basis around nrm.
                    Vec3 t0 = (std::abs(nrm.z) < 0.9) ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
                    Vec3 u_axis{nrm.y * t0.z - nrm.z * t0.y,
                                nrm.z * t0.x - nrm.x * t0.z,
                                nrm.x * t0.y - nrm.y * t0.x};
                    const double um = mag(u_axis);
                    if (um < 1e-30)
                        break;
                    u_axis.x /= um;
                    u_axis.y /= um;
                    u_axis.z /= um;
                    Vec3 v_axis{nrm.y * u_axis.z - nrm.z * u_axis.y,
                                nrm.z * u_axis.x - nrm.x * u_axis.z,
                                nrm.x * u_axis.y - nrm.y * u_axis.x};
                    dir = {ref.x * u_axis.x + ref.y * v_axis.x + ref.z * nrm.x,
                           ref.x * u_axis.y + ref.y * v_axis.y + ref.z * nrm.y,
                           ref.x * u_axis.z + ref.y * v_axis.z + ref.z * nrm.z};
                    continue;
                }
                // Interior face: move to the neighbour cell.
                cur = (Ff.owner[fHit] == cur) ? Ff.neighbor[fHit] : Ff.owner[fHit];
            }
        }
    }

    // Convert absorption per cell to S_rad source (W/m³).
    for (std::size_t c = 0; c < nC; ++c) {
        const double Tc = std::max((*T)[c], 1.0);
        const double Eemit = 4.0 * p_.absorption * kSigmaSB * std::pow(Tc, 4) * C.volume[c];
        const double V = std::max(C.volume[c], 1e-30);
        (*S)[c] += (Eabs[c] - Eemit) / V;
    }
    // Convert wall energies to per-area flux.
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.neighbor[f] != meshing::kBoundaryCell)
            continue;
        const double aMag = std::sqrt(Ff.areaX[f] * Ff.areaX[f] + Ff.areaY[f] * Ff.areaY[f]
                                      + Ff.areaZ[f] * Ff.areaZ[f]);
        if (aMag > 1e-30)
            q_wall_[f] /= aMag;
    }
    return Etotal;
}

} // namespace simall::radiation
