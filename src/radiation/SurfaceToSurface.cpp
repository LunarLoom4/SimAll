// =============================================================================
// SimAll Beta - Radiation Subsystem
// File   : src/radiation/SurfaceToSurface.cpp
// =============================================================================
#include "radiation/SurfaceToSurface.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace simall::radiation
{

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kSigmaSB = 5.670374419e-8;
}

bool SurfaceToSurface::initialize(const meshing::Mesh& mesh, const S2SProps& props)
{
    mesh_ = &mesh;
    p_ = props;
    SIMALL_LOG_INFO(
        "Radiation", "S2S init: nRays=", p_.nRaysPerFace, " occlusion=", p_.checkOcclusion);
    return true;
}

const S2SFaceSpec* SurfaceToSurface::spec_for_zone(meshing::ZoneId z) const
{
    for (const auto& s : zones_)
        if (s.zone == z)
            return &s;
    return nullptr;
}

double SurfaceToSurface::face_flux(std::size_t f) const
{
    for (std::size_t k = 0; k < faceList_.size(); ++k)
        if (faceList_[k] == f)
            return q_[k];
    return 0.0;
}

double SurfaceToSurface::view_factor(std::size_t i, std::size_t j) const
{
    const std::size_t n = faceList_.size();
    if (i >= n || j >= n)
        return 0.0;
    return F_[i * n + j];
}

std::size_t SurfaceToSurface::build_view_factors()
{
    if (!mesh_)
        return 0;
    const auto& Ff = mesh_->faces();
    faceList_.clear();
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.neighbor[f] != meshing::kBoundaryCell)
            continue;
        if (spec_for_zone(Ff.boundaryZone[f]) == nullptr)
            continue;
        faceList_.push_back(f);
    }
    const std::size_t n = faceList_.size();
    F_.assign(n * n, 0.0);
    q_.assign(n, 0.0);

    std::mt19937_64 rng(p_.rngSeed);
    std::uniform_real_distribution<double> U(0.0, 1.0);

    // Stratified hemisphere sampling: cast Nray rays from face i centroid
    // along the inward (-area) direction; check which other face j the
    // ray reaches first via centroid-projection LOS.  Accumulate
    //   F_ij ≈ (Σ_hit cos θ_j_hit / π) / Nray   in the limit (analytic
    // hemisphere cos-weighted = 1/π · cos θ).
    std::vector<double> nx(n), ny(n), nz(n), Aarea(n), cx(n), cy(n), cz(n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t f = faceList_[k];
        double ax = Ff.areaX[f], ay = Ff.areaY[f], az = Ff.areaZ[f];
        const double am = std::sqrt(ax * ax + ay * ay + az * az);
        Aarea[k] = am;
        // Inward normal (toward fluid cell).
        nx[k] = -ax / am;
        ny[k] = -ay / am;
        nz[k] = -az / am;
        cx[k] = Ff.centroidX[f];
        cy[k] = Ff.centroidY[f];
        cz[k] = Ff.centroidZ[f];
    }

    for (std::size_t i = 0; i < n; ++i) {
        // Build local basis (u,v) around n_i.
        const double Nx = nx[i], Ny = ny[i], Nz = nz[i];
        double tx = (std::abs(Nz) < 0.9) ? 0.0 : 1.0;
        double ty = 0.0;
        double tz = (std::abs(Nz) < 0.9) ? 1.0 : 0.0;
        double ux = Ny * tz - Nz * ty;
        double uy = Nz * tx - Nx * tz;
        double uz = Nx * ty - Ny * tx;
        const double um = std::sqrt(ux * ux + uy * uy + uz * uz);
        if (um < 1e-30)
            continue;
        ux /= um;
        uy /= um;
        uz /= um;
        const double vx = Ny * uz - Nz * uy;
        const double vy = Nz * ux - Nx * uz;
        const double vz = Nx * uy - Ny * ux;

        for (std::size_t r = 0; r < p_.nRaysPerFace; ++r) {
            // Cosine-weighted hemisphere sample.
            const double r1 = U(rng), r2 = U(rng);
            const double sR = std::sqrt(r1);
            const double th = 2 * kPi * r2;
            const double lu = sR * std::cos(th);
            const double lv = sR * std::sin(th);
            const double lw = std::sqrt(std::max(1.0 - r1, 0.0));
            const double dx = lu * ux + lv * vx + lw * Nx;
            const double dy = lu * uy + lv * vy + lw * Ny;
            const double dz = lu * uz + lv * vz + lw * Nz;

            // Find first face j (j!=i) whose centroid plane the ray crosses
            // with positive t and which faces back to it (n_j · -d > 0).
            double tBest = std::numeric_limits<double>::infinity();
            std::size_t jBest = static_cast<std::size_t>(-1);
            for (std::size_t j = 0; j < n; ++j) {
                if (j == i)
                    continue;
                const double njdotd = nx[j] * dx + ny[j] * dy + nz[j] * dz;
                if (njdotd >= -1e-12)
                    continue; // away-facing
                const double num =
                    (cx[j] - cx[i]) * nx[j] + (cy[j] - cy[i]) * ny[j] + (cz[j] - cz[i]) * nz[j];
                const double t = num / njdotd;
                if (t > 1e-9 && t < tBest) {
                    tBest = t;
                    jBest = j;
                }
            }
            if (jBest != static_cast<std::size_t>(-1)) {
                // Each cos-weighted hit accumulates 1/Nray of F_ij.
                F_[i * n + jBest] += 1.0 / static_cast<double>(p_.nRaysPerFace);
            }
        }
        // Row normalisation: clamp Σ_j ≤ 1 (open enclosure: rest escapes).
        double sum = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            sum += F_[i * n + j];
        if (sum > 1.0)
            for (std::size_t j = 0; j < n; ++j)
                F_[i * n + j] /= sum;
    }
    return n;
}

double SurfaceToSurface::solve()
{
    if (!mesh_)
        return 0.0;
    const auto& Ff = mesh_->faces();
    const std::size_t n = faceList_.size();
    if (n == 0)
        return 0.0;

    // Build (I - (1-ε) F)·J = ε E_b and solve via Gauss-Seidel.
    std::vector<double> eps(n), Eb(n), J(n);
    for (std::size_t k = 0; k < n; ++k) {
        const auto* s = spec_for_zone(Ff.boundaryZone[faceList_[k]]);
        const double e = std::clamp(s ? s->emissivity : 1.0, 1e-3, 1.0);
        const double T = s ? s->temperature : 300.0;
        eps[k] = e;
        Eb[k] = kSigmaSB * std::pow(T, 4);
        J[k] = Eb[k];
    }
    for (int it = 0; it < 200; ++it) {
        double maxDelta = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            double acc = 0.0;
            for (std::size_t j = 0; j < n; ++j)
                acc += F_[i * n + j] * J[j];
            const double Jnew = eps[i] * Eb[i] + (1.0 - eps[i]) * acc;
            maxDelta = std::max(maxDelta, std::abs(Jnew - J[i]));
            J[i] = Jnew;
        }
        if (maxDelta < 1e-6 * (1.0 + Eb[0]))
            break;
    }
    double qmax = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        double acc = 0.0;
        for (std::size_t j = 0; j < n; ++j)
            acc += F_[i * n + j] * J[j];
        q_[i] = (eps[i] >= 0.999) ? (Eb[i] - acc) : eps[i] / (1.0 - eps[i]) * (Eb[i] - J[i]);
        qmax = std::max(qmax, std::abs(q_[i]));
    }
    return qmax;
}

} // namespace simall::radiation
