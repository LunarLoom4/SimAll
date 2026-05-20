// =============================================================================
// SimAll Beta - Acoustics Subsystem
// File   : src/acoustics/FwhSurface.cpp
// =============================================================================
#include "acoustics/FwhSurface.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::acoustics
{

namespace
{

inline double dot3(const util::Vec3d& a, const util::Vec3d& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline util::Vec3d sub(const util::Vec3d& a, const util::Vec3d& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

/// Linear interpolation in the recorded time history (returns the
/// snapshot index `i` and the fraction `α` such that t lies between
/// times_[i] and times_[i+1]). Returns false if outside the recorded range.
struct Interp
{
    int i = 0;
    double a = 0.0;
    bool ok = false;
};
Interp locate_time(const std::vector<double>& T, double t)
{
    Interp r;
    if (T.size() < 2)
        return r;
    if (t <= T.front() || t >= T.back())
        return r;
    auto it = std::upper_bound(T.begin(), T.end(), t);
    const int i = static_cast<int>(it - T.begin()) - 1;
    r.i = i;
    r.a = (t - T[i]) / (T[i + 1] - T[i]);
    r.ok = true;
    return r;
}

} // namespace

void FwhSurface::initialize(const meshing::Mesh& m,
                            const std::vector<meshing::ZoneId>& zones,
                            FwhSurfaceProps p)
{
    mesh_ = &m;
    props_ = p;
    const auto& F = m.faces();
    panels_.clear();
    for (std::size_t f = 0; f < F.size(); ++f) {
        if (std::find(zones.begin(), zones.end(), F.boundaryZone[f]) == zones.end())
            continue;
        Panel pn;
        pn.centroid = {F.centroidX[f], F.centroidY[f], F.centroidZ[f]};
        const util::Vec3d A{F.areaX[f], F.areaY[f], F.areaZ[f]};
        pn.area = std::sqrt(dot3(A, A));
        if (pn.area < 1e-30)
            continue;
        pn.normal = {A.x / pn.area, A.y / pn.area, A.z / pn.area};
        panels_.push_back(pn);
    }
    times_.clear();
    snaps_.clear();
    SIMALL_LOG_INFO("FWH", "initialized: ", panels_.size(), " panels");
}

void FwhSurface::record_sample(double t, const solver::FieldRegistry& F, const util::Vec3d& Usurf)
{
    if (!mesh_)
        return;
    auto& Fnc = const_cast<solver::FieldRegistry&>(F);
    const auto* Uf = Fnc.find_vector("U");
    const auto* pf = Fnc.find_scalar("p");
    const auto* rhof = Fnc.find_scalar("rho");
    if (!Uf || !pf)
        return;

    Snapshot s;
    s.u.resize(panels_.size());
    s.p.resize(panels_.size());
    s.rho.resize(panels_.size());
    s.surfaceVel = Usurf;
    const auto& Mesh = *mesh_;
    const auto& Ff = Mesh.faces();
    // For each panel, look up the owner-cell velocity / pressure.
    // (Panel index is dense; rebuild via centroid match would be wasteful —
    // we cache the boundary face id once.)
    std::size_t pIdx = 0;
    for (std::size_t f = 0; f < Ff.size() && pIdx < panels_.size(); ++f) {
        if (Ff.areaX[f] * 0 != 0)
            continue;
        const util::Vec3d c{Ff.centroidX[f], Ff.centroidY[f], Ff.centroidZ[f]};
        if (std::abs(c.x - panels_[pIdx].centroid.x) > 1e-12
            || std::abs(c.y - panels_[pIdx].centroid.y) > 1e-12
            || std::abs(c.z - panels_[pIdx].centroid.z) > 1e-12)
            continue;
        const meshing::CellId oc = Ff.owner[f];
        s.u[pIdx] = {Uf->x[oc], Uf->y[oc], Uf->z[oc]};
        s.p[pIdx] = (*pf)[oc];
        s.rho[pIdx] = rhof ? (*rhof)[oc] : props_.rho0;
        ++pIdx;
    }
    times_.push_back(t);
    snaps_.push_back(std::move(s));
}

std::vector<double> FwhSurface::emit(const std::vector<FwhObserver>& obs,
                                     double t0,
                                     double t1,
                                     double dt) const
{
    const int Nt = std::max(1, static_cast<int>(std::round((t1 - t0) / dt)) + 1);
    std::vector<double> out(obs.size() * Nt, 0.0);
    if (snaps_.size() < 2)
        return out;
    const double rho0 = props_.rho0;
    const double a0 = props_.a0;
    const double inv4pi = 1.0 / (4.0 * M_PI);

    // Loop over (observer, observer time) and integrate over surface panels.
    for (std::size_t oi = 0; oi < obs.size(); ++oi) {
        for (int it = 0; it < Nt; ++it) {
            const double t = t0 + it * dt;
            double pT = 0.0, pL = 0.0;
            for (std::size_t pi = 0; pi < panels_.size(); ++pi) {
                const util::Vec3d& y = panels_[pi].centroid;
                const util::Vec3d r = sub(obs[oi].position, y);
                const double rmag = std::sqrt(dot3(r, r));
                if (rmag < 1e-12)
                    continue;
                const util::Vec3d rh{r.x / rmag, r.y / rmag, r.z / rmag};
                // Retarded time τ ≈ t - r/a0  (panel is stationary in this
                // permeable formulation; rigid-body motion would augment y(τ)).
                const double tau = t - rmag / a0;
                const auto in = locate_time(times_, tau);
                if (!in.ok)
                    continue;
                const Snapshot& A = snaps_[in.i];
                const Snapshot& B = snaps_[in.i + 1];
                const double a = in.a;
                const util::Vec3d u{(1 - a) * A.u[pi].x + a * B.u[pi].x,
                                    (1 - a) * A.u[pi].y + a * B.u[pi].y,
                                    (1 - a) * A.u[pi].z + a * B.u[pi].z};
                const double pP = (1 - a) * A.p[pi] + a * B.p[pi];
                const double rhoP = (1 - a) * A.rho[pi] + a * B.rho[pi];
                const util::Vec3d Us{(1 - a) * A.surfaceVel.x + a * B.surfaceVel.x,
                                     (1 - a) * A.surfaceVel.y + a * B.surfaceVel.y,
                                     (1 - a) * A.surfaceVel.z + a * B.surfaceVel.z};
                const util::Vec3d& n = panels_[pi].normal;
                const double area = panels_[pi].area;

                // Surface-relative mass flux through panel:
                //   ṁ/area = ρ (u - U) · n
                const double Un_s = dot3(Us, n);
                const double un = dot3(u, n);
                const double mdotA = rhoP * (un - Un_s);
                // Mach component along observer direction (panel at rest)
                const double Mr = 0.0; // permeable: surface stationary
                const double oneMmr = 1.0 - Mr;
                // Thickness term: ρ0 (U_n + relative flux) / (r(1-Mr)^2)
                if (props_.permeable) {
                    pT += area * (rho0 * (Un_s + mdotA / rhoP)) / (rmag * oneMmr * oneMmr);
                }
                // Loading vector L_i = (p δ_ij + ρ u_i (u_j - U_j)) n_j
                util::Vec3d L{pP * n.x + rhoP * u.x * (un - Un_s),
                              pP * n.y + rhoP * u.y * (un - Un_s),
                              pP * n.z + rhoP * u.z * (un - Un_s)};
                const double Lr = dot3(L, rh);
                pL += area * Lr / (rmag * rmag * oneMmr * oneMmr);
            }
            out[oi * Nt + it] = inv4pi * (pT + pL / a0);
        }
    }
    return out;
}

} // namespace simall::acoustics
