// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/ForcesReport.cpp
// =============================================================================
#include "io/ForcesReport.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>

namespace simall::io
{

namespace
{

inline double dot3(const util::Vec3d& a, const util::Vec3d& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline util::Vec3d cross3(const util::Vec3d& a, const util::Vec3d& b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline util::Vec3d normalize(util::Vec3d v)
{
    const double n = std::sqrt(dot3(v, v));
    return n > 1e-30 ? util::Vec3d{v.x / n, v.y / n, v.z / n} : v;
}

} // namespace

ForceResult ForcesReport::compute(const meshing::Mesh& m,
                                  const solver::FieldRegistry& F,
                                  const std::vector<meshing::ZoneId>& zones,
                                  const ForceReference& ref,
                                  const std::string& muEffName)
{
    ForceResult R;
    auto& Fnc = const_cast<solver::FieldRegistry&>(F);
    const auto* pField = Fnc.find_scalar("p");
    const auto* UF = Fnc.find_vector("U");
    const auto* muE = Fnc.find_scalar(muEffName);
    if (!muE)
        muE = Fnc.find_scalar("mu");
    if (!pField || !UF)
        return R;

    const auto& Ff = m.faces();
    const auto& C = m.cells();

    for (std::size_t f = 0; f < Ff.size(); ++f) {
        const meshing::ZoneId z = Ff.boundaryZone[f];
        if (std::find(zones.begin(), zones.end(), z) == zones.end())
            continue;
        const meshing::CellId c = Ff.owner[f];
        if (c == meshing::kBoundaryCell)
            continue;

        util::Vec3d Avec{Ff.areaX[f], Ff.areaY[f], Ff.areaZ[f]};
        const double Amag = std::sqrt(dot3(Avec, Avec));
        if (Amag < 1e-30)
            continue;
        const util::Vec3d n{Avec.x / Amag, Avec.y / Amag, Avec.z / Amag};
        R.area += Amag;

        // Pressure force = -p · n · A   (outward on body when wall is no-slip)
        const double pf = (*pField)[c];
        const util::Vec3d Fp{-pf * Avec.x, -pf * Avec.y, -pf * Avec.z};

        // Viscous (skin friction) — wall-tangent shear. Off-wall sample is
        // the owner cell-centre at distance d_n = (x_c - x_f) · n .
        util::Vec3d Fv{0, 0, 0};
        if (muE) {
            const double mu = (*muE)[c];
            const util::Vec3d xc{C.centroidX[c], C.centroidY[c], C.centroidZ[c]};
            const util::Vec3d xf{Ff.centroidX[f], Ff.centroidY[f], Ff.centroidZ[f]};
            const util::Vec3d d{xc.x - xf.x, xc.y - xf.y, xc.z - xf.z};
            const double dn = std::abs(dot3(d, n));
            // Tangential cell velocity (subtract normal projection)
            const util::Vec3d uC{UF->x[c], UF->y[c], UF->z[c]};
            const double uN = dot3(uC, n);
            const util::Vec3d uT{uC.x - uN * n.x, uC.y - uN * n.y, uC.z - uN * n.z};
            // Wall shear stress vector τ_w = μ · u_t / d_n  (no-slip wall ⇒
            // u_wall = 0); force = -τ_w · A (reaction on body)
            const double s = -mu * Amag / std::max(dn, 1e-30);
            Fv = {s * uT.x, s * uT.y, s * uT.z};
        }

        R.Fpressure.x += Fp.x;
        R.Fpressure.y += Fp.y;
        R.Fpressure.z += Fp.z;
        R.Fviscous.x += Fv.x;
        R.Fviscous.y += Fv.y;
        R.Fviscous.z += Fv.z;

        const util::Vec3d xf{Ff.centroidX[f], Ff.centroidY[f], Ff.centroidZ[f]};
        const util::Vec3d r{
            xf.x - ref.momentCenter.x, xf.y - ref.momentCenter.y, xf.z - ref.momentCenter.z};
        const util::Vec3d Mp = cross3(r, Fp);
        const util::Vec3d Mv = cross3(r, Fv);
        R.Mpressure.x += Mp.x;
        R.Mpressure.y += Mp.y;
        R.Mpressure.z += Mp.z;
        R.Mviscous.x += Mv.x;
        R.Mviscous.y += Mv.y;
        R.Mviscous.z += Mv.z;
    }

    const util::Vec3d Ftot{
        R.Fpressure.x + R.Fviscous.x, R.Fpressure.y + R.Fviscous.y, R.Fpressure.z + R.Fviscous.z};
    const util::Vec3d Mtot{
        R.Mpressure.x + R.Mviscous.x, R.Mpressure.y + R.Mviscous.y, R.Mpressure.z + R.Mviscous.z};
    const util::Vec3d lh = normalize(ref.liftAxis);
    const util::Vec3d dh = normalize(ref.dragAxis);
    const util::Vec3d sh = normalize(ref.sideAxis);
    const double q = 0.5 * ref.rhoRef * ref.velRef * ref.velRef * ref.areaRef;
    const double qL = q * ref.lengthRef;
    if (q > 1e-30) {
        R.CL = dot3(Ftot, lh) / q;
        R.CD = dot3(Ftot, dh) / q;
        R.CS = dot3(Ftot, sh) / q;
    }
    if (qL > 1e-30) {
        R.CMl = dot3(Mtot, lh) / qL;
        R.CMd = dot3(Mtot, dh) / qL;
        R.CMs = dot3(Mtot, sh) / qL;
    }
    return R;
}

} // namespace simall::io
