// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/SolutionMonitors.cpp
// =============================================================================
#include "io/SolutionMonitors.hpp"

#include "core/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace simall::io
{

namespace
{

meshing::CellId nearest_cell(const meshing::Mesh& m, const util::Vec3d& p)
{
    const auto& C = m.cells();
    const std::size_t nC = C.size();
    meshing::CellId best = static_cast<meshing::CellId>(-1);
    double bestD = std::numeric_limits<double>::max();
    for (std::size_t c = 0; c < nC; ++c) {
        const double dx = C.centroidX[c] - p.x;
        const double dy = C.centroidY[c] - p.y;
        const double dz = C.centroidZ[c] - p.z;
        const double d = dx * dx + dy * dy + dz * dz;
        if (d < bestD) {
            bestD = d;
            best = c;
        }
    }
    return best;
}

} // namespace

// ----------------------------- PointProbe -----------------------------------

PointProbe::PointProbe(std::string n, std::string field, util::Vec3d loc)
    : name_(std::move(n)), field_(std::move(field)), loc_(loc)
{
}

MonitorSample PointProbe::sample(const meshing::Mesh& m,
                                 const solver::FieldRegistry& F,
                                 int it,
                                 double t)
{
    if (hostCell_ == static_cast<meshing::CellId>(-1))
        hostCell_ = nearest_cell(m, loc_);
    MonitorSample s;
    s.iteration = it;
    s.timeStamp = t;
    auto& Fnc = const_cast<solver::FieldRegistry&>(F);
    if (auto* sf = Fnc.find_scalar(field_); sf && hostCell_ < sf->size())
        s.value = (*sf)[hostCell_];
    else if (auto* vf = Fnc.find_vector(field_); vf && hostCell_ < vf->size()) {
        const double x = vf->x[hostCell_], y = vf->y[hostCell_], z = vf->z[hostCell_];
        s.value = std::sqrt(x * x + y * y + z * z);
    }
    return s;
}

// ----------------------------- SurfaceFlux ----------------------------------

SurfaceFlux::SurfaceFlux(std::string n, meshing::ZoneId z, double rho, std::string sf)
    : name_(std::move(n)), scalar_(std::move(sf)), zone_(z), rho_(rho)
{
}

MonitorSample SurfaceFlux::sample(const meshing::Mesh& m,
                                  const solver::FieldRegistry& F,
                                  int it,
                                  double t)
{
    MonitorSample s;
    s.iteration = it;
    s.timeStamp = t;
    auto& Fnc = const_cast<solver::FieldRegistry&>(F);
    const auto* U = Fnc.find_vector("U");
    if (!U)
        return s;
    const auto* phi = scalar_.empty() ? nullptr : Fnc.find_scalar(scalar_);
    const auto& Ff = m.faces();
    double sum = 0.0;
    for (std::size_t f = 0; f < Ff.size(); ++f) {
        if (Ff.boundaryZone[f] != zone_)
            continue;
        const meshing::CellId c = Ff.owner[f];
        const double Uxf = U->x[c], Uyf = U->y[c], Uzf = U->z[c];
        const double Fm = rho_ * (Uxf * Ff.areaX[f] + Uyf * Ff.areaY[f] + Uzf * Ff.areaZ[f]);
        sum += phi ? Fm * (*phi)[c] : Fm;
    }
    s.value = sum;
    return s;
}

// ----------------------------- VolumeStat -----------------------------------

VolumeStat::VolumeStat(std::string n, std::string field, meshing::ZoneId z, Kind k)
    : name_(std::move(n)), field_(std::move(field)), zone_(z), kind_(k)
{
}

MonitorSample VolumeStat::sample(const meshing::Mesh& m,
                                 const solver::FieldRegistry& F,
                                 int it,
                                 double t)
{
    MonitorSample s;
    s.iteration = it;
    s.timeStamp = t;
    auto& Fnc = const_cast<solver::FieldRegistry&>(F);
    auto* sf = Fnc.find_scalar(field_);
    if (!sf)
        return s;
    const auto& C = m.cells();
    double vol = 0.0, vw = 0.0;
    double vmin = std::numeric_limits<double>::infinity();
    double vmax = -std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < C.size(); ++c) {
        // NOTE: zone filtering not enforced here because per-cell zone tags
        // are an external attribute (see Porous/MRF/AMR for examples).
        // Pass a cellZone array via a future overload when domain-restricted
        // stats are required.
        (void) zone_;
        const double v = (*sf)[c];
        vol += C.volume[c];
        vw += v * C.volume[c];
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
    }
    switch (kind_) {
    case Kind::Mean:
        s.value = (vol > 0) ? vw / vol : 0.0;
        break;
    case Kind::Min:
        s.value = vmin;
        break;
    case Kind::Max:
        s.value = vmax;
        break;
    case Kind::IntegralVolume:
        s.value = vw;
        break;
    }
    return s;
}

// ----------------------------- MonitorManager -------------------------------

void MonitorManager::add(std::unique_ptr<IMonitor> m)
{
    monitors_.push_back(std::move(m));
    headerWritten_ = false;
}

void MonitorManager::open_csv(const std::string& path)
{
    csv_.open(path, std::ios::out | std::ios::trunc);
    if (!csv_)
        SIMALL_LOG_WARN("Monitors", "could not open '", path, "' for write");
    headerWritten_ = false;
}

void MonitorManager::advance(const meshing::Mesh& m,
                             const solver::FieldRegistry& F,
                             int it,
                             double t)
{
    if (monitors_.empty())
        return;
    if (csv_.is_open() && !headerWritten_) {
        csv_ << "iter,time";
        for (auto& mn : monitors_)
            csv_ << "," << mn->name();
        csv_ << "\n";
        headerWritten_ = true;
    }
    if (csv_.is_open())
        csv_ << it << "," << t;
    for (auto& mn : monitors_) {
        const auto s = mn->sample(m, F, it, t);
        if (csv_.is_open())
            csv_ << "," << s.value;
    }
    if (csv_.is_open()) {
        csv_ << "\n";
        csv_.flush();
    }
}

} // namespace simall::io
