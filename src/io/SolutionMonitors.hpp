// =============================================================================
// SimAll Beta - I/O Subsystem
// File   : src/io/SolutionMonitors.hpp
// Phase  : 17 — Runtime monitors / probes. Three kinds:
//
//   1. PointProbe   — Cell-centred trilinear sample at a fixed (x,y,z).
//   2. SurfaceFlux  — ∫ (ρ U · n) φ dA over a boundary zone (e.g. mass
//                      flow, heat flux, species mass flow).
//   3. VolumeStat   — Volume-averaged value, min, max over a cell zone.
//
// Each monitor writes a CSV time history to disk and emits an event over
// the EventBus so the GUI can plot it live. Monitors are advanced once
// per outer iteration of the SIMPLE / coupled solver.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "solver/FieldRegistry.hpp"
#include "utilities/MathTypes.hpp"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace simall::io
{

struct MonitorSample
{
    int iteration = 0;
    double timeStamp = 0.0;
    double value = 0.0;
};

class IMonitor
{
public:
    virtual ~IMonitor() = default;
    virtual const std::string& name() const noexcept = 0;
    virtual MonitorSample sample(const meshing::Mesh& mesh,
                                 const solver::FieldRegistry& fields,
                                 int iteration,
                                 double t) = 0;
};

/// Cell-centred sample using nearest-cell lookup (could be upgraded to
/// trilinear inside the host cell if richer interpolation is needed).
class PointProbe : public IMonitor
{
public:
    PointProbe(std::string name, std::string fieldName, util::Vec3d location);
    const std::string& name() const noexcept override { return name_; }
    MonitorSample sample(const meshing::Mesh& m,
                         const solver::FieldRegistry& F,
                         int iter,
                         double t) override;

private:
    std::string name_, field_;
    util::Vec3d loc_;
    meshing::CellId hostCell_ = static_cast<meshing::CellId>(-1);
};

/// ∫ (ρ U · n) φ dA over a boundary zone (φ defaults to 1 → mass flow).
class SurfaceFlux : public IMonitor
{
public:
    SurfaceFlux(std::string name, meshing::ZoneId zone, double rho, std::string scalarField = {});
    const std::string& name() const noexcept override { return name_; }
    MonitorSample sample(const meshing::Mesh& m,
                         const solver::FieldRegistry& F,
                         int iter,
                         double t) override;

private:
    std::string name_, scalar_;
    meshing::ZoneId zone_;
    double rho_;
};

/// Volume-averaged value over a cell zone (zone == 0 → entire domain).
class VolumeStat : public IMonitor
{
public:
    enum class Kind
    {
        Mean,
        Min,
        Max,
        IntegralVolume
    };
    VolumeStat(std::string name, std::string fieldName, meshing::ZoneId zone, Kind kind);
    const std::string& name() const noexcept override { return name_; }
    MonitorSample sample(const meshing::Mesh& m,
                         const solver::FieldRegistry& F,
                         int iter,
                         double t) override;

private:
    std::string name_, field_;
    meshing::ZoneId zone_;
    Kind kind_;
};

class MonitorManager
{
public:
    void add(std::unique_ptr<IMonitor> m);
    void open_csv(const std::string& path);
    void advance(const meshing::Mesh& mesh,
                 const solver::FieldRegistry& fields,
                 int iteration,
                 double t);
    const std::vector<std::unique_ptr<IMonitor>>& monitors() const noexcept { return monitors_; }

private:
    std::vector<std::unique_ptr<IMonitor>> monitors_;
    std::ofstream csv_;
    bool headerWritten_ = false;
};

} // namespace simall::io
