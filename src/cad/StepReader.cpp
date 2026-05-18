// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/StepReader.cpp
// =============================================================================
#include "cad/StepReader.hpp"
#include "cad/ShapeHandleInternal.hpp"
#include "cad/Healing.hpp"
#include "core/Logger.hpp"

#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>

#include <filesystem>
#include <stdexcept>

namespace simall::cad {

ShapeHandle StepReader::read(const std::string& path, const StepReadOptions& opts) {
    report_ = {};

    if (!std::filesystem::exists(path))
        throw std::runtime_error("StepReader: file not found: " + path);

    STEPControl_Reader reader;

    // Precision mode: 0=average, 1=least, 2=greatest. Maps to OCC enum.
    Interface_Static::SetIVal("read.precision.mode", opts.precisionMode);
    if (opts.readNames)  Interface_Static::SetIVal("read.stepcaf.subshapes.name", 1);
    if (opts.readColors) Interface_Static::SetIVal("read.color", 1);
    if (opts.readLayers) Interface_Static::SetIVal("read.layer", 1);

    IFSelect_ReturnStatus status = reader.ReadFile(path.c_str());
    if (status != IFSelect_RetDone)
        throw std::runtime_error("StepReader: read failed for " + path +
                                 " (status=" + std::to_string(int(status)) + ")");

    Standard_Integer nRoots = reader.NbRootsForTransfer();
    report_.rootCount = std::size_t(nRoots);
    Standard_Integer transferred = reader.TransferRoots();
    report_.transferred = std::size_t(transferred);
    if (transferred == 0)
        throw std::runtime_error("StepReader: no roots transferred from " + path);
    report_.failedRoots = std::size_t(nRoots - transferred);

    TopoDS_Shape shape = reader.OneShape();
    if (shape.IsNull())
        throw std::runtime_error("StepReader: assembled shape is null");

    // Optional unit override — STEP files default to mm; OCC reads internal
    // units honoring the file header. If caller forces a unit, scale.
    if (opts.lengthUnitToMeters > 0.0 && opts.lengthUnitToMeters != 1.0) {
        gp_Trsf t;
        t.SetScale(gp_Pnt(0, 0, 0), opts.lengthUnitToMeters);
        BRepBuilderAPI_Transform xf(shape, t, /*copy=*/true);
        shape = xf.Shape();
    }

    auto h = makeHandle(std::move(shape));
    if (opts.heal) Healing{}.repair(h);

    SIMALL_LOG_INFO("CAD/STEP", "Loaded ", path, " — ", transferred, "/", nRoots,
                    " roots, ", h.topology().size(), " topo nodes");
    return h;
}

}  // namespace simall::cad
