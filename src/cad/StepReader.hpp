// =============================================================================
// SimAll Beta - CAD Subsystem
// File   : src/cad/StepReader.hpp
// Phase  : 4.1 (CAD readers)
//
// ISO 10303 (STEP) importer wrapping OpenCASCADE's STEPControl_Reader plus
// the STEPCAFControl_Reader product-structure pass. Returns one ShapeHandle
// per call; geometry is healed only when `options.heal == true`.
// =============================================================================
#pragma once

#include "cad/CadKernel.hpp"

#include <string>
#include <vector>

namespace simall::cad {

struct StepReadOptions {
    double lengthUnitToMeters = 0.0;  // 0 → honour STEP file's declared units
    bool   readNames           = true;
    bool   readColors          = true;
    bool   readLayers          = false;
    bool   heal                = false;
    int    precisionMode       = 0;    // 0 = average, 1 = least, 2 = greatest
};

struct StepReadReport {
    std::size_t rootCount     = 0;
    std::size_t transferred   = 0;
    std::size_t failedRoots   = 0;
    std::vector<std::string> warnings;
};

class StepReader {
public:
    ShapeHandle    read(const std::string& path, const StepReadOptions& opts = {});
    StepReadReport lastReport() const noexcept { return report_; }

private:
    StepReadReport report_;
};

}  // namespace simall::cad
