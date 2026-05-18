// =============================================================================
// SimAll Beta - Meshing Subsystem
// File   : src/meshing/MeshQuality.hpp
// Phase  : 4.6 — Mesh quality metrics for unstructured polyhedral meshes.
//
//   - skewness            : (centroid-line offset from face centroid) / R
//   - non-orthogonality   : angle (deg) between face normal and the
//                            cell-to-cell vector (interior faces only)
//   - aspect ratio        : per-cell ratio of max face area to min face area
//   - volume positivity   : flag for negative/zero volume cells
//
// Outputs per-face and per-cell statistics, plus a histogram report
// suitable for export to the GUI or CSV.
// =============================================================================
#pragma once

#include "meshing/MeshStorage.hpp"
#include "utilities/AlignedAllocator.hpp"

#include <array>
#include <string>
#include <vector>

namespace simall::meshing {

struct QualityHistogram {
    std::array<std::size_t, 10> bins{};
    double                      vmin = 0.0;
    double                      vmax = 0.0;
};

struct MeshQualityReport {
    // Per-face / per-cell scalar arrays (size = nFaces / nCells).
    util::aligned_vector<double> skewness;
    util::aligned_vector<double> nonOrthoDeg;
    util::aligned_vector<double> aspectRatio;
    std::vector<std::uint8_t>    negativeVolume;

    QualityHistogram histSkewness;
    QualityHistogram histNonOrtho;
    QualityHistogram histAspect;

    double maxSkewness  = 0.0;
    double maxNonOrtho  = 0.0;
    double maxAspect    = 0.0;
    std::size_t negativeCount = 0;
    std::size_t nFaces = 0, nCells = 0;
};

class MeshQuality {
public:
    static MeshQualityReport evaluate(const Mesh& mesh);
    static std::string       format(const MeshQualityReport& r);
};

}  // namespace simall::meshing
