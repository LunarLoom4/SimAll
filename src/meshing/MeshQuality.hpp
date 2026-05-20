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
#include <string_view>
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
    QualityHistogram histVolume;

    double maxSkewness  = 0.0;
    double maxNonOrtho  = 0.0;
    double maxAspect    = 0.0;
    double minVolume    = 0.0;
    double maxVolume    = 0.0;
    double meanVolume   = 0.0;
    double totalVolume  = 0.0;
    std::size_t negativeCount = 0;
    std::size_t nFaces = 0, nCells = 0;
};

class MeshQuality {
public:
    static MeshQualityReport evaluate(const Mesh& mesh);
    static std::string       format(const MeshQualityReport& r);

    /// Render a single histogram as an ASCII bar chart (10 bins of `#`
    /// characters scaled so the tallest bin gets `barWidth` chars).
    static std::string format_histogram(const QualityHistogram& h,
                                        std::string_view        label,
                                        std::size_t             barWidth = 40);

    /// CSV with one row per face: faceId,skewness,nonOrthoDeg.  First line
    /// is the header.
    static std::string to_csv_faces(const MeshQualityReport& r);

    /// CSV with one row per cell: cellId,aspectRatio,negativeVolume.  First
    /// line is the header.
    static std::string to_csv_cells(const MeshQualityReport& r);

    /// Per-boundary-zone aggregation of face-based metrics.  One entry per
    /// distinct non-interior zone id encountered on the mesh boundary.
    static std::vector<struct ZoneQualityStats>
    per_zone_stats(const Mesh& mesh, const MeshQualityReport& r);

    /// Human-readable summary table of the per-zone stats.
    static std::string format_per_zone(const std::vector<struct ZoneQualityStats>& s);
};

struct ZoneQualityStats {
    ZoneId      id          = 0;
    std::string name;
    std::size_t nFaces      = 0;
    double      maxSkewness = 0.0;
    double      meanSkewness = 0.0;
    double      maxNonOrtho = 0.0;
    double      meanNonOrtho = 0.0;
};

}  // namespace simall::meshing
