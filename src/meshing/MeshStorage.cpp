#include "meshing/MeshStorage.hpp"

#include <numeric>

namespace simall::meshing
{

ZoneId Mesh::add_zone(std::string name, bool b)
{
    ZoneId id = next_zone_++;
    zones_[id] = ZoneInfo{id, std::move(name), b};
    return id;
}

ZoneInfo* Mesh::find_zone(ZoneId id)
{
    auto it = zones_.find(id);
    return it == zones_.end() ? nullptr : &it->second;
}

void Mesh::compute_geometry()
{
    const std::size_t nF = faces_.size();
    faces_.areaX.assign(nF, 0);
    faces_.areaY.assign(nF, 0);
    faces_.areaZ.assign(nF, 0);
    faces_.centroidX.assign(nF, 0);
    faces_.centroidY.assign(nF, 0);
    faces_.centroidZ.assign(nF, 0);

    for (std::size_t f = 0; f < nF; ++f) {
        const std::int32_t s = faces_.nodeOffsets[f];
        const std::int32_t e = faces_.nodeOffsets[f + 1];
        if (e - s < 3)
            continue;

        // Centroid as nodes average (planar approximation; adequate for FV
        // metric pre-pass before exact polygon decomposition).
        util::Vec3d c{};
        for (std::int32_t i = s; i < e; ++i) {
            NodeId n = faces_.nodeIndices[i];
            c.x += nodes_.x[n];
            c.y += nodes_.y[n];
            c.z += nodes_.z[n];
        }
        const double inv = 1.0 / double(e - s);
        c.x *= inv;
        c.y *= inv;
        c.z *= inv;
        faces_.centroidX[f] = c.x;
        faces_.centroidY[f] = c.y;
        faces_.centroidZ[f] = c.z;

        // Sum of triangle fan areas (Newell's formula for arbitrary polygon).
        util::Vec3d area{};
        for (std::int32_t i = s; i < e; ++i) {
            NodeId a = faces_.nodeIndices[i];
            NodeId b = faces_.nodeIndices[(i + 1 < e) ? i + 1 : s];
            util::Vec3d pa{nodes_.x[a], nodes_.y[a], nodes_.z[a]};
            util::Vec3d pb{nodes_.x[b], nodes_.y[b], nodes_.z[b]};
            util::Vec3d cross = (pa - c).cross(pb - c);
            area = area + cross * 0.5;
        }
        faces_.areaX[f] = area.x;
        faces_.areaY[f] = area.y;
        faces_.areaZ[f] = area.z;
    }

    // Cell volume via divergence theorem: V = (1/3) Σ (x_face · n_face)
    const std::size_t nC = cells_.size();
    cells_.volume.assign(nC, 0);
    cells_.centroidX.assign(nC, 0);
    cells_.centroidY.assign(nC, 0);
    cells_.centroidZ.assign(nC, 0);

    for (std::size_t c = 0; c < nC; ++c) {
        double V = 0;
        util::Vec3d Cw{};
        const std::int32_t s = cells_.faceOffsets[c];
        const std::int32_t e = cells_.faceOffsets[c + 1];
        for (std::int32_t i = s; i < e; ++i) {
            FaceId f = cells_.faceIndices[i];
            util::Vec3d xf{faces_.centroidX[f], faces_.centroidY[f], faces_.centroidZ[f]};
            util::Vec3d af{faces_.areaX[f], faces_.areaY[f], faces_.areaZ[f]};
            const double sgn = (faces_.owner[f] == c) ? 1.0 : -1.0;
            double dv = sgn * xf.dot(af) / 3.0;
            V += dv;
            Cw = Cw + xf * dv * (3.0 / 4.0);
        }
        cells_.volume[c] = V;
        if (V != 0) {
            cells_.centroidX[c] = Cw.x / V;
            cells_.centroidY[c] = Cw.y / V;
            cells_.centroidZ[c] = Cw.z / V;
        }
    }
}

} // namespace simall::meshing
