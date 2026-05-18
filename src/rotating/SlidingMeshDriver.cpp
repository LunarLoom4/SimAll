// =============================================================================
// SimAll Beta - Rotating Subsystem
// File   : src/rotating/SlidingMeshDriver.cpp
// =============================================================================
#include "rotating/SlidingMeshDriver.hpp"

#include <algorithm>
#include <cmath>

namespace simall::rotating {

namespace {

std::array<double, 3> normalise(const std::array<double, 3>& v) {
    const double n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    return n > 0.0 ? std::array<double, 3>{v[0]/n, v[1]/n, v[2]/n}
                   : std::array<double, 3>{0.0, 0.0, 1.0};
}

std::array<double, 3> rotate_about(const std::array<double, 3>& p,
                                    const std::array<double, 3>& origin,
                                    const std::array<double, 3>& axisN,
                                    double theta) {
    const double c = std::cos(theta), s = std::sin(theta), C = 1.0 - c;
    const double R[9] = {
        c + axisN[0]*axisN[0]*C,        axisN[0]*axisN[1]*C - axisN[2]*s, axisN[0]*axisN[2]*C + axisN[1]*s,
        axisN[1]*axisN[0]*C + axisN[2]*s, c + axisN[1]*axisN[1]*C,        axisN[1]*axisN[2]*C - axisN[0]*s,
        axisN[2]*axisN[0]*C - axisN[1]*s, axisN[2]*axisN[1]*C + axisN[0]*s, c + axisN[2]*axisN[2]*C };
    const double dx = p[0] - origin[0];
    const double dy = p[1] - origin[1];
    const double dz = p[2] - origin[2];
    return {
        origin[0] + R[0]*dx + R[1]*dy + R[2]*dz,
        origin[1] + R[3]*dx + R[4]*dy + R[5]*dz,
        origin[2] + R[6]*dx + R[7]*dy + R[8]*dz
    };
}

double sqr_dist(const std::array<double, 3>& a, const std::array<double, 3>& b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return dx*dx + dy*dy + dz*dz;
}

}  // namespace

SlidingMeshDriver::SlidingMeshDriver(std::array<double, 3> origin,
                                       std::array<double, 3> dir,
                                       double omega)
    : axisOrigin_(origin), axisDir_(normalise(dir)), omega_(omega) {}

void SlidingMeshDriver::set_rotor_faces (std::vector<InterfaceFace> f) { rotor_  = std::move(f); }
void SlidingMeshDriver::set_stator_faces(std::vector<InterfaceFace> f) { stator_ = std::move(f); }

double SlidingMeshDriver::advance_rotor(double dt) {
    angle_ += omega_ * dt;
    return angle_;
}

InterfaceFace SlidingMeshDriver::rotated_rotor_face(const InterfaceFace& f) const {
    InterfaceFace r = f;
    r.centroid = rotate_about(f.centroid, axisOrigin_, axisDir_, angle_);
    r.normal   = rotate_about(f.normal,
                               {0.0, 0.0, 0.0},      // pure direction → about origin
                               axisDir_, angle_);
    return r;
}

SlidingInterfaceMap SlidingMeshDriver::build_interface_map() const {
    SlidingInterfaceMap m;
    if (rotor_.empty() || stator_.empty()) return m;
    constexpr std::size_t K = 3;
    m.entries.reserve(rotor_.size() * K);

    for (const auto& rf : rotor_) {
        const auto rfR = rotated_rotor_face(rf);
        // Find K nearest stator face centroids.
        std::vector<std::pair<double, FaceIdx>> dists;
        dists.reserve(stator_.size());
        for (const auto& sf : stator_) dists.emplace_back(sqr_dist(rfR.centroid, sf.centroid), sf.id);
        const std::size_t nk = std::min(K, dists.size());
        std::partial_sort(dists.begin(), dists.begin() + nk, dists.end(),
                          [](const auto& a, const auto& b){ return a.first < b.first; });
        double wsum = 0.0;
        for (std::size_t i = 0; i < nk; ++i) wsum += 1.0 / std::max(dists[i].first, 1e-30);
        for (std::size_t i = 0; i < nk; ++i) {
            const double w = (1.0 / std::max(dists[i].first, 1e-30)) / wsum;
            m.entries.push_back({ rfR.id, dists[i].second, w });
        }
    }
    return m;
}

}  // namespace simall::rotating
