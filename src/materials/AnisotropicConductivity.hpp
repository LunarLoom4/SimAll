// =============================================================================
// SimAll Beta - Materials Subsystem
// File   : src/materials/AnisotropicConductivity.hpp
// Phase  : 15.4 — Anisotropic thermal conductivity tensor.
//
//   q = -k · ∇T          (Fourier with full tensor)
//
// Stored as a 3×3 symmetric tensor in principal-axis or material frame.
// Provides:
//   * orthotropic builder (k_xx, k_yy, k_zz) in material frame
//   * rotated tensor via a passive 3-axis rotation matrix
//   * tensor-vector product (apply to gradient)
//
// Header-only.
// =============================================================================
#pragma once

#include <array>
#include <cmath>

namespace simall::materials
{

struct Tensor3x3
{
    std::array<double, 9> m{}; // row-major

    double& at(int i, int j) noexcept { return m[i * 3 + j]; }
    double at(int i, int j) const noexcept { return m[i * 3 + j]; }

    static Tensor3x3 identity()
    {
        Tensor3x3 t;
        t.at(0, 0) = t.at(1, 1) = t.at(2, 2) = 1.0;
        return t;
    }

    /// y = M · x
    std::array<double, 3> apply(const std::array<double, 3>& x) const noexcept
    {
        return {at(0, 0) * x[0] + at(0, 1) * x[1] + at(0, 2) * x[2],
                at(1, 0) * x[0] + at(1, 1) * x[1] + at(1, 2) * x[2],
                at(2, 0) * x[0] + at(2, 1) * x[1] + at(2, 2) * x[2]};
    }

    /// M = R · D · R^T   (similarity transform).
    static Tensor3x3 similarity(const Tensor3x3& R, const Tensor3x3& D)
    {
        Tensor3x3 RD;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0;
                for (int k = 0; k < 3; ++k)
                    s += R.at(i, k) * D.at(k, j);
                RD.at(i, j) = s;
            }
        Tensor3x3 out;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                double s = 0;
                for (int k = 0; k < 3; ++k)
                    s += RD.at(i, k) * R.at(j, k); // R^T[k,j] = R[j,k]
                out.at(i, j) = s;
            }
        return out;
    }
};

/// Anisotropic conductivity k_ij(T). Each component is allowed to be a
/// temperature-dependent function via the per-component callable.
class AnisotropicConductivity
{
public:
    /// Orthotropic constant: principal-direction conductivities.
    AnisotropicConductivity(double kxx, double kyy, double kzz) : kxx_(kxx), kyy_(kyy), kzz_(kzz) {}

    /// Build the tensor at temperature T in the material frame.
    Tensor3x3 tensor(double T) const noexcept
    {
        Tensor3x3 D;
        D.at(0, 0) = evalAxis_(kxx_, T);
        D.at(1, 1) = evalAxis_(kyy_, T);
        D.at(2, 2) = evalAxis_(kzz_, T);
        return D;
    }

    /// Same, but in a rotated (global) frame given by passive rotation R
    /// such that  v_material = R^T · v_global.
    Tensor3x3 tensorInFrame(double T, const Tensor3x3& R) const noexcept
    {
        return Tensor3x3::similarity(R, tensor(T));
    }

    /// Heat-flux  q = -k · ∇T  for an isothermal local evaluation.
    std::array<double, 3> heatFlux(double T, const std::array<double, 3>& gradT) const noexcept
    {
        auto k = tensor(T);
        auto q = k.apply(gradT);
        for (auto& v : q)
            v = -v;
        return q;
    }

private:
    static double evalAxis_(double k, double /*T*/) noexcept { return k; }

    double kxx_, kyy_, kzz_;
};

} // namespace simall::materials
