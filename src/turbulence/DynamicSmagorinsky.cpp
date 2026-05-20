// =============================================================================
// SimAll Beta - Turbulence Subsystem
// File   : src/turbulence/DynamicSmagorinsky.cpp
// =============================================================================
#include "turbulence/DynamicSmagorinsky.hpp"

#include "core/Logger.hpp"
#include "solver/Gradient.hpp"

#include <algorithm>
#include <cmath>

namespace simall::turbulence
{

namespace
{
// Test-filter a per-cell scalar by averaging over the cell and its
// face-neighbour cells (volume-weighted box filter ≈ Δ̂ = 2Δ).
template <class F, class OUT> void test_filter(const meshing::Mesh& m, F fn, OUT& out)
{
    const auto& C = m.cells();
    const auto& Fc = m.faces();
    const std::size_t nC = C.size();
    out.assign(nC, 0.0);
    util::aligned_vector<double> wsum(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c) {
        out[c] += C.volume[c] * fn(c);
        wsum[c] += C.volume[c];
    }
    const std::size_t nF = Fc.size();
    for (std::size_t f = 0; f < nF; ++f) {
        const auto o = Fc.owner[f];
        const auto n = Fc.neighbor[f];
        if (n == meshing::kBoundaryCell)
            continue;
        const double vo = C.volume[o], vn = C.volume[n];
        out[o] += vn * fn(n);
        out[n] += vo * fn(o);
        wsum[o] += vn;
        wsum[n] += vo;
    }
    for (std::size_t c = 0; c < nC; ++c)
        out[c] /= std::max(wsum[c], 1e-30);
}
} // namespace

void DynamicSmagorinsky_LES::initialize(meshing::Mesh& m, solver::FieldRegistry& f)
{
    mesh_ = &m;
    const std::size_t nC = m.cells().size();
    f.scalar("mut", nC);
    f.scalar("Cs2", nC);
    mut_.assign(nC, 0.0);
    Cs2_.assign(nC, 0.0);
    delta_.assign(nC, 0.0);
    for (std::size_t c = 0; c < nC; ++c)
        delta_[c] = std::cbrt(std::max(m.cells().volume[c], 1.0e-30));
    SIMALL_LOG_INFO(
        "Turbulence", "Dynamic Smagorinsky initialised (cells=", nC, ", Cs_max=", csMax_, ")");
}

void DynamicSmagorinsky_LES::solve(double /*dt*/, solver::FieldRegistry& f)
{
    if (!mesh_)
        return;
    const std::size_t nC = mesh_->cells().size();
    const auto& U = *f.find_vector("U");
    auto& mut = *f.find_scalar("mut");
    auto& csF = *f.find_scalar("Cs2");

    // 1) Velocity gradient → strain rate S_ij, |S|.
    solver::LeastSquaresGradient G(*mesh_);
    solver::VectorField gUx, gUy, gUz;
    G.evaluate(U.x, gUx);
    G.evaluate(U.y, gUy);
    G.evaluate(U.z, gUz);

    util::aligned_vector<double> S11(nC), S22(nC), S33(nC), S12(nC), S13(nC), S23(nC), Smag(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        S11[c] = gUx.x[c];
        S22[c] = gUy.y[c];
        S33[c] = gUz.z[c];
        S12[c] = 0.5 * (gUx.y[c] + gUy.x[c]);
        S13[c] = 0.5 * (gUx.z[c] + gUz.x[c]);
        S23[c] = 0.5 * (gUy.z[c] + gUz.y[c]);
        const double SS = 2.0
                          * (S11[c] * S11[c] + S22[c] * S22[c] + S33[c] * S33[c]
                             + 2.0 * (S12[c] * S12[c] + S13[c] * S13[c] + S23[c] * S23[c]));
        Smag[c] = std::sqrt(std::max(SS, 0.0));
    }

    // 2) Test-filtered fields.  Δ̂² / Δ² = 4.
    util::aligned_vector<double> ux_t, uy_t, uz_t;
    test_filter(*mesh_, [&](std::size_t c) { return U.x[c]; }, ux_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.y[c]; }, uy_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.z[c]; }, uz_t);

    util::aligned_vector<double> uxx_t, uyy_t, uzz_t, uxy_t, uxz_t, uyz_t;
    test_filter(*mesh_, [&](std::size_t c) { return U.x[c] * U.x[c]; }, uxx_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.y[c] * U.y[c]; }, uyy_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.z[c] * U.z[c]; }, uzz_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.x[c] * U.y[c]; }, uxy_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.x[c] * U.z[c]; }, uxz_t);
    test_filter(*mesh_, [&](std::size_t c) { return U.y[c] * U.z[c]; }, uyz_t);

    // L_ij = (u_i u_j)̃ - ũ_i ũ_j      (deviatoric — subtract trace later)
    util::aligned_vector<double> L11(nC), L22(nC), L33(nC), L12(nC), L13(nC), L23(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        L11[c] = uxx_t[c] - ux_t[c] * ux_t[c];
        L22[c] = uyy_t[c] - uy_t[c] * uy_t[c];
        L33[c] = uzz_t[c] - uz_t[c] * uz_t[c];
        L12[c] = uxy_t[c] - ux_t[c] * uy_t[c];
        L13[c] = uxz_t[c] - ux_t[c] * uz_t[c];
        L23[c] = uyz_t[c] - uy_t[c] * uz_t[c];
        const double tr = (L11[c] + L22[c] + L33[c]) / 3.0;
        L11[c] -= tr;
        L22[c] -= tr;
        L33[c] -= tr;
    }

    // Test-filter |S|S_ij  (each component).
    util::aligned_vector<double> SS11_t(nC), SS22_t(nC), SS33_t(nC), SS12_t(nC), SS13_t(nC),
        SS23_t(nC);
    test_filter(*mesh_, [&](std::size_t c) { return Smag[c] * S11[c]; }, SS11_t);
    test_filter(*mesh_, [&](std::size_t c) { return Smag[c] * S22[c]; }, SS22_t);
    test_filter(*mesh_, [&](std::size_t c) { return Smag[c] * S33[c]; }, SS33_t);
    test_filter(*mesh_, [&](std::size_t c) { return Smag[c] * S12[c]; }, SS12_t);
    test_filter(*mesh_, [&](std::size_t c) { return Smag[c] * S13[c]; }, SS13_t);
    test_filter(*mesh_, [&](std::size_t c) { return Smag[c] * S23[c]; }, SS23_t);

    // Test-filtered strain rate from filtered velocity → S̃_ij.
    // Approximation: rebuild S̃_ij via filtering of S_ij (compatible with
    // smooth flows; full re-gradient on Ũ would be more accurate but costlier).
    util::aligned_vector<double> s11t(nC), s22t(nC), s33t(nC), s12t(nC), s13t(nC), s23t(nC);
    test_filter(*mesh_, [&](std::size_t c) { return S11[c]; }, s11t);
    test_filter(*mesh_, [&](std::size_t c) { return S22[c]; }, s22t);
    test_filter(*mesh_, [&](std::size_t c) { return S33[c]; }, s33t);
    test_filter(*mesh_, [&](std::size_t c) { return S12[c]; }, s12t);
    test_filter(*mesh_, [&](std::size_t c) { return S13[c]; }, s13t);
    test_filter(*mesh_, [&](std::size_t c) { return S23[c]; }, s23t);

    util::aligned_vector<double> LM(nC), MM(nC);
    for (std::size_t c = 0; c < nC; ++c) {
        const double SmagT = std::sqrt(
            std::max(0.0,
                     2.0
                         * (s11t[c] * s11t[c] + s22t[c] * s22t[c] + s33t[c] * s33t[c]
                            + 2.0 * (s12t[c] * s12t[c] + s13t[c] * s13t[c] + s23t[c] * s23t[c]))));
        const double D2 = delta_[c] * delta_[c];
        const double Dh2 = 4.0 * D2; // Δ̂² = (2Δ)²
        const double M11 = 2.0 * (Dh2 * SmagT * s11t[c] - D2 * SS11_t[c]);
        const double M22 = 2.0 * (Dh2 * SmagT * s22t[c] - D2 * SS22_t[c]);
        const double M33 = 2.0 * (Dh2 * SmagT * s33t[c] - D2 * SS33_t[c]);
        const double M12 = 2.0 * (Dh2 * SmagT * s12t[c] - D2 * SS12_t[c]);
        const double M13 = 2.0 * (Dh2 * SmagT * s13t[c] - D2 * SS13_t[c]);
        const double M23 = 2.0 * (Dh2 * SmagT * s23t[c] - D2 * SS23_t[c]);
        LM[c] = L11[c] * M11 + L22[c] * M22 + L33[c] * M33
                + 2.0 * (L12[c] * M12 + L13[c] * M13 + L23[c] * M23);
        MM[c] = M11 * M11 + M22 * M22 + M33 * M33 + 2.0 * (M12 * M12 + M13 * M13 + M23 * M23);
    }
    // Average LM and MM over neighbours (stabiliser).
    util::aligned_vector<double> LMb, MMb;
    test_filter(*mesh_, [&](std::size_t c) { return LM[c]; }, LMb);
    test_filter(*mesh_, [&](std::size_t c) { return MM[c]; }, MMb);

    const double csMax2 = csMax_ * csMax_;
    for (std::size_t c = 0; c < nC; ++c) {
        double cs2 = (MMb[c] > 1.0e-30) ? (LMb[c] / MMb[c]) : 0.0;
        cs2 = std::clamp(cs2, 0.0, csMax2);
        Cs2_[c] = cs2;
        csF[c] = cs2;
        const double Ls = std::sqrt(cs2) * delta_[c];
        mut[c] = rho_ * Ls * Ls * Smag[c];
        mut_[c] = mut[c];
    }
}

} // namespace simall::turbulence
