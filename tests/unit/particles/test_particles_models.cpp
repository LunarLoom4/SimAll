// =============================================================================
// SimAll Beta — Particles Unit Tests (Week 9)
// File   : tests/unit/particles/test_particles_models.cpp
//
// Sanity tests for the eight new W9 modules:
//   1. RosinRammlerInjector   (size distribution sampling)
//   2. BlobInjector           (Reitz 1987)
//   3. ConeInjector           (solid / hollow cone)
//   4. ParticleEvaporation    (Spalding + Abramzon-Sirignano)
//   5. ParticleHeating        (Ranz-Marshall + radiation)
//   6. ParcelCoalescence      (O'Rourke 1981)
//   7. WallFilm               (Stanton-Rutland + Bai-Gosman)
//   8. ChargedParticleField   (Coulomb / Lorentz)
// =============================================================================
#include "meshing/MeshStorage.hpp"
#include "particles/BlobInjector.hpp"
#include "particles/ChargedParticleField.hpp"
#include "particles/ConeInjector.hpp"
#include "particles/LagrangianTracker.hpp"
#include "particles/ParcelCoalescence.hpp"
#include "particles/ParticleEvaporation.hpp"
#include "particles/ParticleHeating.hpp"
#include "particles/RosinRammlerInjector.hpp"
#include "particles/WallFilm.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

using namespace simall;
using Catch::Approx;

namespace
{
meshing::Mesh build_single_cell_mesh()
{
    meshing::Mesh m;
    auto& C = m.cells();
    auto& F = m.faces();
    C.volume.assign(1, 1.0);
    C.centroidX.assign(1, 0.5);
    C.centroidY.assign(1, 0.5);
    C.centroidZ.assign(1, 0.5);
    C.faceOffsets.assign(2, 0);
    C.faceOffsets[1] = 2;
    C.faceIndices.assign(2, 0);
    C.faceIndices[0] = 0;
    C.faceIndices[1] = 1;
    F.owner.assign(2, 0);
    F.neighbor.assign(2, meshing::kBoundaryCell);
    F.areaX.assign(2, 0.0);
    F.areaY.assign(2, 0.0);
    F.areaZ.assign(2, 1.0);
    F.areaZ[0] = -1.0;
    F.centroidX.assign(2, 0.5);
    F.centroidY.assign(2, 0.5);
    F.centroidZ.assign(2, 0.0);
    F.centroidZ[1] = 1.0;
    F.boundaryZone.assign(2, 1);
    F.boundaryZone[1] = 2;
    return m;
}
} // namespace

TEST_CASE("RosinRammlerInjector: diameter samples lie in bounds",
          "[particles][injector][rosin][week9]")
{
    particles::RosinRammlerInjector inj;
    particles::RosinRammlerProps p;
    p.X = 5e-5;
    p.n = 3.5;
    p.d_min = 1e-6;
    p.d_max = 5e-4;
    inj.initialize(p);
    for (int i = 0; i < 1000; ++i) {
        const double d = inj.sample_diameter();
        REQUIRE(d >= p.d_min);
        REQUIRE(d <= p.d_max);
    }
}

TEST_CASE("BlobInjector: U_inj = Cd √(2 Δp / ρ_l)", "[particles][injector][blob][week9]")
{
    particles::BlobInjector inj;
    particles::BlobInjectorProps p;
    p.Cd = 0.8;
    p.p_inj = 1.0e8;
    p.p_amb = 0.0;
    p.rho_l = 1000.0;
    inj.initialize(p);
    const double expected = 0.8 * std::sqrt(2.0 * 1.0e8 / 1000.0);
    REQUIRE(inj.injection_velocity() == Approx(expected).epsilon(1e-12));
}

TEST_CASE("BlobInjector: injects only inside duration window", "[particles][injector][blob][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.225, 1.81e-5, {});
    particles::BlobInjector inj;
    particles::BlobInjectorProps p;
    p.duration = 1e-3;
    inj.initialize(p);
    REQUIRE(inj.inject(1e-4, 10, 0.5e-3, tr) == 10); // inside window
    REQUIRE(inj.inject(1e-4, 10, 2.0e-3, tr) == 0);  // past window
}

TEST_CASE("ConeInjector: parcels lie within outer cone angle", "[particles][injector][cone][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.225, 1.81e-5, {});
    particles::ConeInjector inj;
    particles::ConeInjectorProps p;
    p.kind = particles::ConeKind::Solid;
    p.theta_outer = 0.4;
    p.speed = 25.0;
    p.axis = {1, 0, 0};
    inj.initialize(p);
    const std::size_t N = inj.inject(1e-4, 200, tr);
    REQUIRE(N == 200);
    const double cosmin = std::cos(p.theta_outer);
    for (const auto& pt : tr.particles()) {
        const double mag = std::sqrt(pt.v.x * pt.v.x + pt.v.y * pt.v.y + pt.v.z * pt.v.z);
        REQUIRE(mag == Approx(p.speed).epsilon(1e-9));
        const double cosA = pt.v.x / mag; // axis is +x
        REQUIRE(cosA >= cosmin - 1e-9);
    }
}

TEST_CASE("ParticleEvaporation: parcel mass decreases when gas is dry & hot",
          "[particles][evaporation][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.0, 1.81e-5, {});
    particles::ParticleState pt;
    pt.x = {0.5, 0.5, 0.5};
    pt.v = {0, 0, 0};
    pt.cell = 0;
    pt.mass = 4.0 / 3.0 * 3.14159265 * std::pow(5e-5, 3) * 998.2;
    pt.active = true;
    tr.mutable_particles().push_back(pt);

    particles::ParticleEvaporation evap;
    particles::EvaporationProps ep;
    ep.d_min = 1e-9;
    evap.initialize(m, ep);
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1);
    T[0] = 600.0;
    auto& P = F.scalar("p", 1);
    P[0] = 101325.0;
    auto& Y = F.scalar("Y_vapor", 1);
    Y[0] = 0.0;
    const double m_before = tr.mutable_particles()[0].mass;
    const double m_evap = evap.apply(1e-3, tr, F, 350.0);
    REQUIRE(m_evap > 0.0);
    REQUIRE(tr.mutable_particles()[0].mass < m_before);
    REQUIRE(F.find_scalar("S_vapor")->at(0) > 0.0);
}

TEST_CASE("ParticleHeating: parcel temperature rises toward gas", "[particles][heating][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.0, 1.81e-5, {});
    particles::ParticleState pt;
    pt.x = {0.5, 0.5, 0.5};
    pt.v = {10, 0, 0};
    pt.cell = 0;
    pt.mass = 1e-9;
    pt.active = true;
    tr.mutable_particles().push_back(pt);

    particles::ParticleHeating heat;
    heat.initialize(m, {});
    heat.resize_to(1);
    heat.set_parcel_temperature(0, 300.0);
    solver::FieldRegistry F;
    auto& T = F.scalar("T", 1);
    T[0] = 1500.0;
    heat.apply(1e-3, tr, F);
    REQUIRE(heat.parcel_temperature(0) > 300.0);
    REQUIRE(F.find_scalar("S_particle_energy")->at(0) < 0.0); // gas loses heat
}

TEST_CASE("ParcelCoalescence: stochastic algorithm decreases active count",
          "[particles][coalescence][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.0, 1.81e-5, {});
    // Two heavy collinear parcels, low Weber → coalesce.
    auto& parts = tr.mutable_particles();
    for (int k = 0; k < 100; ++k) {
        particles::ParticleState s;
        s.x = {0.5, 0.5, 0.5};
        s.v = (k % 2) ? util::Vec3d{0.001, 0, 0} : util::Vec3d{-0.001, 0, 0};
        s.cell = 0;
        s.mass = 1e-5;
        s.active = true;
        parts.push_back(s);
    }
    particles::ParcelCoalescence co;
    particles::CoalescenceProps cp;
    cp.We_crit = 1e9; // force coalesce
    co.initialize(m, cp);
    const std::size_t events = co.apply(1e-2, tr);
    REQUIRE(events > 0);
    std::size_t live = 0;
    for (const auto& s : tr.particles())
        if (s.active)
            ++live;
    REQUIRE(live < 100);
}

TEST_CASE("WallFilm: deposits mass on a sticking impact", "[particles][wallfilm][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.0, 1.81e-5, {});
    particles::ParticleState pt;
    pt.x = {0.5, 0.5, 0.99};
    pt.v = {0, 0, 1.0};
    pt.cell = 0;
    pt.mass = 1e-7;
    pt.active = true;
    tr.mutable_particles().push_back(pt);
    particles::WallFilm film;
    particles::WallFilmProps fp;
    fp.wallZones = {2}; // top face zone
    film.initialize(m, fp);
    const double dep = film.deposit(1e-3, tr);
    REQUIRE(dep > 0.0);
    REQUIRE(film.face_count() == 2);
}

TEST_CASE("ChargedParticleField: F = qE accelerates parcel along E", "[particles][charged][week9]")
{
    auto m = build_single_cell_mesh();
    particles::LagrangianTracker tr;
    tr.initialize(m, 1.0, 1.81e-5, {});
    particles::ParticleState pt;
    pt.x = {0.5, 0.5, 0.5};
    pt.v = {0, 0, 0};
    pt.cell = 0;
    pt.mass = 1e-12;
    pt.active = true;
    tr.mutable_particles().push_back(pt);

    particles::ChargedParticleField cpf;
    particles::ChargedParticleProps cp;
    cp.charge_per_parcel = 1e-15;
    cp.E_uniform = {1e6, 0, 0};
    cpf.initialize(m, cp);
    solver::FieldRegistry F;
    const double dvMax = cpf.apply(1e-6, tr, F);
    REQUIRE(dvMax > 0.0);
    REQUIRE(tr.particles()[0].v.x > 0.0); // pushed along +E
}
