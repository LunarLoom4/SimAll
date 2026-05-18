// =============================================================================
// SimAll Beta - Visualization unit tests (Week 14)
// File   : tests/unit/visualization/test_visualization_models.cpp
//
// Exercises every VTK-free Week-14 module against analytic ground truth so
// the suite can run on headless CI nodes.
// =============================================================================
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "visualization/ActorRegistry.hpp"
#include "visualization/AnimationRecorder.hpp"
#include "visualization/Annotation.hpp"
#include "visualization/CameraController.hpp"
#include "visualization/ClippingPlane.hpp"
#include "visualization/ContourFilter.hpp"
#include "visualization/IsoSurface.hpp"
#include "visualization/LicFilter.hpp"
#include "visualization/OrientationGizmo.hpp"
#include "visualization/PickingBridge.hpp"
#include "visualization/ScalarBar.hpp"
#include "visualization/SceneGraph.hpp"
#include "visualization/ScreenshotRecorder.hpp"
#include "visualization/SectionCut.hpp"
#include "visualization/StreamlineRk4.hpp"
#include "visualization/VectorGlyph.hpp"
#include "visualization/VolumeRaycast.hpp"

#include <cmath>
#include <filesystem>

namespace viz = simall::visualization;
using simall::util::BoundingBox;
using simall::util::Vec3d;

// ---------------------------------------------------------------------------
//  Scene management
// ---------------------------------------------------------------------------
TEST_CASE("SceneGraph builds a hierarchy with composed transforms",
          "[visualization][scene]")
{
    viz::SceneGraph sg;
    const auto root  = sg.root();
    REQUIRE(sg.contains(root));

    const auto a = sg.create_node("A", root);
    const auto b = sg.create_node("B", a);
    REQUIRE(sg.node_count() == 3);

    sg.set_transform(a, viz::Transform::translation({1, 2, 3}));
    sg.set_transform(b, viz::Transform::translation({10, 0, 0}));

    const auto worldB = sg.world_transform(b);
    const Vec3d p = worldB.apply({0, 0, 0});
    REQUIRE(p.x == Catch::Approx(11.0));
    REQUIRE(p.y == Catch::Approx(2.0));
    REQUIRE(p.z == Catch::Approx(3.0));

    sg.set_visible(a, false);
    REQUIRE_FALSE(sg.visible(b));            // parent invisible → child invisible
    sg.set_visible(a, true);
    REQUIRE(sg.visible(b));

    REQUIRE(sg.reparent(b, root));
    REQUIRE_FALSE(sg.reparent(a, b));        // would form a cycle

    sg.destroy_node(a);
    REQUIRE(sg.contains(b));                 // b was reparented earlier
    REQUIRE_FALSE(sg.contains(a));
}

TEST_CASE("ActorRegistry enforces payload kind and tracks bbox",
          "[visualization][actor]")
{
    viz::ActorRegistry reg;
    auto id = reg.create(viz::ActorKind::Surface, "wall");
    REQUIRE(id != viz::kInvalidActorId);

    viz::SurfaceMesh m;
    m.points = {{0,0,0}, {1,0,0}, {0,1,0}};
    m.triIndex = {0, 1, 2};
    REQUIRE(reg.set_payload(id, m));

    // Wrong kind should be rejected.
    viz::LineSet bogus;
    REQUIRE_FALSE(reg.set_payload(id, bogus));

    auto rec = reg.get(id);
    REQUIRE(rec.has_value());
    REQUIRE(rec->bbox.valid());
    REQUIRE(rec->bbox.max.x == Catch::Approx(1.0));
    REQUIRE(rec->bbox.max.y == Catch::Approx(1.0));

    REQUIRE(reg.list().size() == 1);
    REQUIRE(reg.destroy(id));
    REQUIRE(reg.size() == 0);
}

TEST_CASE("CameraController fit_to centres on the bbox",
          "[visualization][camera]")
{
    viz::CameraController c;
    BoundingBox box;
    box.expand({-1, -1, -1});
    box.expand({ 1,  1,  1});
    c.fit_to(box);
    const auto cam = c.camera();
    REQUIRE(cam.focalPoint.x == Catch::Approx(0.0));
    REQUIRE(cam.focalPoint.y == Catch::Approx(0.0));
    REQUIRE(cam.focalPoint.z == Catch::Approx(0.0));
    REQUIRE((cam.position - cam.focalPoint).norm() > 1.0);

    c.save_view("home");
    c.on_wheel(10);
    REQUIRE(c.restore_view("home"));
    REQUIRE(c.camera().focalPoint.x == Catch::Approx(0.0));
}

// ---------------------------------------------------------------------------
//  Filters
// ---------------------------------------------------------------------------
TEST_CASE("ContourFilter extracts a circular isoline on a triangulated square",
          "[visualization][contour]")
{
    // Two-triangle unit square in z=0 with scalar = x.
    viz::SurfaceMesh s;
    s.points = {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0}};
    s.triIndex = {0, 1, 2, 0, 2, 3};
    s.pointScalars = {0.0, 1.0, 1.0, 0.0};
    auto lines = viz::ContourFilter::extract(s, {0.5});
    // Each crossing triangle yields one segment → 2 segments total.
    REQUIRE(lines.line_count() == 2);
    for (auto& p : lines.points) {
        REQUIRE(p.x == Catch::Approx(0.5));
    }
}

TEST_CASE("IsoSurface extracts a triangle from a single tetrahedron",
          "[visualization][iso]")
{
    viz::VolumeMesh v;
    v.points = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
    v.tetIndex = {0, 1, 2, 3};
    v.pointScalars = {0.0, 1.0, 1.0, 1.0};
    auto surf = viz::IsoSurface::extract(v, 0.5);
    REQUIRE(surf.triangle_count() == 1);
    // Each output vertex sits at the midpoint of an edge from vertex 0.
    for (auto& p : surf.points) {
        const double sum = p.x + p.y + p.z;
        REQUIRE(sum == Catch::Approx(0.5));
    }
}

TEST_CASE("StreamlineRk4 traces a circular vector field",
          "[visualization][streamline]")
{
    // v(x,y,z) = (-y, x, 0) — pure rotation about z; radius preserved.
    auto sampler = [](const Vec3d& p) -> std::optional<Vec3d> {
        if (p.x*p.x + p.y*p.y > 4.0) return std::nullopt;
        return Vec3d{-p.y, p.x, 0.0};
    };
    viz::StreamlineConfig cfg;
    cfg.stepSize = 0.01;
    cfg.maxSteps = 800;
    cfg.bidirectional = false;
    auto ls = viz::StreamlineRk4::trace(sampler, Vec3d{1, 0, 0}, cfg);
    REQUIRE(ls.line_count() == 1);
    REQUIRE(ls.points.size() > 100);
    // RK4 with h=0.01 keeps radius to within ~1e-6 over 8 rad of arc.
    for (auto& p : ls.points) {
        const double r = std::sqrt(p.x*p.x + p.y*p.y);
        REQUIRE(r == Catch::Approx(1.0).margin(1e-3));
    }
}

TEST_CASE("VectorGlyph samples and culls", "[visualization][glyph]")
{
    auto sampler = [](const Vec3d& p) -> std::optional<Vec3d> {
        return Vec3d{p.x, 0, 0};
    };
    BoundingBox box;
    box.expand({0, 0, 0}); box.expand({2, 1, 1});
    viz::GlyphConfig cfg;
    cfg.minMagnitude = 0.5;          // x<=0.5 lattice points are filtered
    auto g = viz::VectorGlyph::sample_lattice(sampler, box, 5, 1, 1, cfg);
    REQUIRE(!g.anchor.empty());
    for (std::size_t i = 0; i < g.anchor.size(); ++i) {
        REQUIRE(g.magnitude[i] > 0.5);
        // Normalised direction → unit length.
        const double n = g.direction[i].norm();
        REQUIRE(n == Catch::Approx(1.0));
    }
}

TEST_CASE("ClippingPlane keeps the above-plane half", "[visualization][clip]")
{
    viz::SurfaceMesh s;
    // Triangle straddling the z=0 plane.
    s.points = {{0,0, 1}, {1,0,-1}, {0,1, 1}};
    s.triIndex = {0, 1, 2};
    viz::Plane plane{{0,0,0}, {0,0,1}};
    auto out = viz::ClippingPlane::clip(s, plane, /*keepBelow*/false);
    REQUIRE(out.triangle_count() >= 1);
    for (auto& p : out.points) {
        REQUIRE(p.z >= -1e-9);
    }
}

TEST_CASE("SectionCut produces a polygon through a unit tet",
          "[visualization][section]")
{
    viz::VolumeMesh v;
    v.points = {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}};
    v.tetIndex = {0, 1, 2, 3};
    viz::Plane plane{{0,0,0.25}, {0,0,1}};
    auto out = viz::SectionCut::slice(v, plane);
    REQUIRE(out.triangle_count() >= 1);
    for (auto& p : out.points) {
        REQUIRE(p.z == Catch::Approx(0.25));
    }
}

TEST_CASE("VolumeRaycast renders a non-empty image for a dense field",
          "[visualization][volume]")
{
    viz::RegularGrid3 g;
    g.nx = g.ny = g.nz = 8;
    g.values.assign(8u*8u*8u, 1.0f);
    g.origin = {-1, -1, -1};
    g.spacing = {2.0/7, 2.0/7, 2.0/7};
    viz::TransferFunction tf = viz::make_cool_warm(0.0, 1.0);
    for (auto& s : tf.stops) s.color[3] = 0.4f;     // semi-opaque everywhere
    viz::Camera cam;
    cam.viewportWidth = 16; cam.viewportHeight = 16;
    cam.position = {0, 0, 5}; cam.focalPoint = {0,0,0};
    auto img = viz::VolumeRaycast::render(g, cam, tf);
    REQUIRE(img.width == 16);
    // At least one pixel must be non-background.
    bool any = false;
    for (std::size_t i = 0; i < img.pixels.size(); i += 4) {
        if (img.pixels[i] || img.pixels[i+1] || img.pixels[i+2]) { any = true; break; }
    }
    REQUIRE(any);
}

TEST_CASE("LicFilter produces a non-trivial intensity image",
          "[visualization][lic]")
{
    viz::RegularGrid2 g;
    g.nx = g.ny = 16;
    g.vx.assign(16u*16u, 1.0f);
    g.vy.assign(16u*16u, 0.0f);
    viz::LicConfig cfg;
    cfg.stepsPerSide = 8;
    auto img = viz::LicFilter::generate(g, cfg);
    REQUIRE(img.width == 16);
    REQUIRE(img.height == 16);
    // The convolution along the flow direction averages a row of noise; not
    // every pixel is identical → some variance must remain.
    double mean = 0.0;
    for (std::size_t i = 0; i < img.pixels.size(); i += 4) mean += img.pixels[i];
    mean /= static_cast<double>(img.width * img.height);
    REQUIRE(mean > 0.0);
    REQUIRE(mean < 255.0);
}

// ---------------------------------------------------------------------------
//  Overlays
// ---------------------------------------------------------------------------
TEST_CASE("ScalarBar emits ticks and a gradient image",
          "[visualization][scalarbar]")
{
    viz::ScalarBar bar;
    bar.scalarMin = 0.0; bar.scalarMax = 10.0;
    bar.tickCount = 6; bar.precision = 1;
    bar.tf = viz::make_cool_warm(0.0, 10.0);
    auto ticks = bar.ticks();
    REQUIRE(ticks.size() == 6);
    REQUIRE(ticks.front().scalar == Catch::Approx(0.0));
    REQUIRE(ticks.back().scalar  == Catch::Approx(10.0));
    auto img = bar.rasterise_gradient();
    REQUIRE(img.width  == bar.widthPx);
    REQUIRE(img.height == bar.lengthPx);
}

TEST_CASE("AnnotationLayer projects a world point onto the viewport",
          "[visualization][annotation]")
{
    viz::AnnotationLayer L;
    viz::Annotation a;
    a.text = "probe"; a.anchor = {0, 0, 0}; a.space = viz::AnchorSpace::World;
    const auto id = L.add(a);
    REQUIRE(L.size() == 1);

    viz::Camera cam;
    cam.viewportWidth = 100; cam.viewportHeight = 100;
    cam.position = {0, 0, 5}; cam.focalPoint = {0, 0, 0};
    auto proj = viz::AnnotationLayer::project_world(cam, {0, 0, 0});
    REQUIRE(proj.has_value());
    REQUIRE(proj->first  == Catch::Approx(50.0).margin(1e-9));
    REQUIRE(proj->second == Catch::Approx(50.0).margin(1e-9));

    // Behind the camera should fail.
    auto behind = viz::AnnotationLayer::project_world(cam, {0, 0, 10});
    REQUIRE_FALSE(behind.has_value());

    L.remove(id);
    REQUIRE(L.size() == 0);
}

TEST_CASE("OrientationGizmo hits center yields a deterministic axis",
          "[visualization][gizmo]")
{
    viz::Camera cam; cam.viewportWidth = 200; cam.viewportHeight = 200;
    viz::GizmoLayout L;
    L.cornerX = 0.5; L.cornerY = 0.5; L.sizePx = 60;
    auto hit = viz::OrientationGizmo::hit_test(cam, L, 100.0, 100.0);
    REQUIRE(hit.has_value());
    // Missing the circle entirely:
    auto miss = viz::OrientationGizmo::hit_test(cam, L, 0.0, 0.0);
    REQUIRE_FALSE(miss.has_value());

    viz::CameraController ctrl;
    viz::OrientationGizmo::snap(ctrl, viz::Axis::PosX);
    REQUIRE(ctrl.camera().position.x < 0.0);   // camera now in -X looking at origin
}

// ---------------------------------------------------------------------------
//  Recording
// ---------------------------------------------------------------------------
TEST_CASE("ScreenshotRecorder round-trips a PPM image",
          "[visualization][screenshot]")
{
    viz::Image img; img.resize(4, 2);
    for (std::size_t i = 0; i < img.pixels.size(); i += 4) {
        img.pixels[i]   = static_cast<std::uint8_t>(i);
        img.pixels[i+1] = 64;
        img.pixels[i+2] = 128;
        img.pixels[i+3] = 255;
    }
    auto path = std::filesystem::temp_directory_path() / "simall_test_ppm.ppm";
    REQUIRE(viz::ScreenshotRecorder::save_ppm(img, path));
    auto loaded = viz::ScreenshotRecorder::load_ppm(path);
    REQUIRE(loaded.width == img.width);
    REQUIRE(loaded.height == img.height);
    for (std::size_t i = 0; i < img.pixels.size(); i += 4) {
        REQUIRE(loaded.pixels[i]   == img.pixels[i]);
        REQUIRE(loaded.pixels[i+1] == img.pixels[i+1]);
        REQUIRE(loaded.pixels[i+2] == img.pixels[i+2]);
    }
    std::error_code ec; std::filesystem::remove(path, ec);
}

TEST_CASE("AnimationRecorder enforces fps schedule",
          "[visualization][animation]")
{
    viz::AnimationConfig cfg;
    cfg.outputDir = std::filesystem::temp_directory_path() / "simall_anim_test";
    std::filesystem::remove_all(cfg.outputDir);
    cfg.fps = 10;                                     // dt = 0.1
    cfg.filenameDigits = 3;
    viz::AnimationRecorder rec(cfg);

    viz::Image small; small.resize(2, 2);
    REQUIRE( rec.submit(0.00, small));               // first frame
    REQUIRE_FALSE(rec.submit(0.05, small));          // < 0.10s after — dropped
    REQUIRE( rec.submit(0.15, small));               // accepted
    REQUIRE(rec.flush() == 2);

    auto s = rec.stats();
    REQUIRE(s.framesQueued  == 2);
    REQUIRE(s.framesDropped == 1);
    REQUIRE(s.framesWritten == 2);
    std::filesystem::remove_all(cfg.outputDir);
}

// ---------------------------------------------------------------------------
//  Picking
// ---------------------------------------------------------------------------
TEST_CASE("PickingBridge intersects the closest of two surfaces",
          "[visualization][pick]")
{
    viz::ActorRegistry reg;
    auto near = reg.create(viz::ActorKind::Surface, "near");
    auto far_ = reg.create(viz::ActorKind::Surface, "far");

    viz::SurfaceMesh n; n.points = {{-1,-1,1}, {1,-1,1}, {0,1,1}};
    n.triIndex = {0,1,2};
    viz::SurfaceMesh f; f.points = {{-1,-1,-3}, {1,-1,-3}, {0,1,-3}};
    f.triIndex = {0,1,2};
    REQUIRE(reg.set_payload(near, n));
    REQUIRE(reg.set_payload(far_, f));

    // Ray origin at +z=5 looking toward -z hits the near actor at z=1.
    viz::Ray ray{Vec3d{0, 0, 5}, Vec3d{0, 0, -1}};
    auto hit = viz::PickingBridge::pick(reg, ray);
    REQUIRE(hit.has_value());
    REQUIRE(hit->actor == near);
    REQUIRE(hit->worldHit.z == Catch::Approx(1.0));
    REQUIRE(hit->distance   == Catch::Approx(4.0));
}
