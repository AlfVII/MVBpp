// Real-winding geometry ([realwinding]): ONE continuous conductor per (winding, parallel)
// replacing the per-turn closed loops, with every MKF turn position honoured exactly.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "mvb/MagneticBuilder.h"
#include "mvb/Utils.h"
#include "constructive_models/Magnetic.h"
#include "json.hpp"
#include <BRepAlgoAPI_Common.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
#include <cmath>
#include <fstream>
#include <numbers>

using json = nlohmann::json;

namespace {

json loadFixture(const std::string& name) {
    std::ifstream f("testData/" + name);
    if (!f.good()) f = std::ifstream("tests/realwinding_fixtures/" + name);
    REQUIRE(f.good());
    json j = json::parse(f);
    return j.contains("magnetic") ? j.at("magnetic") : j;
}

double shapeVolume(const TopoDS_Shape& s) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(s, props);
    return props.Mass();
}

double commonVolume(const TopoDS_Shape& a, const TopoDS_Shape& b) {
    BRepAlgoAPI_Common common(a, b);
    if (!common.IsDone()) return -1.0;
    return shapeVolume(common.Shape());
}

bool pointInsideShape(const TopoDS_Shape& shape, const gp_Pnt& p, double tol = 1e-9) {
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        BRepClass3d_SolidClassifier cls(TopoDS::Solid(exp.Current()), p, tol);
        if (cls.State() == TopAbs_IN || cls.State() == TopAbs_ON) return true;
    }
    return false;
}

} // namespace

TEST_CASE("Real winding: single-parallel PQ33 becomes one continuous conductor",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_round_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, /*includeBobbin=*/true, /*symmetryPlanes=*/0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       /*paintCoating=*/true, /*emitCoatingShells=*/false,
                                       /*includeInsulation=*/false, /*coreCoatingThickness=*/0.0,
                                       /*useRealWindingGeometry=*/true);

    // Exactly one conductor for the single (winding, parallel); no per-turn ring solids.
    const mvb::NamedShape* conductor = nullptr;
    for (const auto& ns : named) {
        INFO(ns.name);
        REQUIRE_THAT(ns.name, !Catch::Matchers::ContainsSubstring(" turn "));
        if (ns.name == "Primary parallel 0") {
            REQUIRE(conductor == nullptr);
            conductor = &ns;
        }
    }
    REQUIRE(conductor != nullptr);
    REQUIRE(!conductor->shape.IsNull());
    REQUIRE(shapeVolume(conductor->shape) > 0.0);

    // The "nothing moves" regression guard: every MKF turn station lies INSIDE the
    // conductor's copper. Helical wraps pass the exact station at their start phase, so
    // probe the full station RING (a thin torus at (r, y) around the column axis) — it
    // must intersect the copper.
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    REQUIRE(!turnsOpt->empty());
    // The station RING (circle at the turn's exact (r, y)) must enter the copper: its
    // minimum distance to the conductor is 0 when the wire centreline passes through the
    // station somewhere (extrema is far more robust than boolean common on compounds).
    for (const auto& turn : *turnsOpt) {
        const auto& c = turn.get_coordinates();
        REQUIRE(c.size() >= 2);
        gp_Circ stationRing(gp_Ax2(gp_Pnt(0.0, c[1], 0.0), gp_Dir(0, 1, 0)), c[0]);
        TopoDS_Edge ringEdge = BRepBuilderAPI_MakeEdge(stationRing).Edge();
        BRepExtrema_DistShapeShape dist(ringEdge, conductor->shape);
        REQUIRE(dist.IsDone());
        INFO("turn " << turn.get_name() << " station (r=" << c[0] << ", y=" << c[1]
                     << ") ring-distance=" << dist.Value());
        REQUIRE(dist.Value() <= 1e-9);
    }

    // Conductor never intersects the core. Tolerance 1e-10 m^3: the cores are polygon-
    // faceted approximations (the facets dip inside the true round window bore), so
    // micron-sliver grazes at the window border are render artifacts, not collisions —
    // same class the [battery] suite tolerates at 1e-7.
    for (const auto& ns : named) {
        if (&ns == conductor || ns.name.find("Bobbin") != std::string::npos) continue;
        if (ns.name == "Primary parallel 0") continue;
        double v = commonVolume(conductor->shape, ns.shape);
        INFO("overlap of conductor vs " << ns.name << " = " << v);
        REQUIRE(v <= 1e-10);
    }
}

TEST_CASE("Real winding: flag off keeps the per-turn loops unchanged", "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_round_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/false);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched);

    size_t turnSolids = 0;
    bool sawConductor = false;
    for (const auto& ns : named) {
        if (ns.name.find(" turn ") != std::string::npos) ++turnSolids;
        if (ns.name == "Primary parallel 0") sawConductor = true;
    }
    REQUIRE(turnSolids == 34);   // 34 turns x 1 parallel, one closed loop each
    REQUIRE(!sawConductor);
}

TEST_CASE("Real winding: two parallels become two independent conductors", "[realwinding]") {
    auto magneticJson = loadFixture("round_2p_1layer.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true);

    const mvb::NamedShape* p0 = nullptr;
    const mvb::NamedShape* p1 = nullptr;
    for (const auto& ns : named) {
        if (ns.name == "Primary parallel 0") p0 = &ns;
        if (ns.name == "Primary parallel 1") p1 = &ns;
    }
    REQUIRE(p0 != nullptr);
    REQUIRE(p1 != nullptr);
    REQUIRE(shapeVolume(p0->shape) > 0.0);
    REQUIRE(shapeVolume(p1->shape) > 0.0);

    // Parallel conductors are fully independent copper bodies: zero overlap (contact only).
    double v = commonVolume(p0->shape, p1->shape);
    INFO("parallel-parallel common volume = " << v);
    REQUIRE(v <= 1e-12);
}

TEST_CASE("Real winding: multi-parallel LAYER JUMPS throw until MKF reserves a corridor",
          "[realwinding]") {
    // Known, documented limitation: parallels jump layers in lockstep and their crossing
    // jump chords cannot be routed collision-free while every turn stays at its exact MKF
    // position — MKF must first reserve an angular crossover corridor (MKF ABT #187). The
    // collision gate must refuse loudly, never bend the geometry.
    auto magneticJson = loadFixture("realwinding_round_2p.json");   // 8t x 2p -> multi-layer
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);

    mvb::MagneticBuilder builder;
    REQUIRE_THROWS_WITH(
        builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                              mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                              /*useRealWindingGeometry=*/true),
        Catch::Matchers::ContainsSubstring("ABT #187"));
}

TEST_CASE("Real winding: pre-enriched input with the flag on throws", "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_round_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, false);

    // Round-trip to a MAS::Magnetic that carries the geometricalDescription.
    json enrichedJson;
    to_json(enrichedJson, enriched);
    MAS::Magnetic masMagnetic = enrichedJson.get<MAS::Magnetic>();
    REQUIRE(masMagnetic.get_core().value().get_geometrical_description().has_value());

    mvb::MagneticBuilder builder;
    REQUIRE_THROWS_WITH(
        builder.buildAllNamed(masMagnetic, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                              mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                              /*useRealWindingGeometry=*/true),
        Catch::Matchers::ContainsSubstring("geometricalDescription"));
}
