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
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <TopoDS_Vertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BOPAlgo_CheckerSI.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt.hxx>
#include <cmath>
#include <fstream>
#include <numbers>
#include <variant>

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

// Exact OCCT self-intersection check on the emitted body. Consecutive-wrap CONTACT is by
// design (a spring); only genuine face-face interference counts. When the conductor is a
// fused/swept single solid, any residual self-intersection is a modelling defect.
bool hasSelfIntersections(const TopoDS_Shape& s) {
    BOPAlgo_CheckerSI checker;
    TopTools_ListOfShape args;
    args.Append(s);
    checker.SetArguments(args);
    checker.Perform();
    return checker.HasErrors();
}

// Facet-wedge bound: the cores are polygon-faceted (n-gon) approximations whose flats dip
// inside the true round window bore by sag = r*(1-cos(pi/n)); a lead legitimately ending
// at the true window border can therefore interpenetrate a facet by up to a wedge of
// volume ~ pi*wireRadius^2 * sag, at each of the two lead ends.
double coreFacetWedgeBound(double wireRadius, double borderRadius, int coreSegments) {
    double sag = borderRadius * (1.0 - std::cos(std::numbers::pi / coreSegments));
    return 2.0 * std::numbers::pi * wireRadius * wireRadius * sag;
}

// All-pairs boolean interference among named bodies (skipping the bobbin, which is
// deliberately cut to yield). Tolerance covers polygon-facet slivers of the cores.
void requireNoPairwiseOverlap(const std::vector<mvb::NamedShape>& named, double tol) {
    for (size_t i = 0; i < named.size(); ++i) {
        if (named[i].name.find("Bobbin") != std::string::npos) continue;
        for (size_t j = i + 1; j < named.size(); ++j) {
            if (named[j].name.find("Bobbin") != std::string::npos) continue;
            double v = commonVolume(named[i].shape, named[j].shape);
            INFO("pairwise overlap '" << named[i].name << "' vs '" << named[j].name
                                      << "' = " << v);
            REQUIRE(v <= tol);
        }
    }
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
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

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

    // "Absolutely no body collides with another or itself":
    // (a) all-pairs boolean interference across every emitted body (tolerance 1e-10 m^3
    //     for the cores' polygon-facet slivers — same class [battery] tolerates at 1e-7);
    // (b) exact OCCT self-intersection check when the conductor fused to a single solid.
    //     When OCCT's tangent-contact fuse defect forces the per-run compound fallback,
    //     the consecutive pieces legitimately overlap at their junctions (the wire's own
    //     crossovers) — for that case the capsule gate + station probes + (a) are the
    //     collision guarantee.
    // Tolerance: the facet-wedge bound for this fixture (wire radius 0.4795 mm, window
    // border ~13.75 mm, 16-gon cores).
    requireNoPairwiseOverlap(named, coreFacetWedgeBound(0.0004795, 0.01375, 16));
    int conductorSolids = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_SOLID); exp.More(); exp.Next())
        ++conductorSolids;
    if (conductorSolids == 1) {
        REQUIRE(!hasSelfIntersections(conductor->shape));
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
    REQUIRE(turnSolids == 16);   // 16 turns x 1 parallel, one closed loop each (2 layers)
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
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);

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

    // Parallel conductors are fully independent copper bodies: zero overlap (contact only),
    // no pairwise interference anywhere, and no self-intersections when fused single.
    double v = commonVolume(p0->shape, p1->shape);
    INFO("parallel-parallel common volume = " << v);
    REQUIRE(v <= 1e-12);
    requireNoPairwiseOverlap(named, 1e-10);
    for (const mvb::NamedShape* c : {p0, p1}) {
        int solids = 0;
        for (TopExp_Explorer exp(c->shape, TopAbs_SOLID); exp.More(); exp.Next()) ++solids;
        if (solids == 1) REQUIRE(!hasSelfIntersections(c->shape));
    }
}

TEST_CASE("Real winding: multi-parallel throws on MKF's overlapping lead rows",
          "[realwinding]") {
    // Under the N+1-crossing model every wrap (including layer transitions) is a spiral
    // between consecutive crossings, so parallel conductors advance in lockstep and the
    // former corridor problem is gone. What still blocks multi-parallel builds is that
    // MKF draws every parallel's terminal lead at the SAME edge row (identical
    // rectangles) — two physical wires on one line. The collision gate must refuse
    // loudly until MKF allocates one row per lead (MKF ABT #229).
    auto magneticJson = loadFixture("realwinding_round_2p.json");   // 8t x 2p -> multi-layer
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);

    mvb::MagneticBuilder builder;
    REQUIRE_THROWS_WITH(
        builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                              mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                              /*useRealWindingGeometry=*/true),
        Catch::Matchers::ContainsSubstring("collision"));
}

namespace {

// Test-side mirror of MagneticBuilder's bobbin resolution (getBobbinProcessed +
// patchBobbinDimensions) so probes use the same column dimensions as the builder.
MAS::CoreBobbinProcessedDescription bobbinPdOf(const OpenMagnetics::Magnetic& m) {
    auto bobbinVar = m.get_coil().get_bobbin();
    auto* b = std::get_if<OpenMagnetics::Bobbin>(&bobbinVar);
    REQUIRE(b != nullptr);
    auto pd = b->get_processed_description();
    REQUIRE(pd.has_value());
    MAS::CoreBobbinProcessedDescription bobbinPd = *pd;
    if (bobbinPd.get_column_width().value_or(0.0) <= 0.0) {
        auto corePd = m.get_core().get_processed_description();
        REQUIRE(corePd.has_value());
        REQUIRE(!corePd->get_columns().empty());
        const auto& col = corePd->get_columns()[0];
        double wall = bobbinPd.get_wall_thickness();
        if (std::isnan(wall) || wall < 0.0) wall = 0.0;
        bobbinPd.set_column_width(col.get_width() / 2.0 + wall);
        bobbinPd.set_column_depth(col.get_depth() / 2.0 + wall);
        bobbinPd.set_column_shape(col.get_shape());
    }
    return bobbinPd;
}

// A crossing must lie ON THE CENTERLINE of the copper: the crossing itself and four
// probes at ±0.99·wireRadius along the two directions perpendicular to the wire's travel
// must all be inside the solid. If the centerline missed the crossing by more than
// 0.01·wireRadius, the probe opposite the offset would fall outside the pipe. (A
// boundary-distance check is unsound here: at the conductor's free ends and at run
// junctions internal cap faces pass exactly through the crossing.)
void requireCrossingOnCenterline(const TopoDS_Shape& conductor, const gp_Pnt& crossing,
                                 double wireRadius, const gp_Dir& perpA, const gp_Dir& perpB,
                                 const std::string& what) {
    INFO(what << " at (" << crossing.X() << "," << crossing.Y() << "," << crossing.Z()
              << "), wire radius " << wireRadius);
    REQUIRE(pointInsideShape(conductor, crossing, 1e-9));
    for (const gp_Dir* d : {&perpA, &perpB}) {
        for (double sgn : {1.0, -1.0}) {
            gp_Pnt probe(crossing.XYZ() + d->XYZ() * (sgn * 0.99 * wireRadius));
            INFO("perpendicular probe at (" << probe.X() << "," << probe.Y() << ","
                                            << probe.Z() << ")");
            REQUIRE(pointInsideShape(conductor, probe, 1e-9));
        }
    }
}

// Outer-footprint wire radius of an MKF-enriched turn (turn.dimensions = outer w/h).
double turnWireRadius(const MAS::Turn& turn) {
    auto dims = turn.get_dimensions();
    REQUIRE(dims.has_value());
    REQUIRE(dims->size() >= 2);
    return std::min((*dims)[0], (*dims)[1]) / 2.0;
}

const mvb::NamedShape* findConductor(const std::vector<mvb::NamedShape>& named,
                                     const std::string& name) {
    const mvb::NamedShape* found = nullptr;
    for (const auto& ns : named) {
        REQUIRE_THAT(ns.name, !Catch::Matchers::ContainsSubstring(" turn "));
        if (ns.name == name) {
            REQUIRE(found == nullptr);
            found = &ns;
        }
    }
    REQUIRE(found != nullptr);
    REQUIRE(!found->shape.IsNull());
    REQUIRE(shapeVolume(found->shape) > 0.0);
    return found;
}

} // namespace

TEST_CASE("Real winding: rectangular-column E core zigzag racetrack conductor",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_rect_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");

    // Every MKF crossing must lie ON the copper. For a rectangular column the crossing
    // sits at the -Z transition-face centre: (0, y, -(x + columnDepth - columnWidth)).
    auto bobbinPd = bobbinPdOf(enriched);
    REQUIRE(bobbinPd.get_column_shape() == MAS::ColumnShape::RECTANGULAR);
    double zoff = bobbinPd.get_column_depth() - bobbinPd.get_column_width().value();
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    for (const auto& turn : *turnsOpt) {
        const auto& c = turn.get_coordinates();
        REQUIRE(c.size() >= 2);
        requireCrossingOnCenterline(conductor->shape, gp_Pnt(0.0, c[1], -(c[0] + zoff)),
                                    turnWireRadius(turn), gp_Dir(0, 1, 0), gp_Dir(0, 0, 1),
                                    "crossing " + turn.get_name());
    }

    // E-core window walls are planar (no facet sag): tangent contact only.
    requireNoPairwiseOverlap(named, 1e-10);
    // A rectangular column's round-wire conductor is built from per-primitive cylinders + torus
    // elbows and fused into ONE FEM-ready solid (the swept-pipe compound would sliver on the fuse).
    int solids = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_SOLID); exp.More(); exp.Next())
        ++solids;
    REQUIRE(solids == 1);
    REQUIRE(!hasSelfIntersections(conductor->shape));
}

TEST_CASE("Real winding: oblong-column EP core stadium conductor", "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_oblong_U.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");

    // Oblong crossing sits at the -Z cap apex: same z = -(x + D - W) mapping.
    auto bobbinPd = bobbinPdOf(enriched);
    REQUIRE(bobbinPd.get_column_shape() == MAS::ColumnShape::OBLONG);
    double zoff = std::max(0.0, bobbinPd.get_column_depth() -
                                    bobbinPd.get_column_width().value());
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    for (const auto& turn : *turnsOpt) {
        const auto& c = turn.get_coordinates();
        REQUIRE(c.size() >= 2);
        requireCrossingOnCenterline(conductor->shape, gp_Pnt(0.0, c[1], -(c[0] + zoff)),
                                    turnWireRadius(turn), gp_Dir(0, 1, 0), gp_Dir(0, 0, 1),
                                    "crossing " + turn.get_name());
    }

    requireNoPairwiseOverlap(named, coreFacetWedgeBound(0.00028, 0.0089, 16));
    int solids = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_SOLID); exp.More(); exp.Next())
        ++solids;
    if (solids == 1) REQUIRE(!hasSelfIntersections(conductor->shape));
}

TEST_CASE("Real winding: toroidal conductor threads the exact inner and outer crossings",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_toroid.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    // Exact core (segments=0): the wall-adjacent ring's tubes touch the bore tangentially
    // along their whole length, so any polygon-faceted bore would interpenetrate them by
    // its facet sag; the exact annulus makes true tangency testable at boolean tolerance.
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");

    // The assembly is counter-rotated to the MAS frame (ring in XY, hole axis Z), so a
    // hole-plane crossing (cx, cy) must lie on the copper at (cx, cy, 0) — the inner
    // crossing of every turn, the outer (additionalCoordinates) of every wrapped turn
    // (the last crossing entry's outer is not wrapped through).
    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    const auto& turns = *turnsOpt;
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& c = turns[i].get_coordinates();
        REQUIRE(c.size() >= 2);
        double wr = turnWireRadius(turns[i]);
        requireCrossingOnCenterline(conductor->shape, gp_Pnt(c[0], c[1], 0.0), wr,
                                    gp_Dir(1, 0, 0), gp_Dir(0, 1, 0),
                                    "turn " + turns[i].get_name() + " inner crossing");
        if (i + 1 < turns.size()) {
            auto add = turns[i].get_additional_coordinates();
            REQUIRE(add.has_value());
            REQUIRE(!add->empty());
            requireCrossingOnCenterline(conductor->shape,
                                        gp_Pnt((*add)[0][0], (*add)[0][1], 0.0), wr,
                                        gp_Dir(1, 0, 0), gp_Dir(0, 1, 0),
                                        "turn " + turns[i].get_name() + " outer crossing");
        }
    }

    // Exact bore: wall contact is true tangency, zero interference within OCCT booleans.
    requireNoPairwiseOverlap(named, 1e-12);
    // This 12-turn toroid sweeps to ONE FEM-ready solid via the simple MakePipe (its framing keeps
    // the round section centred on the high-torsion hole-threading spine, where MakePipeShell's
    // corrected Frenet drifts a crossing off-centre). Denser toroids whose spine MakePipe cannot
    // close fall back to the exact per-run compound -- still crossing-exact, just multi-solid.
    int solids = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_SOLID); exp.More(); exp.Next())
        ++solids;
    REQUIRE(solids == 1);
    REQUIRE(!hasSelfIntersections(conductor->shape));

    // FEM terminal faces: the two free ends of the conductor are flat PLANAR discs (the
    // swept lead cylinder's end cap, no sphere), so downstream FEM can assign a current
    // BC on a planar surface. The wire is round and the whole conductor is otherwise
    // cylinders/tori/revolves/spheres, so the only planar faces are the two terminals.
    int planarFaces = 0;
    for (TopExp_Explorer exp(conductor->shape, TopAbs_FACE); exp.More(); exp.Next()) {
        BRepAdaptor_Surface sa(TopoDS::Face(exp.Current()));
        if (sa.GetType() == GeomAbs_Plane) ++planarFaces;
    }
    INFO("planar (terminal) faces on the conductor = " << planarFaces);
    REQUIRE(planarFaces >= 2);
}

TEST_CASE("Real winding: toroidal RECTANGULAR wire threads the crossings", "[realwinding]") {
    // Rectangular wire on a toroid (single layer): each turn is built from per-primitive rect
    // solids (prisms + revolved poloidal elbows) oriented on the local AZIMUTHAL axis, then fused.
    // Every MKF inner/outer crossing must still lie on the copper (the section's inscribed circle
    // is min(w,h)/2 = turnWireRadius, so the 0.99*r probes stay inside whatever the orientation).
    auto magneticJson = loadFixture("realwinding_toroid_rect.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0,
                                       true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);

    auto turnsOpt = enriched.get_coil().get_turns_description();
    REQUIRE(turnsOpt.has_value());
    const auto& turns = *turnsOpt;
    for (size_t i = 0; i < turns.size(); ++i) {
        const auto& c = turns[i].get_coordinates();
        REQUIRE(c.size() >= 2);
        double wr = turnWireRadius(turns[i]);
        requireCrossingOnCenterline(conductor->shape, gp_Pnt(c[0], c[1], 0.0), wr,
                                    gp_Dir(1, 0, 0), gp_Dir(0, 1, 0),
                                    "turn " + turns[i].get_name() + " inner crossing");
        if (i + 1 < turns.size()) {
            auto add = turns[i].get_additional_coordinates();
            REQUIRE(add.has_value());
            REQUIRE(!add->empty());
            requireCrossingOnCenterline(conductor->shape,
                                        gp_Pnt((*add)[0][0], (*add)[0][1], 0.0), wr,
                                        gp_Dir(1, 0, 0), gp_Dir(0, 1, 0),
                                        "turn " + turns[i].get_name() + " outer crossing");
        }
    }
    // Flat-wire-on-round-bore sagitta: a rectangular wire's FLAT inner face can't sit flush against
    // the round bore the way a round wire's tangent does -- placed tangent at its centre, its
    // corners dip into the core by ~(height/2)^2 / (2*boreRadius). That is a real, tiny (< 1e-3
    // mm^3) geometric artifact of flat wire on a curved bore, not an interference to fix, so the
    // core<->conductor tolerance here is looser than the round-wire toroid's exact-tangency 1e-12.
    requireNoPairwiseOverlap(named, 1e-9);
}

// Count the solids in a named conductor body.
static int conductorSolidCount(const TopoDS_Shape& shape) {
    int n = 0;
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) ++n;
    return n;
}

TEST_CASE("Real winding: LITZ wire builds ONE continuous body", "[realwinding]") {
    // Litz flows through the round path as a bare bundle. Round column -> ONE single solid.
    auto magneticJson = loadFixture("realwinding_litz_round.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    REQUIRE(conductorSolidCount(conductor->shape) == 1);          // FEM-ready single body
    REQUIRE(BRepCheck_Analyzer(conductor->shape).IsValid());
}

TEST_CASE("Real winding: round-column RECTANGULAR wire is ONE body", "[realwinding]") {
    // 03_buck (PQ32, 3x0.5 mm rectangular wire) -> fixed-binormal single solid.
    auto magneticJson = loadFixture("realwinding_rect_wire_round.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    REQUIRE(conductorSolidCount(conductor->shape) == 1);
    REQUIRE(BRepCheck_Analyzer(conductor->shape).IsValid());
}

TEST_CASE("Real winding: rectangular-column RECTANGULAR wire builds", "[realwinding]") {
    // 18_stacked (E70, 5x1 mm): the per-turn racetrack solids. Its copper turns TOUCH, so the fuse
    // would short them into a brick -- the per-turn compound is kept (correct), so it is multi-solid
    // BY DESIGN. Assert only that it builds valid positive-volume copper for every turn.
    auto magneticJson = loadFixture("realwinding_rect_wire_rect.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);
    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       mvb::DEFAULT_CORE_POLYGON_SEGMENTS, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(shapeVolume(conductor->shape) > 0.0);
    REQUIRE(conductorSolidCount(conductor->shape) >= 1);
}

TEST_CASE("Real winding: MULTI-LAYER spread 3-winding toroidal CMC builds clean",
          "[realwinding]") {
    // A spread 3-winding toroidal CMC with TWO rings per 120-degree section
    // (realwinding_cmc_3w_2layer: 18 turns of 1.4 mm wire per winding). Multilayer
    // toroidal builds once MKF's non-physical outer crossings are corrected on the
    // builder side (MKF ABT #231): each outer ring's outer crossing is re-placed at the
    // physical radial stack (ring 0 outer + ring*OD) AND at the inner crossing's azimuth
    // — MKF staggers outer-ring outer angles out of sequence, which would cross
    // consecutive turns' top chords (and the gate exempts consecutive wraps, so it slips
    // through). Ring returns are depth-staggered under the core. The three windings must
    // be fully independent bodies with two layers each.
    auto magneticJson = loadFixture("realwinding_cmc_3w_2layer.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);
    // Confirm the fixture really is 2 rings per section (else the test is vacuous).
    {
        auto layers = enriched.get_coil().get_layers_description();
        REQUIRE(layers.has_value());
        int primaryConductionLayers = 0;
        for (const auto& L : *layers)
            if (L.get_type() == MAS::ElectricalType::CONDUCTION &&
                !L.get_partial_windings().empty() &&
                L.get_partial_windings()[0].get_winding() == "Primary")
                ++primaryConductionLayers;
        REQUIRE(primaryConductionLayers == 2);
    }

    mvb::MagneticBuilder builder;
    // CMC spread windings stay on the fast drawing compound (femReady=false): the test only asserts
    // the three windings don't overlap EACH OTHER, which the per-run compound already satisfies.
    auto named = builder.buildAllNamed(enriched, false, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/false);
    const char* names[] = {"Primary parallel 0", "Secondary parallel 0",
                           "Tertiary parallel 0"};
    std::vector<const mvb::NamedShape*> conductors;
    for (const char* n : names) conductors.push_back(findConductor(named, n));
    for (size_t i = 0; i < conductors.size(); ++i)
        for (size_t j = i + 1; j < conductors.size(); ++j) {
            double v = commonVolume(conductors[i]->shape, conductors[j]->shape);
            INFO("winding-winding overlap '" << names[i] << "' vs '" << names[j]
                                             << "' = " << v);
            REQUIRE(v <= 1e-12);
        }
    requireNoPairwiseOverlap(named, 1e-12);
}

TEST_CASE("Real winding: single-layer spread 3-winding toroidal CMC builds clean",
          "[realwinding]") {
    auto magneticJson = loadFixture("realwinding_cmc_3w_1layer.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, true);

    mvb::MagneticBuilder builder;
    // CMC spread windings stay on the fast drawing compound (femReady=false): the test only asserts
    // the three windings don't overlap EACH OTHER, which the per-run compound already satisfies.
    auto named = builder.buildAllNamed(enriched, false, 0,
                                       mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/false);

    // One continuous conductor per winding (three windings, each spread over its own
    // ~120-degree arc), all fully independent copper bodies.
    const char* names[] = {"Primary parallel 0", "Secondary parallel 0",
                           "Tertiary parallel 0"};
    std::vector<const mvb::NamedShape*> conductors;
    for (const char* n : names) conductors.push_back(findConductor(named, n));
    for (size_t i = 0; i < conductors.size(); ++i)
        for (size_t j = i + 1; j < conductors.size(); ++j) {
            double v = commonVolume(conductors[i]->shape, conductors[j]->shape);
            INFO("winding-winding overlap '" << names[i] << "' vs '" << names[j]
                                             << "' = " << v);
            REQUIRE(v <= 1e-12);
        }
    requireNoPairwiseOverlap(named, 1e-12);
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

TEST_CASE("Real winding: FEM dense toroid is a conformal (non-overlapping) mitre compound",
          "[realwinding]") {
    // A 60-turn toroid is too dense for the single-solid MakePipe sweep to close on its packed hole
    // spine, so femReady=true builds the CONFORMAL mitre-jointed compound instead: each primitive is
    // its own round solid, and neighbours are sliced on their shared angle-bisector plane so they
    // ABUT on a coincident elliptical face rather than interpenetrating. This is the FEM-meshable
    // (no double material) form of a winding that cannot be a single solid.
    auto magneticJson = loadFixture("realwinding_toroid_3in.json");
    auto enriched = mvb::magnetic_autocomplete_safe(magneticJson, /*useRealWindingGeometry=*/true);

    mvb::MagneticBuilder builder;
    auto named = builder.buildAllNamed(enriched, true, 0, mvb::DEFAULT_WIRE_POLYGON_SEGMENTS,
                                       /*corePolygonSegments=*/0, true, false, false, 0.0,
                                       /*useRealWindingGeometry=*/true, /*femReady=*/true);
    const auto* conductor = findConductor(named, "Primary parallel 0");
    REQUIRE(conductor != nullptr);

    std::vector<TopoDS_Solid> solids;
    for (TopExp_Explorer e(conductor->shape, TopAbs_SOLID); e.More(); e.Next())
        solids.push_back(TopoDS::Solid(e.Current()));
    // MakePipe can't close one body here, so the conductor is a multi-solid conformal compound.
    REQUIRE(solids.size() > 1);
    // Every mitred solid is watertight/valid -- the per-solid ShapeFix makes the torus cut-curves
    // survive a STEP round-trip (they were valid in memory but degraded on reload without it).
    for (const auto& s : solids) REQUIRE(BRepCheck_Analyzer(s).IsValid());

    // Conformal: neighbouring (consecutive-primitive) solids share a mitre face and do NOT overlap.
    // Sample consecutive pairs across the whole winding and require ~zero common volume; a plain
    // overlapping per-run compound would carry ~0.5-1 mm^3 (5e-10 m^3) at each joint here.
    int checked = 0;
    int stride = std::max<int>(1, static_cast<int>(solids.size()) / 15);
    for (size_t i = 0; i + 1 < solids.size(); i += static_cast<size_t>(stride)) {
        double v = commonVolume(solids[i], solids[i + 1]);
        INFO("mitre neighbour overlap [" << i << "," << (i + 1) << "] = " << v);
        REQUIRE(v <= 1e-12);
        ++checked;
    }
    REQUIRE(checked >= 10);
}
