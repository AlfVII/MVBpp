#include "mvb/WireAssembler.h"
#include "mvb/TurnBuilder.h"
#include "mvb/Utils.h"
#include <array>
#include <functional>
#include <set>
#include "constructive_models/Coil.h"
#include "constructive_models/Wire.h"
#include "support/Utils.h"
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRep_Tool.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepCheck_Result.hxx>
#include <BRepCheck_ListOfStatus.hxx>
#include <Geom_BezierCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <cstdio>
#include <BRepAlgoAPI_Cut.hxx>
#include <ShapeFix_ShapeTolerance.hxx>
#include <Precision.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepAlgoAPI_Common.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepOffsetAPI_MakePipe.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepLib.hxx>
#include <GCE2d_MakeSegment.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_Plane.hxx>
#include <Geom_Surface.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <gp_Vec.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Circ.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Dir.hxx>
#include <gp_XY.hxx>
#include <gp_XYZ.hxx>
#include <BRepAlgoAPI_Check.hxx>
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <BOPAlgo_GlueEnum.hxx>
#include <TopTools_ListOfShape.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Wire.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <BRepBuilderAPI_Copy.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <iostream>
#include <Standard_Failure.hxx>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>
#include <numbers>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace mvb {

int curveSampleCount(double radius, double azSpan, double wireRadius) {
    double maxSag = samplingSag(wireRadius);
    double stepAz = (radius > 1e-12)
                        ? 2.0 * std::acos(std::clamp(1.0 - maxSag / radius, 0.0, 1.0))
                        : std::max(azSpan, 1e-3);
    stepAz = std::clamp(stepAz, 1e-3, 0.2);
    return std::max(2, static_cast<int>(std::ceil(azSpan / stepAz)) + 1);
}

// A SPIRAL'S SAMPLING IS NOT A CIRCLE'S. ABT #373 (:818, Alf 2026-08-23). curveSampleCount
// sizes its step from the AZIMUTHAL geometry alone -- radius x azimuth -- which is exactly
// right for an arc or a helix and badly wrong for a spiral whose RADIUS runs while its azimuth
// barely turns. Measured on the FEM toroid's face crossing (realwinding_toroid, wrap
// 'turn 6' -> 'turn 7' bottom chord): 7.98 mm of copper travelled across 0.0377 rad, radius
// 19.978 -> 12.023 mm. The circle rule saw a 0.86 mm arc and returned TWO samples, so the
// "spiral" was emitted as a straight two-point chord -- the very shape faceSpiral exists to
// replace (ABT #685) -- and, worse, its swept end tangents were then the spiral's START
// tangent at BOTH ends while the assembler mitred the joint against the analytic END tangent.
// The junction reported 3.6e-11 deg (tangent, bridged, uncut) and was really a 4.21 deg corner:
// the chord and the following bottom-inner-corner arc interpenetrated by (2/3) r^3 tan(4.21 deg)
// = 0.0051 mm^3 (measured 0.00512 mm^3, 3050/216000 strict-IN grid points).
//
// The step therefore comes from the curve's OWN speed and curvature. For
// P(u) = (cx + r cos u, y, cz - r sin u) with r' = k and y' = m,
//     |P'| = sqrt(k^2 + r^2 + m^2)
//     kappa = sqrt(m^2 (4k^2 + r^2) + (2k^2 + r^2)^2) / (k^2 + r^2 + m^2)^(3/2)
// and a chord of length L sags by L^2 kappa / 8, so the admissible parameter step is
// sqrt(8 sag / kappa) / |P'|. For k = m = 0 that reduces to sqrt(8 sag / r), which is the
// small-angle form of curveSampleCount's own 2 acos(1 - sag/r): arcs and helices are sampled
// EXACTLY as before. The circle rule is kept as a floor so no primitive anywhere can come out
// coarser than it does today; the same [1e-3, 0.2] step clamp applies, for the same reason.
int spiralSampleCount(const Spiral& sp, double wireRadius) {
    const double dAz = sp.az1 - sp.az0;
    const double azSpan = std::abs(dAz);
    int n = curveSampleCount(std::max(sp.r0, sp.r1), azSpan, wireRadius);
    if (azSpan < 1e-12) return n;
    const double maxSag = samplingSag(wireRadius);
    // The blend's r and y advance with f'(t), so k and m are read LOCALLY at each station.
    constexpr int kStations = 9;
    double stepAz = std::numeric_limits<double>::max();
    for (int i = 0; i < kStations; ++i) {
        const double t = static_cast<double>(i) / (kStations - 1);
        const double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
        const double fp = sp.blend ? 0.5 * kPi * std::sin(kPi * t) : 1.0;
        const double r = sp.r0 + (sp.r1 - sp.r0) * f;
        const double k = (sp.r1 - sp.r0) * fp / dAz;
        const double m = (sp.y1 - sp.y0) * fp / dAz;
        const double speed2 = k * k + r * r + m * m;
        if (speed2 < 1e-24) continue;
        const double quad = 2.0 * k * k + r * r;
        const double kappa =
            std::sqrt(m * m * (4.0 * k * k + r * r) + quad * quad) / (speed2 * std::sqrt(speed2));
        if (kappa < 1e-12) continue;   // locally straight: no sag to bound
        stepAz = std::min(stepAz, std::sqrt(8.0 * maxSag / kappa) / std::sqrt(speed2));
    }
    if (stepAz == std::numeric_limits<double>::max()) return n;
    stepAz = std::clamp(stepAz, 1e-3, 0.2);
    return std::max(n, static_cast<int>(std::ceil(azSpan / stepAz)) + 1);
}

std::vector<gp_Pnt> samplePrim(const Primitive& p, double wireRadius) {
    if (p.kind == Primitive::SEG) return {p.seg.a, p.seg.b};
    if (p.kind == Primitive::ARC3) {
        const Arc3& a = p.arc;
        double radius = a.v0.Modulus();
        int n = curveSampleCount(radius, std::abs(a.sweep), wireRadius);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double t = a.sweep * i / (n - 1);
            pts.push_back(gp_Pnt(a.c.XYZ() + rotateXYZ(a.v0, a.axis, t)));
        }
        return pts;
    }
    if (p.kind == Primitive::BLEND) {
        const Blend& bl = p.blendc;
        gp_XYZ d = bl.b.XYZ() - bl.a.XYZ();
        double L = d.Dot(bl.u);
        gp_XYZ perp = d - bl.u * L;
        double perpM = perp.Modulus();
        if (perpM < 1e-12) return {bl.a, bl.b};
        // Peak curvature of the cosine blend: kappa = perp * pi^2 / (2 L^2); sample with
        // the same sag rule as the arcs, on the tightest osculating radius.
        double rMin = 2.0 * L * L / (kPi * kPi * perpM);
        double maxSag = samplingSag(wireRadius);
        double step = 2.0 * std::sqrt(std::max(1e-18, 2.0 * rMin * maxSag));
        double arcLen = std::sqrt(L * L + perpM * perpM) * 1.1;
        int n = std::clamp(static_cast<int>(std::ceil(arcLen / step)) + 1, 8, 512);
        std::vector<gp_Pnt> pts;
        pts.reserve(static_cast<size_t>(n));
        for (int i = 0; i < n; ++i) {
            double t = static_cast<double>(i) / (n - 1);
            pts.push_back(gp_Pnt(bl.a.XYZ() + bl.u * (L * t) +
                                 perp * (0.5 * (1.0 - std::cos(kPi * t)))));
        }
        return pts;
    }
    const Spiral& sp = p.spiral;
    int n = spiralSampleCount(sp, wireRadius);
    std::vector<gp_Pnt> pts;
    pts.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        double t = static_cast<double>(i) / (n - 1);
        double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
        pts.push_back(azPointC(sp.cx, sp.cz, sp.r0 + (sp.r1 - sp.r0) * f,
                               sp.y0 + (sp.y1 - sp.y0) * f, sp.az0 + (sp.az1 - sp.az0) * t));
    }
    return pts;
}

std::pair<gp_Pnt, gp_Pnt> primEndpoints(const Primitive& p) {
    if (p.kind == Primitive::ARC3) {
        return {gp_Pnt(p.arc.c.XYZ() + p.arc.v0),
                gp_Pnt(p.arc.c.XYZ() + rotateXYZ(p.arc.v0, p.arc.axis, p.arc.sweep))};
    }
    if (p.kind == Primitive::SPIRAL) {
        return {azPointC(p.spiral.cx, p.spiral.cz, p.spiral.r0, p.spiral.y0, p.spiral.az0),
                azPointC(p.spiral.cx, p.spiral.cz, p.spiral.r1, p.spiral.y1, p.spiral.az1)};
    }
    if (p.kind == Primitive::BLEND) return {p.blendc.a, p.blendc.b};
    return {p.seg.a, p.seg.b};
}

// --- geometry emission ------------------------------------------------------------------
TopoDS_Wire wireProfileWire(const gp_Pnt& center, const gp_Dir& normal, double radius,
                            int segments) {
    // DETERMINISTIC section orientation. gp_Ax2(center, normal) lets OCC choose the X
    // direction, so two primitives meeting with the SAME tangent could still get n-gons rotated
    // relative to each other: their abutting faces then only partially coincide, and at a mitre
    // junction the pieces can miss entirely (measured: a 56-primitive rect-column conductor fell
    // into 7 components with no endpoint gaps and nothing dropped). Deriving X from a fixed world
    // reference makes the orientation a pure function of the tangent, so tangent neighbours share
    // an exactly coincident section -- which is what conformal abutment requires.
    gp_XYZ ref(0, 1, 0);
    if (std::abs(normal.XYZ().Dot(ref)) > 0.9) ref = gp_XYZ(1, 0, 0);
    gp_XYZ xdir = ref - normal.XYZ() * normal.XYZ().Dot(ref);
    gp_Ax2 plane(center, normal, gp_Dir(xdir));
    if (segments <= 0) {
        gp_Circ circ(plane, radius);
        return BRepBuilderAPI_MakeWire(BRepBuilderAPI_MakeEdge(circ).Edge()).Wire();
    }
    BRepBuilderAPI_MakePolygon poly;
    gp_Dir dx = plane.XDirection();
    gp_Dir dy = plane.YDirection();
    double offset = kPi / segments;
    // INSCRIBED n-gon (vertices on the nominal circle): the faceted wire never protrudes past
    // the wire's real envelope, which is what the zero-tolerance pairwise-overlap contracts rely
    // on. Its flats consequently sit at the apothem, r*cos(pi/n) -- 0.981*r at n=16 -- so any
    // clearance or coverage check against a faceted section must use the apothem, not r.
    for (int i = 0; i < segments; ++i) {
        double ang = kTwoPi * i / segments + offset;
        gp_XYZ p = center.XYZ() + dx.XYZ() * (radius * std::cos(ang)) +
                   dy.XYZ() * (radius * std::sin(ang));
        poly.Add(gp_Pnt(p));
    }
    poly.Close();
    return poly.Wire();
}

TopoDS_Face wireProfile(const gp_Pnt& center, const gp_Dir& normal, double radius,
                        int segments) {
    return BRepBuilderAPI_MakeFace(wireProfileWire(center, normal, radius, segments)).Face();
}

// The same round profile, but as TWO half-circle edges instead of one closed circle.
// Sweeping a CLOSED circular edge produces a single PERIODIC B-spline surface carrying a seam
// curve, and gmsh cannot order the boundary loop of such a face -- it aborts the 2D stage with
// "the 1D mesh seems not to be forming a closed loop" (the offending face reported simply as
// "BSpline surface"). That is what made the real-winding conductor unmeshable: every other
// primitive is an analytic cylinder or torus, which mesh fine; only the swept axial transitions
// carried a seam. Split in two and the sweep yields ordinary patches with real boundaries -- the
// identical solid, now meshable. (ABT #332)
TopoDS_Wire wireProfileWireSplit(const gp_Pnt& center, const gp_Dir& normal, double radius,
                                 int pieces = 2) {
    gp_Circ circ(gp_Ax2(center, normal), radius);
    BRepBuilderAPI_MakeWire w;
    for (int i = 0; i < pieces; ++i)
        w.Add(BRepBuilderAPI_MakeEdge(circ, kTwoPi * i / pieces, kTwoPi * (i + 1) / pieces).Edge());
    return w.Wire();
}

// Oriented RECTANGLE profile for rectangular/planar wire, in the plane perpendicular to the wire's
// travel (tangent). The AXIAL in-plane axis is the column axis projected perpendicular to the
// tangent; the RADIAL axis is tangent x axial. width spans radial, height spans axial -- matching
// TurnBuilder::build_rect_profile so the real-winding section is consistent with the per-turn one.
TopoDS_Wire rectProfileWire(const gp_Pnt& center, const gp_Dir& tangent, const gp_Dir& axialAxis,
                            double width, double height) {
    gp_XYZ t = tangent.XYZ();
    gp_XYZ a = axialAxis.XYZ() - t * (axialAxis.XYZ().Dot(t));  // axial, ⊥ tangent
    if (a.Modulus() < 1e-9) {                                   // tangent ∥ axial: pick any ⊥ axis
        gp_XYZ ref = std::abs(t.X()) < 0.9 ? gp_XYZ(1, 0, 0) : gp_XYZ(0, 0, 1);
        a = ref - t * (ref.Dot(t));
    }
    a.Normalize();
    gp_XYZ w = t.Crossed(a);                                    // radial
    w.Normalize();
    double hw = width / 2.0, hh = height / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(center.XYZ() - w * hw - a * hh));
    poly.Add(gp_Pnt(center.XYZ() + w * hw - a * hh));
    poly.Add(gp_Pnt(center.XYZ() + w * hw + a * hh));
    poly.Add(gp_Pnt(center.XYZ() - w * hw + a * hh));
    poly.Close();
    return poly.Wire();
}

// One analytic edge per primitive, so the whole conductor is ONE path:
//  - SEG: straight edge.
//  - ARC3: exact circle edge trimmed to [0, sweep].
//  - SPIRAL at constant radius: TRUE HELIX — a 2D line in the (U = azimuth, V = height)
//    parametric space of a cylinder about its vertical axis (the canonical OCCT
//    construction: Geom2d line on Geom_CylindricalSurface + BRepLib::BuildCurves3d).
//  - SPIRAL with varying radius (conical wraps, oblong cap ramps): BSpline through the
//    sampled centreline.
// A spiral drawn on an analytic surface can absorb its mitre growth INTO THE PCURVE, so the
// grown spine stays ONE edge on ONE surface. ABT #685 (Alf, 2026-08-18): the alternative --
// splicing a straight tangent edge onto each end, which is what every other primitive does --
// makes MakePipeShell emit a lateral face PER SPINE EDGE, and the mitre plane then falls exactly
// on the seam between two of them (the plane passes through the primitive's own endpoint, which
// is precisely where the splice sits). OCC's boolean cannot cut a pipe on a face boundary: every
// failure measured had faces=5 (3 lateral + 2 caps) and returned debris -- 06_llc lost 3.0 of
// 3.87 mm3 to a 0.59 mm3 knife. With the growth inside the pcurve there is no seam to cut on.
bool primEdgeGrowsAnalytically(const Primitive& pr) {
    if (pr.kind != Primitive::SPIRAL) return false;
    const Spiral& sp = pr.spiral;
    return std::abs(sp.r1 - sp.r0) < 1e-12 || std::abs(sp.y1 - sp.y0) > 1e-12;
}

TopoDS_Edge primEdge(const Primitive& pr, double wireRadius, double overA, double overB) {
    if (pr.kind == Primitive::SEG) {
        if (pr.seg.a.Distance(pr.seg.b) < 1e-12) return TopoDS_Edge();
        return BRepBuilderAPI_MakeEdge(pr.seg.a, pr.seg.b).Edge();
    }
    if (pr.kind == Primitive::ARC3) {
        double radius = pr.arc.v0.Modulus();
        if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return TopoDS_Edge();
        try {
            // gp_Circ parametrisation P(u) = c + r (cos u XDir + sin u YDir) with
            // YDir = axis x XDir matches rotateXYZ (axis is perpendicular to v0).
            gp_Ax2 ax2(pr.arc.c, gp_Dir(pr.arc.axis), gp_Dir(pr.arc.v0));
            gp_Circ circ(ax2, radius);
            return BRepBuilderAPI_MakeEdge(circ, 0.0, pr.arc.sweep).Edge();
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    if (pr.kind == Primitive::SPIRAL) {
        // Every spiral lives on an ANALYTIC surface of revolution about its vertical
        // axis: a cylinder when the radius is constant, a cone otherwise (radius and
        // height share the blend factor, so the (r, y) trace is a straight meridian).
        // Building the edge as a 2D curve in the surface's (U = azimuth, V) parameter
        // space (+ BRepLib::BuildCurves3d) is the one construction OCCT's pipe-shell
        // sweeps reliably — free-form interpolated 3D splines, though positionally
        // tight, produce unbounded surface flares when swept.
        const Spiral& sp = pr.spiral;
        double dr = sp.r1 - sp.r0;
        double dy = sp.y1 - sp.y0;
        bool constantRadius = std::abs(dr) < 1e-12;
        if (constantRadius || std::abs(dy) > 1e-12) {
            try {
                // Frame with U measured from +X toward -Z — exactly the azPointC()
                // convention: P(U,V) matches (cx + r cos U, y, cz - r sin U).
                gp_Ax3 axis(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                            gp_Dir(1, 0, 0));
                Handle(Geom_Surface) surf;
                double v1;   // V at the end (V = 0 at the start)
                if (constantRadius) {
                    axis.SetLocation(gp_Pnt(sp.cx, 0, sp.cz));
                    surf = new Geom_CylindricalSurface(axis, sp.r0);
                    // Cylinder V is plain height.
                    v1 = dy;
                } else {
                    double alpha = std::atan(dr / dy);
                    surf = new Geom_ConicalSurface(axis, alpha, sp.r0);
                    v1 = dy / std::cos(alpha);
                }
                double v0 = constantRadius ? sp.y0 : 0.0;
                // Mitre growth, in arc length, converted to azimuth on this surface: a step du
                // advances sqrt(radius^2 + (dv/du)^2) of arc. Extending the pcurve's range is an
                // EXACT continuation of the same helix/conical spiral for the linear case, and a
                // tangent (constant-V) continuation for the blend, whose end tangents are purely
                // azimuthal by construction. Either way it stays one edge -- see the note above.
                const double dAz = sp.az1 - sp.az0;
                const double sgn = (dAz >= 0.0) ? 1.0 : -1.0;
                const double dvdu = (std::abs(dAz) > 1e-12) ? v1 / dAz : 0.0;
                const double dsduA = std::sqrt(sp.r0 * sp.r0 + dvdu * dvdu);
                const double dsduB = std::sqrt(sp.r1 * sp.r1 + dvdu * dvdu);
                const double az0e =
                    sp.az0 - (overA > 0.0 && dsduA > 1e-12 ? sgn * overA / dsduA : 0.0);
                const double az1e =
                    sp.az1 + (overB > 0.0 && dsduB > 1e-12 ? sgn * overB / dsduB : 0.0);
                TopoDS_Edge e;
                if (!sp.blend) {
                    auto vAt = [&](double az) { return v0 + (az - sp.az0) * dvdu; };
                    Handle(Geom2d_TrimmedCurve) seg2d =
                        GCE2d_MakeSegment(gp_Pnt2d(az0e, vAt(az0e)),
                                          gp_Pnt2d(az1e, vAt(az1e))).Value();
                    e = BRepBuilderAPI_MakeEdge(seg2d, surf).Edge();
                } else {
                    // Cosine blend of V over U — end tangents purely azimuthal, so the
                    // junctions into the neighbouring on-station geometry are exact.
                    int n = curveSampleCount(std::max(sp.r0, sp.r1),
                                             std::abs(az1e - az0e), wireRadius);
                    Handle(TColgp_HArray1OfPnt2d) arr =
                        new TColgp_HArray1OfPnt2d(1, n);
                    for (int i = 0; i < n; ++i) {
                        double az = az0e + (az1e - az0e) * (static_cast<double>(i) / (n - 1));
                        // Clamped: outside the primitive's own range the blend holds its end
                        // value, which IS the tangent continuation (the cosine's slope is zero
                        // at both ends), so the imposed (1,0) end tangents stay exact.
                        const double t = std::clamp(
                            std::abs(dAz) > 1e-12 ? (az - sp.az0) / dAz : 0.0, 0.0, 1.0);
                        double f = 0.5 * (1.0 - std::cos(kPi * t));
                        arr->SetValue(i + 1, gp_Pnt2d(az, v0 + v1 * f));
                    }
                    Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
                    interp.Load(gp_Vec2d(1, 0), gp_Vec2d(1, 0), Standard_True);
                    interp.Perform();
                    if (!interp.IsDone() || interp.Curve().IsNull())
                        return TopoDS_Edge();
                    e = BRepBuilderAPI_MakeEdge(interp.Curve(), surf).Edge();
                }
                BRepLib::BuildCurves3d(e);
                return e;
            } catch (const Standard_Failure&) {
                return TopoDS_Edge();
            }
        }
        // FLAT SPIRAL (radius varies at exactly constant height) -- blended OR linear.
        // A 2D curve in the horizontal plane's parameter space, interpolated with the EXACT
        // analytic end tangents.
        //
        // ABT #373 (:818): the LINEAR case used to fall through to a 3D BSpline through
        // samplePrim's points, and on a toroid's face crossing samplePrim returned TWO of them
        // (see spiralSampleCount) -- so the emitted piece was a straight chord swept by
        // MakePipeShell along a straight spine, which carries the START tangent's section to
        // BOTH ends. The assembler, meanwhile, mitres against spiralTangent's analytic END
        // tangent. On realwinding_toroid's 'turn 6' -> 'turn 7' bottom chord those differ by the
        // spiral's whole turning, 4.21 deg: the joint was read as tangent (3.6e-11 deg, bridged,
        // uncut) and the chord and the bottom-inner-corner arc interpenetrated by 0.00512 mm^3.
        // Interpolating in the plane with the analytic tangents IMPOSED makes the swept end
        // tangent equal the one the junction logic uses, by construction -- which is what
        // "neighbours meet on the plane bisecting their tangents" requires.
        {
            const double dAz = sp.az1 - sp.az0;
            // No azimuthal advance at all: the locus IS the straight radial segment between the
            // endpoints (r linear, y and az constant), and that is the exact edge for it --
            // faceSpiralTangent reports the same purely radial direction for this case.
            if (std::abs(dAz) < 1e-12) {
                const auto ends = primEndpoints(pr);
                if (ends.first.Distance(ends.second) < 1e-12) return TopoDS_Edge();
                return BRepBuilderAPI_MakeEdge(ends.first, ends.second).Edge();
            }
            try {
                gp_Ax3 frame(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                             gp_Dir(1, 0, 0));
                Handle(Geom_Plane) plane = new Geom_Plane(frame);
                // Plane P(U,V) = O + U XDir + V YDir with YDir = -Z, so the azPointC
                // trace maps to (U, V) = (r cos az, r sin az).
                int n = spiralSampleCount(sp, wireRadius);
                Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
                for (int i = 0; i < n; ++i) {
                    double t = static_cast<double>(i) / (n - 1);
                    double f = sp.blend ? 0.5 * (1.0 - std::cos(kPi * t)) : t;
                    double az = sp.az0 + dAz * t;
                    double r = sp.r0 + dr * f;
                    arr->SetValue(i + 1, gp_Pnt2d(r * std::cos(az), r * std::sin(az)));
                }
                // d/dt (r cos az, r sin az) with az' = dAz and r' = dr f'(t). A cosine blend has
                // f'(0) = f'(1) = 0, so this reduces to dAz * r * (-sin az, cos az) -- the purely
                // azimuthal tangent the blend has always been given, now carrying dAz's SIGN as
                // well (a spiral running the other way round was previously handed a reversed
                // end tangent).
                auto tangentAt = [&](double az, double r, double rp) {
                    return gp_Vec2d(rp * std::cos(az) - r * dAz * std::sin(az),
                                    rp * std::sin(az) + r * dAz * std::cos(az));
                };
                const double rp = sp.blend ? 0.0 : dr;
                gp_Vec2d ta = tangentAt(sp.az0, sp.r0, rp);
                gp_Vec2d tb = tangentAt(sp.az1, sp.r1, rp);
                Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
                interp.Load(ta, tb, Standard_True);
                interp.Perform();
                if (!interp.IsDone() || interp.Curve().IsNull()) return TopoDS_Edge();
                TopoDS_Edge e = BRepBuilderAPI_MakeEdge(interp.Curve(), plane).Edge();
                BRepLib::BuildCurves3d(e);
                return e;
            } catch (const Standard_Failure&) {
                return TopoDS_Edge();
            }
        }
    }
    if (pr.kind == Primitive::BLEND) {
        // The S-blend is PLANAR (the longitudinal direction and the offset span one
        // plane), so it too builds as a 2D cosine in an analytic surface's parameter
        // space — the pipe-reliable construction (see the spiral comment above).
        const Blend& bl = pr.blendc;
        gp_XYZ d = bl.b.XYZ() - bl.a.XYZ();
        double L = d.Dot(bl.u);
        gp_XYZ perp = d - bl.u * L;
        double perpM = perp.Modulus();
        if (L < 1e-12) return TopoDS_Edge();
        if (perpM < 1e-12) return BRepBuilderAPI_MakeEdge(bl.a, bl.b).Edge();
        try {
            gp_Ax3 frame(bl.a, gp_Dir(bl.u.Crossed(perp)), gp_Dir(bl.u));
            Handle(Geom_Plane) plane = new Geom_Plane(frame);
            // Plane P(U,V) = O + U*XDir + V*YDir with XDir = u, YDir = axis x XDir =
            // unit(perp).
            double stepL = 2.0 * std::sqrt(std::max(
                1e-18, 4.0 * L * L / (kPi * kPi * perpM) *
                           samplingSag(wireRadius)));
            int n = std::clamp(static_cast<int>(std::ceil(L / stepL)) + 1, 8, 512);
            Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
            for (int i = 0; i < n; ++i) {
                double t = static_cast<double>(i) / (n - 1);
                arr->SetValue(i + 1, gp_Pnt2d(L * t,
                                              perpM * 0.5 * (1.0 - std::cos(kPi * t))));
            }
            Geom2dAPI_Interpolate interp(arr, Standard_False, 1e-12);
            interp.Load(gp_Vec2d(1, 0), gp_Vec2d(1, 0), Standard_True);
            interp.Perform();
            if (!interp.IsDone() || interp.Curve().IsNull()) return TopoDS_Edge();
            TopoDS_Edge e = BRepBuilderAPI_MakeEdge(interp.Curve(), plane).Edge();
            BRepLib::BuildCurves3d(e);
            return e;
        } catch (const Standard_Failure&) {
            return TopoDS_Edge();
        }
    }
    // Every kind returns above: SEG, ARC3, and all three spiral families (on a cylinder, on a
    // cone, and flat in a plane), plus BLEND. ABT #373: the flat LINEAR spiral used to land here
    // on a 3D BSpline through samplePrim's points -- the construction whose two-point degenerate
    // case produced the :818 mitre interpenetration -- and it now has its own analytic branch.
    // Nothing may fall through to a silent empty edge; a new Primitive kind must bring its own
    // construction.
    throw std::runtime_error(
        "WireAssembler: primEdge has no construction for primitive '" + pr.label +
        "' (unknown kind " + std::to_string(static_cast<int>(pr.kind)) + ")");
}

// Fuse all solids of a conductor's compound into as FEW bodies as possible, in ONE general
// boolean (BRepAlgoAPI_Fuse over the whole set via SetArguments/SetTools). A single BOP with
// OBB culling is both far faster than N sequential pairwise fuses (which is O(N) booleans on
// swept BSpline solids — minutes for a many-piece winding) AND more robust: OCCT resolves all
// the mutual intersections together, so overlapping chunk pipes + elbow/junction spheres merge
// into one solid while any piece that genuinely touches nothing stays its own solid. No glue:
// the pieces OVERLAP volumetrically (junction spheres straddle the pipe caps), which is a real
// intersection, not the coincident-face case glue is for. A gross-corruption volume guard
// rejects the known OCCT "ate a whole body" defect; on any failure the exact (unfused) compound
// is returned rather than losing copper — same guarantee as before, just reached far less often.
// A conductor face whose bounding extent is many wire-radii long but whose area is a sliver:
// the degenerate SEAM faces OCCT's boolean union leaves where two swept pipes (or a pipe and
// its junction sphere) merge -- the thin "line/cone" sheets seen on a multilayer winding. A
// genuine tube-wall face of the same extent has area ~ 2*pi*r*length, orders of magnitude
// larger; end caps and short faces (diag <= a few radii) are exempt.
bool hasDegenerateSheetFace(const TopoDS_Shape& shape, double wireRadius) {
    for (TopExp_Explorer fe(shape, TopAbs_FACE); fe.More(); fe.Next()) {
        Bnd_Box fb;
        BRepBndLib::Add(fe.Current(), fb);
        if (fb.IsVoid()) continue;
        double x0, y0, z0, x1, y1, z1;
        fb.Get(x0, y0, z0, x1, y1, z1);
        double diag = std::sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) +
                                (z1 - z0) * (z1 - z0));
        if (diag <= 4.0 * wireRadius) continue;
        GProp_GProps fp;
        BRepGProp::SurfaceProperties(fe.Current(), fp);
        if (fp.Mass() < 0.2 * wireRadius * diag) return true;
    }
    return false;
}

// Drop degenerate (near-zero-volume) solids from a conductor result. OCCT's fuse can leave
// thin sheet/shell solids behind when it welds a dense multilayer winding: they carry no
// copper, render as stray sheets, and break tetrahedral meshing. Any real conductor piece --
// even a lead-junction sphere -- is orders of magnitude above the 1e-12 m^3 (1e-3 mm^3) floor.
TopoDS_Shape pruneDegenerateSolids(const TopoDS_Shape& shape, double wireRadius) {
    // The threshold is RELATIVE TO THE WIRE, never an absolute volume. This exists to drop
    // boolean SLIVERS -- sheets with no thickness -- so the scale that matters is the wire's
    // own: a piece thinner than a hundredth of a wire radius along the wire is a sliver, any
    // longer piece is real copper. An absolute 1e-12 m^3 cutoff silently deleted every bump
    // riser of 13_current_sense's 0.049 mm-radius secondary (riser volume 8.4e-13 m^3, a
    // legitimate 0.112 mm step), leaving the bump arcs disconnected -- Alf, 2026-08-08:
    // "Secondary parallel 0 16 and 17 are not connected, they are missing the straight
    // segment of the bump".
    const double kMinSolidVolume =
        kPi * wireRadius * wireRadius * (wireRadius / 100.0);
    std::vector<TopoDS_Shape> kept;
    int total = 0;
    double droppedVolume = 0.0;
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        ++total;
        GProp_GProps props;
        BRepGProp::VolumeProperties(exp.Current(), props);
        if (std::abs(props.Mass()) >= kMinSolidVolume) kept.push_back(exp.Current());
        else droppedVolume += std::abs(props.Mass());
    }
    // Dropping copper is never routine: say so.
    if (static_cast<int>(kept.size()) != total)
        std::cerr << "[ConductorBuilder] pruneDegenerateSolids dropped "
                  << (total - static_cast<int>(kept.size())) << " of " << total
                  << " solids as slivers (" << droppedVolume * 1e9 << " mm3, threshold "
                  << kMinSolidVolume * 1e9 << " mm3)\n";
    if (kept.empty() || static_cast<int>(kept.size()) == total) return shape;
    if (kept.size() == 1) return kept.front();
    TopoDS_Compound out;
    BRep_Builder b;
    b.MakeCompound(out);
    for (const auto& s : kept) b.Add(out, s);
    return out;
}

// ---------------------------------------------------------------------------------------
// Conformal per-primitive toroid (mitre joints). A dense toroid's winding cannot become one
// solid (MakePipe self-intersects on the packed hole spine; a boolean fuse of the per-run pipes
// comes out invalid), and the plain per-primitive compound OVERLAPS at every joint (double
// material -> not a valid FEM assembly). Here each primitive is built as its OWN round solid,
// grown a hair past both endpoints, then sliced by the ANGLE-BISECTOR plane through each shared
// point. Two equal round tubes mitred on their bisector leave IDENTICAL elliptical faces (same
// centre, tilt and radius), so neighbouring solids ABUT on a coincident face instead of
// overlapping: a conformal, mesh-shareable assembly. Centrelines/crossings never move.
// ANALYTIC end tangents, not sampled chords. Chord directions from samplePrim run up to half
// a sample step's turn away from the true tangent (a 16-sample wrap chord is ~5-6 deg off),
// which made emitToroidConformal misread every TANGENT wrap junction as a corner and slice it
// on a chord-tilted bisector: the mitre faces then no longer coincide and the 'conformal'
// assembly falls apart into disconnected pieces (measured: PQ33 round column, 18 prims,
// 34 spurious cuts, 14 components).
static gp_Dir spiralTangent(const Spiral& sp, bool atStart) {
    const double daz = sp.az1 - sp.az0;
    const double az = atStart ? sp.az0 : sp.az1;
    double rp = 0.0, yp = 0.0;
    const double r = atStart ? sp.r0 : sp.r1;
    if (std::fabs(daz) > 1e-12 && !sp.blend) {
        rp = (sp.r1 - sp.r0) / daz;
        yp = (sp.y1 - sp.y0) / daz;
    }   // cosine-blend spirals have purely azimuthal end tangents (r' = y' = 0)
    // P(az) = (cx + r cos az, y, cz - r sin az)  [azPointC convention]
    gp_XYZ t(rp * std::cos(az) - r * std::sin(az), yp,
             -rp * std::sin(az) - r * std::cos(az));
    if (daz < 0) t *= -1.0;
    if (t.Modulus() < 1e-15) return gp_Dir(1, 0, 0);
    return gp_Dir(t);
}

gp_Dir primFwdStart(const Primitive& p, double r) {
    switch (p.kind) {
        case Primitive::SEG: {
            gp_Vec v(p.seg.a, p.seg.b);
            if (v.Magnitude() > 1e-12) return gp_Dir(v);
            break;
        }
        case Primitive::ARC3: {
            gp_XYZ t = p.arc.axis.Crossed(p.arc.v0);
            if (p.arc.sweep < 0) t *= -1.0;
            if (t.Modulus() > 1e-15) return gp_Dir(t);
            break;
        }
        case Primitive::SPIRAL:
            return spiralTangent(p.spiral, /*atStart=*/true);
        case Primitive::BLEND:
            if (p.blendc.u.Modulus() > 1e-15) return gp_Dir(p.blendc.u);
            break;
    }
    auto pts = samplePrim(p, r);
    for (size_t i = 1; i < pts.size(); ++i) {
        gp_Vec v(pts.front(), pts[i]);
        if (v.Magnitude() > 1e-12) return gp_Dir(v);
    }
    return gp_Dir(1, 0, 0);
}

gp_Dir primFwdEnd(const Primitive& p, double r) {
    switch (p.kind) {
        case Primitive::SEG: {
            gp_Vec v(p.seg.a, p.seg.b);
            if (v.Magnitude() > 1e-12) return gp_Dir(v);
            break;
        }
        case Primitive::ARC3: {
            gp_XYZ ve = rotateXYZ(p.arc.v0, p.arc.axis, p.arc.sweep);
            gp_XYZ t = p.arc.axis.Crossed(ve);
            if (p.arc.sweep < 0) t *= -1.0;
            if (t.Modulus() > 1e-15) return gp_Dir(t);
            break;
        }
        case Primitive::SPIRAL:
            return spiralTangent(p.spiral, /*atStart=*/false);
        case Primitive::BLEND:
            if (p.blendc.u.Modulus() > 1e-15) return gp_Dir(p.blendc.u);
            break;
    }
    auto pts = samplePrim(p, r);
    for (size_t i = pts.size(); i-- > 1;) {
        gp_Vec v(pts[i - 1], pts.back());
        if (v.Magnitude() > 1e-12) return gp_Dir(v);
    }
    return gp_Dir(1, 0, 0);
}

// A round solid for one primitive, grown by overA/overB past its start/end so a tilted bisector
// plane fully crosses the tube there (SEG -> longer cylinder, ARC3 -> wider revolve, else pipe).
// A tangent junction passes 0 -> the tube ends flush on a perpendicular cap, no growth, no cut.
// segments <= 0 keeps the EXACT round section (analytic cylinders/tori); segments > 0 emits the
// n-gon prism/revolve. This is not cosmetic: an exact round sweep produces PERIODIC surfaces
// (seam-carrying cylinders and tori), and gmsh's "Impossible to mesh periodic surface" is a
// known failure cluster on them -- the faceted section has no seam at all, which is why the
// pipeline has always meshed with segments=12. The conformal assembler must honour the same
// knob the retired sweep path did.
static TopoDS_Shape rawGrownSolid(const Primitive& pr, double r, double overA, double overB,
                                  int segments) {
    if (pr.kind == Primitive::SEG) {
        gp_Vec d(pr.seg.a, pr.seg.b);
        double len = d.Magnitude();
        if (len < 1e-12) return {};
        gp_Dir dir(d);
        gp_Pnt a = pr.seg.a.Translated(gp_Vec(dir) * (-overA));
        const double total = len + overA + overB;
        if (segments <= 0)
            return BRepPrimAPI_MakeCylinder(gp_Ax2(a, dir), r, total).Shape();
        TopoDS_Face prof = wireProfile(a, dir, r, segments);
        return BRepPrimAPI_MakePrism(prof, gp_Vec(dir) * total).Shape();
    }
    if (pr.kind == Primitive::ARC3) {
        double radius = pr.arc.v0.Modulus();
        if (radius < 1e-12 || std::abs(pr.arc.sweep) < 1e-12) return {};
        double sgn = pr.arc.sweep < 0 ? -1.0 : 1.0;
        double ddA = (overA / radius) * sgn, ddB = (overB / radius) * sgn;
        double total = pr.arc.sweep + ddA + ddB;
        if (std::getenv("MVB_GROW_DIAG") && (overA > 0 || overB > 0)) {
            gp_XYZ vs = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA);
            gp_XYZ ve = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA + total);
            std::fprintf(stderr,
                "[grow-arc] '%s' overA=%.4f overB=%.4f grownStart=(%.4f,%.4f,%.4f) grownEnd=(%.4f,%.4f,%.4f)\n",
                pr.label.c_str(), overA * 1e3, overB * 1e3,
                (pr.arc.c.X() + vs.X()) * 1e3, (pr.arc.c.Y() + vs.Y()) * 1e3,
                (pr.arc.c.Z() + vs.Z()) * 1e3, (pr.arc.c.X() + ve.X()) * 1e3,
                (pr.arc.c.Y() + ve.Y()) * 1e3, (pr.arc.c.Z() + ve.Z()) * 1e3);
        }
        if (std::abs(total) > 1.9 * kPi) return {};  // a near-full revolve would self-close: caller retries
        gp_XYZ v0e = rotateXYZ(pr.arc.v0, pr.arc.axis, -ddA);
        gp_Pnt start(pr.arc.c.XYZ() + v0e);
        gp_XYZ tangent = pr.arc.axis.Crossed(v0e);
        TopoDS_Face prof = wireProfile(start, gp_Dir(tangent), r, segments);
        BRepPrimAPI_MakeRevol rev(prof, gp_Ax1(pr.arc.c, gp_Dir(pr.arc.axis)), total);
        return rev.IsDone() ? rev.Shape() : TopoDS_Shape();
    }
    // BLEND / SPIRAL: swept as exact pipes. A pipe has no revolve/extrude axis to overhang
    // along, so mitre growth EXTENDS THE SPINE instead: a short straight run along the analytic
    // end tangent (G1-continuous, so MakePipeShell sweeps it as one smooth body) which the
    // bisector cut then trims back -- identical overhang semantics to the SEG/ARC3 branches.
    // Without this, a corner cut at a spiral/blend end carved a wedge with nothing to fill it
    // (measured: every lead<->wrap corner on PQ33 left a gap -> 5 disconnected components).
    auto pts = samplePrim(pr, r);
    if (pts.size() < 2) return {};
    // Where the growth can live INSIDE the curve (spirals on an analytic surface), put it there:
    // the spine is then a single edge and the pipe a single lateral face, with no seam for the
    // mitre plane to land on. See primEdgeGrowsAnalytically.
    const bool analyticGrowth = primEdgeGrowsAnalytically(pr);
    TopoDS_Edge e = analyticGrowth ? primEdge(pr, r, overA, overB) : primEdge(pr, r);
    if (e.IsNull()) return {};
    try {
        const auto ends = primEndpoints(pr);
        const gp_Dir tA = primFwdStart(pr, r), tB = primFwdEnd(pr, r);
        BRepBuilderAPI_MakeWire mw;
        gp_Pnt spineStart = ends.first;
        gp_Dir spineDir = tA;
        if (analyticGrowth) {
            mw.Add(e);
            // The profile must sit square to the GROWN spine's start, which the extended pcurve
            // has moved; read both straight off the edge rather than re-deriving them.
            BRepAdaptor_Curve bac(e);
            gp_Pnt p0;
            gp_Vec d0;
            bac.D1(bac.FirstParameter(), p0, d0);
            if (d0.Magnitude() > 1e-12) {
                spineStart = p0;
                spineDir = gp_Dir(d0);
                if (spineDir.Dot(tA) < 0.0) spineDir.Reverse();
            }
        } else {
            if (overA > 0.0) {
                spineStart = gp_Pnt(ends.first.XYZ() - tA.XYZ() * overA);
                mw.Add(BRepBuilderAPI_MakeEdge(spineStart, ends.first).Edge());
            }
            mw.Add(e);
            if (overB > 0.0)
                mw.Add(BRepBuilderAPI_MakeEdge(
                           ends.second, gp_Pnt(ends.second.XYZ() + tB.XYZ() * overB)).Edge());
        }
        if (!mw.IsDone()) return {};
        TopoDS_Wire spine = mw.Wire();
        // ONE attempt, at the REQUESTED faceting -- no round-profile retry (Alf: no
        // fallbacks). A failed faceted sweep returns null and the caller throws naming the
        // primitive, so the geometry that produced it gets fixed instead of silently mixing
        // a round (periodic-surface, gmsh-hostile) piece into a faceted conductor.
        //
        // The sweep runs in a MILLIMETRE frame: OCC's MakePipeShell mis-frames
        // sub-millimetre faceted profiles in METRES (measured standalone: the identical
        // helix/profile that fails at x1 -- PipeNotDone on a 0.1 mm 12-gon, 10_emi/
        // 13_current_sense/20_iso/24_margin -- sweeps cleanly at x1000). Pure scaling is an
        // EXACT affine map both ways, so this is a numerical frame choice, not an
        // approximation. TODO: migrate the other MakePipeShell sites to the same frame as
        // they are touched.
        TopoDS_Wire prof = wireProfileWire(spineStart, spineDir, r, segments);
        gp_Trsf up, down;
        up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
        down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
        TopoDS_Wire spineMm =
            TopoDS::Wire(BRepBuilderAPI_Transform(spine, up, Standard_True).Shape());
        TopoDS_Wire profMm =
            TopoDS::Wire(BRepBuilderAPI_Transform(prof, up, Standard_True).Shape());
        BRepOffsetAPI_MakePipeShell ps(spineMm);
        ps.Add(profMm);
        ps.Build();
        if (ps.IsDone() && ps.MakeSolid())
            return BRepBuilderAPI_Transform(ps.Shape(), down, Standard_True).Shape();
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[sweep-fail] MakePipeShell status=" << (int)ps.GetStatus()
                      << " isDone=" << ps.IsDone() << " segments=" << segments
                      << " r=" << r << "\n";
    } catch (const Standard_Failure& f) {
        if (std::getenv("MVB_DIAG"))
            std::cerr << "[sweep-fail] OCC exception: "
                      << (f.GetMessageString() ? f.GetMessageString() : "(null)") << "\n";
    }
    return {};
}

// SCALING A SHAPE SCALES ITS TOLERANCE TOO -- AND A 1000x TOLERANCE BLINDS THE BOOLEAN.
// ABT #860 (2026-08-23). Every mitre boolean runs in the MILLIMETRE frame because OCC mis-frames
// sub-millimetre profiles in metres. BRepBuilderAPI_Transform carries the shape's tolerances
// through the same scale, so a solid holding OCC's own confusion (1e-7) in metres arrives in the
// millimetre frame claiming 1e-4 -- a tenth of a micron of slop on geometry that is analytically
// exact. At that tolerance OCC 7.9.3 declares a TOROIDAL face and the knife's PLANE
// non-intersecting: Cut returns its argument UNCHANGED and Common returns EMPTY, both reporting
// IsDone with no error and no warning. Measured on 24_margin's entrance corner (--segments 0):
// 0.0581 mm3 of the terminal stub stands inside the knife box by point classification, and both
// booleans see nothing; with the tolerance reset to 1e-7 the same Cut removes 0.0583 mm3. That is
// the whole '--segments 0 mitred corner did not abut' failure -- the faceted profile never hit it
// because a plane/plane intersection survives the inflated tolerance that a plane/torus one does
// not.
// A shape does not become less accurate by being measured in millimetres, so the millimetre copy
// is given the millimetre frame's confusion -- and only if it is still VALID with it, because a
// swept BSpline pipe may genuinely need the slop it carries.
static const gp_Trsf& mitreUpTrsf() {
    static const gp_Trsf t = [] { gp_Trsf x; x.SetScale(gp_Pnt(0, 0, 0), 1000.0); return x; }();
    return t;
}
static const gp_Trsf& mitreDownTrsf() {
    static const gp_Trsf t = [] { gp_Trsf x; x.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0); return x; }();
    return t;
}
static TopoDS_Shape toMitreFrame(const TopoDS_Shape& metres) {
    TopoDS_Shape mm = BRepBuilderAPI_Transform(metres, mitreUpTrsf(), Standard_True).Shape();
    if (mm.IsNull()) return mm;
    TopoDS_Shape tight = BRepBuilderAPI_Copy(mm).Shape();
    ShapeFix_ShapeTolerance().SetTolerance(tight, Precision::Confusion());
    if (!tight.IsNull() && BRepCheck_Analyzer(tight).IsValid()) return tight;
    return mm;   // the shape really does need its slop: leave it, the probe below still judges it
}
// Back to metres. The inverse scale divides the tolerances too, which would leave the result
// claiming 1e-10 -- below the kernel's own confusion, which OCC reports as an invalid tolerance.
// Raise the too-small ones back to confusion and leave anything larger alone.
static TopoDS_Shape fromMitreFrame(const TopoDS_Shape& mm) {
    TopoDS_Shape metres = BRepBuilderAPI_Transform(mm, mitreDownTrsf(), Standard_True).Shape();
    if (!metres.IsNull()) ShapeFix_ShapeTolerance().LimitTolerance(metres, Precision::Confusion());
    return metres;
}

// DID THE MITRE ACTUALLY HAPPEN? ABT #860 (2026-08-23). Every previous verdict on a mitre cut came
// from the SAME boolean kernel that performs it -- 'removed' from the Cut, 'shouldRemove' from a
// Common with the knife -- so when the kernel goes blind (above) it reports success AND reports
// nothing left to remove, and the untrimmed overhang ships. The suggested cure of probing with a
// large half-space is blind in exactly the same way: measured on the same corner, Common against a
// 4691 mm3 box returned 0 mm3 while point classification found 0.058 mm3 of the solid inside it.
// So the check must not be a boolean at all. It asks the trimmed solid a question with a
// yes/no answer: is a point on the piece's OWN CENTRELINE, past the bisector plane, still inside
// it? The overhang runs from the piece's un-grown end `stubEnd` along `stubDir` for `grow`; a
// correct mitre leaves none of it. Only points standing at least kProbeClear past the plane are
// asked, so the answer never depends on the tolerance of the mitre face itself; an overhang
// thinner than that is thinner than the modelling tolerance and there is nothing to detect.
// Returns the largest stand-off (in metres) at which copper was found beyond the plane, or -1
// when the piece is clean.
static double overhangBeyondPlane(const TopoDS_Shape& trimmed, const gp_Pnt& stubEnd,
                                  const gp_Dir& stubDir, const gp_Pnt& P, const gp_XYZ& w,
                                  double grow) {
    if (trimmed.IsNull() || grow <= 0.0) return -1.0;
    const double rate = stubDir.XYZ().Dot(w);          // how fast the overhang leaves the plane
    if (rate <= 1e-9) return -1.0;                     // it does not head into the discard side
    const double s0 = (stubEnd.XYZ() - P.XYZ()).Dot(w);
    const double kProbeClear = 1e-6;                   // 1 um: 10x the shapes' own tolerance
    const double tCross = (kProbeClear - s0) / rate;   // first t standing clear of the plane
    if (tCross >= grow) return -1.0;                   // the grown end never reaches past it
    const double tFrom = std::max(0.0, tCross);
    double worst = -1.0;
    for (double f : {0.30, 0.55, 0.80, 0.97}) {
        const double t = tFrom + f * (grow - tFrom);
        const gp_Pnt probe(stubEnd.XYZ() + stubDir.XYZ() * t);
        const double stand = (probe.XYZ() - P.XYZ()).Dot(w);
        if (stand < kProbeClear) continue;
        // Per SOLID: a trim that fragmented the piece is a compound, and the classifier
        // silently misclassifies one of those.
        for (TopExp_Explorer ex(trimmed, TopAbs_SOLID); ex.More(); ex.Next()) {
            BRepClass3d_SolidClassifier c(TopoDS::Solid(ex.Current()), probe, Precision::Confusion());
            if (c.State() == TopAbs_IN) { worst = std::max(worst, stand); break; }
        }
    }
    return worst;
}

// Trim the grown overhang of `solid` that pokes past the mitre plane (through P, material side =
// keepDir): SUBTRACT a small knife box sitting on the discard side, localized to the junction.
// Neither of the two 'obvious' cuts is correct here:
//  - a global half-space keep is WRONG for near-closed wraps -- the wrap legitimately curves back
//    across the junction plane far from the joint, so the far side would be sliced off;
//  - the previous 'tight beam' knife (8r-wide box, depth spanning the solid, kept via Common) was
//    NOT a half-space: it discarded everything outside the beam laterally. Measured on PQ33: it
//    kept 17% of the first wrap and split it into 3 solids -> disconnected slivers.
// Locality bounds: the grown overhang lies within `grow`+r of P along -keepDir and within
// ~2.5r+`grow` of P laterally (bisector tilt <= ~50 deg), so a (3r+grow)-lateral x (grow+2r)-deep
// knife covers it with margin. A volume guard rejects any trim that removed more than a stub's
// worth of material (i.e. the knife reached a neighbouring pass of the same solid).
// ABT #685 (Alf, 2026-08-18): a failed mitre THROWS. "Don't allow stuff to return stuff when
// things fail, make them throw!" -- an untrimmed mitre is not a degraded result, it is copper
// standing proud of the corner, and in FEM that is fatal.
//
// MVB_MITRE_KEEP=1 is the DEBUGGING escape: the same complaint goes to stderr, in full, and the
// untrimmed solid is handed back so the geometry reaches disk and the fault can be LOOKED AT.
// Never a quiet fallback -- it shouts either way, and the loud path is opt-in.
static TopoDS_Shape refuseTrim(const std::string& why, const TopoDS_Shape& untrimmed) {
    if (std::getenv("MVB_MITRE_KEEP")) {
        std::cerr << "\n*** MITRE REFUSED (MVB_MITRE_KEEP set -- geometry is NOT valid) ***\n"
                  << why << "\n*** the piece below keeps its growth and overshoots the joint ***\n"
                  << std::endl;
        return untrimmed;
    }
    throw std::runtime_error(why + " [set MVB_MITRE_KEEP=1 to export it anyway and look at it]");
}

// ABT #685 (Alf, 2026-08-18): " A rejected trim used to return the solid UNCHANGED and silently, so the piece kept the
// growth it was given for the mitre and visibly ran past the corner -- "the turn continues a bit
// further than the corner", on Primary parallel 0's turn 30 against its dragback. An untrimmed
// mitre is not a degraded result, it is copper in the wrong place, and in FEM that is fatal. It
// throws, naming the primitive and the junction, so the geometry that caused it gets fixed.
// A SPIRAL that closes a FULL REVOLUTION is the one piece whose own copper comes back beside
// its corner: its other end sits one pitch away along the column axis, and on a coil wound with
// touching turns that is exactly one wire radius from this end's surface. No knife can sever such
// a corner without biting the piece itself (measured: PQ33 pitch 0.959133 mm against 2r 0.959 mm,
// 0.13 um to spare; the single-switch forward is identical at 0.06 um). ABT #685 (Alf,
// 2026-08-18): the adjacency is an artefact of it being ONE SOLID -- a knife only ever subtracts
// from the solid it is applied to, so a pass belonging to a different primitive is untouchable by
// construction. Sweeping the revolution as two halves puts each half's other end 180 degrees
// away, where no knife reaches, and makes the pass beside the corner somebody else's solid.
// The halves overlap by half a wire radius of arc at a tangent seam and are fused back into ONE
// solid after both corners are trimmed, so a wrap still ships as a single named part.
static TopoDS_Shape localMitreTrim(const TopoDS_Shape& solid, const gp_Pnt& P, const gp_Dir& keepDir,
                                   const gp_Dir& stubDir, const gp_Pnt& stubEnd, double r,
                                   double grow, const std::string& what, double neighbourStep,
                                   const std::string& diagnosis) {
    // HOW FAR MAY EACH KNIFE AXIS TRAVEL? ABT #685 (Alf, 2026-08-18). The same wire's next pass
    // sits `neighbourStep` away ALONG THE COLUMN AXIS (one pitch per revolution -- and a
    // full-revolution wrap has its OWN other end exactly there), so its copper begins at
    // `neighbourStep - r`. A knife face may go that far and no further.
    //
    // The bound is PER AXIS, not one scalar. A single `reach` applied to every axis is wrong in
    // both directions: it clamped the DEPTH, which here runs almost perpendicular to the column
    // axis and never approaches the neighbour, while leaving the one axis that does point along
    // it (v) free to run to 1.6r. On the flyback's Primary that put the knife 0.511 mm into a gap
    // whose neighbouring copper starts at 0.334 mm: it gouged the turn one pitch below (Alf: "a
    // dent below the mitre, as if the mitre had cut into the turn below") and sheared off a
    // fragment that shipped as its own solid ("[2/2] ... shouldn't exist").
    const gp_XYZ kColumnAxis(0, 1, 0);
    const double neighbourGap = neighbourStep > 0.0 ? std::max(0.0, neighbourStep - r) : 0.0;
    auto axisLimit = [&](const gp_XYZ& e) {
        if (neighbourGap <= 0.0) return std::numeric_limits<double>::max();
        const double c = std::abs(e.Dot(kColumnAxis));
        return c < 1e-9 ? std::numeric_limits<double>::max() : neighbourGap / c;
    };
    // Oriented knife: the overhang stub extends along `stubDir` beyond the plane, so the box
    // follows its lateral DRIFT direction and stays tube-radius-thin on the perpendicular
    // (axial/pitch) axis -- a symmetric (3r+grow)-wide box nicked the neighbouring wrap pass
    // one pitch (~3r) away and detached a 0.37 mm fragment (measured, PQ33 junction 9->10).
    const gp_XYZ w = keepDir.XYZ() * -1.0;                              // into the discard side
    const double cosT = std::max(0.2, std::abs(stubDir.XYZ().Dot(w)));  // clamp: tilt <= ~78 deg
    const double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    // Knife extent. It must cover the overhang and NOTHING else: a CURVED piece bends back
    // into the far part of the box and gets a rectangular notch gouged out of its keep side
    // (visible as square voids beside every dragback). Depth is therefore just the overhang plus
    // the tube radius -- what the tilted section actually needs -- and the perpendicular half is
    // barely wider than the tube.
    const double ell = r / cosT + 0.25 * r;   // mitre ellipse semi-extent + margin
    const double drift = (grow + r) * sinT;   // stub tip's lateral reach
    // LOCALITY BY SIZE, NOT BY INTERSECTION. ABT #685 (Alf, 2026-08-18): the knife used to be a
    // box intersected with a sphere centred on the joint, to stop its straight sides biting a
    // curved piece further along. The sphere was a disaster as a *cutting tool*: its face meets
    // the tube almost tangentially, and OCC's boolean returns a shape it then reports as INVALID
    // whose volume is garbage -- on 10_emi the cut of a 0.421 mm3 wrap measured 0.006 mm3, so the
    // guard read a 0.415 mm3 "removal" from a knife that was only 0.0165 mm3, refused the trim,
    // and left the piece standing proud. That refusal was the reported "mitre overshoot".
    // The box alone is local as long as its EXTENTS are: no face of it lies further than `reach`
    // from the joint, which is the same bound the sphere enforced, with flat faces the boolean can
    // actually handle.
    const double depth = std::min(grow + r, axisLimit(w));
    // The lateral faces must CLEAR the tube, never grave it. ABT #685 (Alf, 2026-08-18): at a
    // near-tangent joint the cut normal runs along the wire, so u and v are both lateral and the
    // knife becomes a sleeve around the tube -- and a half-extent of 1.05r puts its side face
    // 5.5 um outside a 0.11 mm cylinder. A plane grazing a cylinder is the classic boolean
    // degeneracy: OCC returned an INVALID cut whose volume read 0.006 mm3 for a 0.42 mm3 wrap.
    // A mitre cuts the WHOLE cross-section, so the box must contain the tube laterally with
    // margin; `reach` still caps it away from the neighbouring pass.
    // The knife must COVER THE WHOLE STUB in the lateral directions. ABT #685 (Alf, 2026-08-18):
    // clamping these by `reach` let the box's far face stop INSIDE the tube, so the cut severed
    // only part of the section and left a sliver -- OCC then returned an invalid shape (measured
    // on the push-pull: the stub needed 3.114 mm of drift, the clamp gave 2.529 mm, and a 418 mm3
    // wrap came back as 7.46 mm3). `reach` guards the direction the knife DIGS (depth, below),
    // which is where a neighbouring pass of the same wire sits; it has no business trimming the
    // knife across its own cut.
    const double lateral = 1.6 * r;
    gp_XYZ u = stubDir.XYZ() - w * stubDir.XYZ().Dot(w);
    if (u.Modulus() > 1e-9) u.Normalize();
    else u = gp_Ax2(P, gp_Dir(w)).XDirection().XYZ();   // stub ~normal to plane: no drift
    const gp_XYZ v = w.Crossed(u);
    const double uPos = std::min(std::max(ell + drift, lateral), axisLimit(u));
    const double uNeg = std::min(std::max(ell, lateral), axisLimit(u));
    const double vHalf = std::min(lateral, axisLimit(v));
    // A mitre SEVERS the section: the knife must span the tube on every lateral axis. If the
    // neighbouring pass is so close that it cannot, there is no knife that both cuts this corner
    // and spares the turn beside it -- say so rather than emit a partial cut (which leaves a
    // sliver and an invalid solid) or a gouged neighbour.
    // The knife needs real clearance outside the tube, not a tangent touch: a face lying ON the
    // cylinder is the grazing degeneracy that makes OCC return debris. Measured: the flyback's
    // Primary has 4.6% of a radius to spare (pitch 0.6536, 2r 0.639) and cuts cleanly, while the
    // single-switch forward's Primary is wound TOUCHING (pitch 0.383553 against 2r 0.3835), which
    // put the knife face 0.06 um outside the tube and returned an invalid shape that removed
    // nothing at all (faces 3 -> 9, volume unchanged). 2% of the radius separates the two.
    const double kSeverClearance = 1.02;
    if (uNeg < kSeverClearance * r || vHalf < kSeverClearance * r) {
        return refuseTrim(
            "ConductorBuilder: at " + what + " the same wire's next pass is only " +
                std::to_string(neighbourStep * 1e3) + " mm away (its copper starts at " +
                std::to_string(neighbourGap * 1e3) + " mm) against a " + std::to_string(r * 1e3) +
                " mm wire radius, so no knife can sever this corner without cutting into it "
                "(the turns are wound touching, and this wrap is a full revolution, so the pass "
                "beside the corner is this very piece's own other end); " +
                diagnosis,
            solid);
    }
    gp_Pnt corner(P.XYZ() - u * uNeg - v * vHalf);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Ax2(corner, gp_Dir(w), gp_Dir(u)),
                                           uNeg + uPos, 2.0 * vHalf, depth).Shape();
    // ZERO REMOVAL = FAILED CUT (ABT #685, 2026-08-20). A mitred end always carries growth
    // (overS/overE >= mitreGrow > 0), so its knife must always remove material; OCC's Cut on
    // one steep spiral pipe (14_dab's steep-exit stub) returned the solid unchanged while
    // reporting Done+valid -- removed -9e-12 mm3 with the overhang tip verified inside the
    // box -- and the untrimmed wedge shipped overlapping the exit lead by 0.148 mm3. When the
    // local box removes nothing and the piece has no neighbouring pass to protect (the only
    // reason the knife is local), retry once with a half-space-scale knife; a cut that still
    // removes nothing is refused loudly.
    const bool canRetryHalfSpace = neighbourStep <= 0.0;
    TopoDS_Shape bigBox;
    if (canRetryHalfSpace) {
        const double big = 20.0 * (r + grow);
        gp_Pnt bigCorner(P.XYZ() - u * big - v * big);
        bigBox = BRepPrimAPI_MakeBox(gp_Ax2(bigCorner, gp_Dir(w), gp_Dir(u)),
                                     2.0 * big, 2.0 * big, big).Shape();
    }
    try {
        // AN INVALID INPUT CANNOT BE CUT: OCC's booleans return an invalid argument unchanged
        // while reporting Done (measured, 14_dab's steep-exit stub: a 2.5-degree steep spiral
        // pipe whose sweep comes out invalid -- the local knife AND a half-space both "removed"
        // -9e-12 mm3). Repair the piece first; if it will not repair, the refusal below names it.
        TopoDS_Shape solidIn = solid;
        if (!BRepCheck_Analyzer(solidIn).IsValid()) {
            if (std::getenv("MVB_MITRE_DIAG"))
                std::cerr << "[trim-input-invalid] " << what << " -- repairing before the cut"
                          << std::endl;
            try {
                ShapeFix_Shape fixIn(solidIn);
                fixIn.Perform();
                if (!fixIn.Shape().IsNull() && BRepCheck_Analyzer(fixIn.Shape()).IsValid())
                    solidIn = fixIn.Shape();
            } catch (const Standard_Failure&) {
            }
        }
        GProp_GProps gpBefore, gpAfter;
        BRepGProp::VolumeProperties(solidIn, gpBefore);
        // CUT IN THE MILLIMETRE FRAME. ABT #685 (Alf, 2026-08-18): rawGrownSolid already sweeps
        // spiral/blend pipes at x1000 because OCC mis-frames sub-millimetre profiles in metres --
        // the boolean has exactly the same problem one step later. Measured on the push-pull's
        // 90-degree layer-link corner, in METRES: Cut(wrap, knife) returned a 7.46 mm3 fragment of
        // a 418 mm3 wrap, and Common(wrap, knife) returned ZERO -- OCC could not even see that the
        // two solids touch, though the joint sits exactly on the wrap's centreline. Scaling is an
        // EXACT affine map both ways, so this is a numerical frame choice, not an approximation.
        const TopoDS_Shape solidMm = toMitreFrame(solidIn);
        const TopoDS_Shape boxMm = toMitreFrame(box);
        BRepAlgoAPI_Cut cut(solidMm, boxMm);
        if (cut.IsDone() && !cut.Shape().IsNull()) {
            TopoDS_Shape trimmed = fromMitreFrame(cut.Shape());
            // REPAIR BEFORE JUDGING. ABT #685 (Alf, 2026-08-18): cutting a swept-pipe wrap leaves a
            // shell OCC reports as invalid -- and the volume of an invalid shell is meaningless
            // (the divergence integral needs a closed one). Measured on the push-pull: a 418 mm3
            // wrap cut by a 31 mm3 knife came back reading 7.46 mm3, so the guard below saw a
            // 411 mm3 "removal" and refused a trim that had actually done the right thing (faces
            // 5 -> 6: the mitre face WAS added). ShapeFix is the same repair the assembler already
            // runs on every trimmed piece for the STEP round-trip; running it here, before the
            // volume is read, means the guard judges the repaired solid instead of the debris.
            const bool rawValid = BRepCheck_Analyzer(trimmed).IsValid();
            double rawVol = -1;
            { GProp_GProps g; BRepGProp::VolumeProperties(trimmed, g); rawVol = g.Mass(); }
            if (!BRepCheck_Analyzer(trimmed).IsValid()) {
                try {
                    ShapeFix_Shape fix(trimmed);
                    fix.Perform();
                    if (!fix.Shape().IsNull() && BRepCheck_Analyzer(fix.Shape()).IsValid())
                        trimmed = fix.Shape();
                } catch (const Standard_Failure&) {
                }
            }
            // A shape OCC reports as invalid has a MEANINGLESS volume, so the guard below cannot
            // be trusted on one -- that is exactly how a good trim got read as a 98.5% removal.
            if (!BRepCheck_Analyzer(trimmed).IsValid()) {
                if (std::getenv("MVB_MITRE_DIAG")) {
                    GProp_GProps gk, ga;
                    BRepGProp::VolumeProperties(box, gk);
                    BRepGProp::VolumeProperties(trimmed, ga);
                    int nf = 0, nsol = 0;
                    for (TopExp_Explorer e(trimmed, TopAbs_FACE); e.More(); e.Next()) ++nf;
                    for (TopExp_Explorer e(trimmed, TopAbs_SOLID); e.More(); e.Next()) ++nsol;
                    int nfIn = 0;
                    for (TopExp_Explorer e(solid, TopAbs_FACE); e.More(); e.Next()) ++nfIn;
                    std::cerr << "[mitre-invalid] " << what << "\n   r=" << r * 1e3
                              << " grow=" << grow * 1e3 << " nstep=" << neighbourStep * 1e3
                              << " depth=" << depth * 1e3 << " uNeg=" << uNeg * 1e3
                              << " uPos=" << uPos * 1e3 << " vHalf=" << vHalf * 1e3
                              << " mm\n   knife=" << gk.Mass() * 1e9 << " before="
                              << gpBefore.Mass() * 1e9 << " after=" << ga.Mass() * 1e9
                              << " mm3, faces " << nfIn << "->" << nf << ", solids " << nsol
                              << "\n   cosT=" << cosT << " sinT=" << sinT << " drift="
                              << drift * 1e3 << " ell=" << ell * 1e3 << " mm" << std::endl;
                }
                return refuseTrim("ConductorBuilder: the mitre cut at " + what +
                                      " produced a shape OCC reports as INVALID, so neither the"
                                      " geometry nor its volume can be trusted; " + diagnosis,
                                  solid);
            }
            BRepGProp::VolumeProperties(trimmed, gpAfter);
            double removed = gpBefore.Mass() - gpAfter.Mass();
            if (removed < 1e-15 && canRetryHalfSpace && !bigBox.IsNull()) {
                // The local knife cut nothing off a grown end: retry with the half-space knife.
                const TopoDS_Shape bigMm = toMitreFrame(bigBox);
                BRepAlgoAPI_Cut cut2(solidMm, bigMm);
                if (cut2.IsDone() && !cut2.Shape().IsNull()) {
                    TopoDS_Shape trimmed2 = fromMitreFrame(cut2.Shape());
                    if (!BRepCheck_Analyzer(trimmed2).IsValid()) {
                        try {
                            ShapeFix_Shape fix2(trimmed2);
                            fix2.Perform();
                            if (!fix2.Shape().IsNull() &&
                                BRepCheck_Analyzer(fix2.Shape()).IsValid())
                                trimmed2 = fix2.Shape();
                        } catch (const Standard_Failure&) {
                        }
                    }
                    if (BRepCheck_Analyzer(trimmed2).IsValid()) {
                        GProp_GProps g2;
                        BRepGProp::VolumeProperties(trimmed2, g2);
                        const double removed2 = gpBefore.Mass() - g2.Mass();
                        if (removed2 > 1e-15) {
                            if (std::getenv("MVB_MITRE_DIAG"))
                                std::cerr << "[trim-retry] " << what
                                          << " half-space knife removed "
                                          << removed2 * 1e9 << " mm3 after the local box "
                                          << "removed nothing" << std::endl;
                            trimmed = trimmed2;
                            gpAfter = g2;
                            removed = removed2;
                        }
                    }
                }
            }
            // THE ONLY VERDICT THAT IS NOT THE BOOLEAN'S OWN. ABT #860 (2026-08-23): both of the
            // old verdicts came from the kernel that does the cutting -- `removed` from the Cut and
            // a Common against the knife -- so a kernel that goes blind passes both while the
            // overhang still stands (it did, on every --segments 0 terminal stub). Ask the trimmed
            // solid directly whether copper of its own is still standing past the bisector plane.
            const double leftBeyond =
                overhangBeyondPlane(trimmed, stubEnd, stubDir, P, w, grow);
            if (leftBeyond > 0.0) {
                std::ostringstream why2;
                why2 << "ConductorBuilder: the mitre cut at " << what << " left the overhang"
                     << " standing -- a point on the piece's own centreline " << leftBeyond * 1e3
                     << " mm past the bisector plane is still INSIDE the trimmed solid (the cut"
                     << " reported success and removed " << removed * 1e9 << " mm3 of "
                     << gpBefore.Mass() * 1e9 << " mm3; growth " << grow * 1e3 << " mm, wire r "
                     << r * 1e3 << " mm, half-space retry "
                     << (canRetryHalfSpace ? "available" : "unavailable") << "); " << diagnosis;
                return refuseTrim(why2.str(), solid);
            }
            if (removed < 1e-15 && std::getenv("MVB_MITRE_DIAG"))
                std::cerr << "[trim-zero-ok] " << what
                          << " grown end reaches but does not cross the plane" << std::endl;
            // stub bound: the overhang is at most a full-radius tube of length grow+2r, doubled
            // for the tilted-ellipse wedge. More than that = the knife ate distant material.
            const double stubMax = 2.0 * kPi * r * r * (grow + 2.0 * r);
            if (removed <= stubMax && gpAfter.Mass() > 0.0) {
                if (std::getenv("MVB_MITRE_DIAG")) {
                    std::cerr << "[trim] " << what << " removed=" << removed * 1e9
                              << " mm3 (grow=" << grow * 1e3 << " mm)" << std::endl;
                    if (removed < 1e-15) {
                        Bnd_Box bs;
                        BRepBndLib::Add(solid, bs);
                        double sx0, sy0, sz0, sx1, sy1, sz1;
                        bs.Get(sx0, sy0, sz0, sx1, sy1, sz1);
                        std::cerr << "[trim-miss] P=(" << P.X() * 1e3 << "," << P.Y() * 1e3
                                  << "," << P.Z() * 1e3 << ") w=(" << w.X() << "," << w.Y()
                                  << "," << w.Z() << ") stubDir=(" << stubDir.X() << ","
                                  << stubDir.Y() << "," << stubDir.Z() << ") depth="
                                  << depth * 1e3 << " uPos=" << uPos * 1e3 << " uNeg="
                                  << uNeg * 1e3 << " vHalf=" << vHalf * 1e3
                                  << " solidBB=[" << sx0 * 1e3 << "," << sy0 * 1e3 << ","
                                  << sz0 * 1e3 << "]..[" << sx1 * 1e3 << "," << sy1 * 1e3
                                  << "," << sz1 * 1e3 << "]" << std::endl;
                    }
                    int ns = 0;
                    for (TopExp_Explorer e(trimmed, TopAbs_SOLID); e.More(); e.Next()) ++ns;
                    if (ns != 1)
                        std::cerr << "[mitre-split] " << what << " -> " << ns << " solids"
                                  << ", removed " << removed * 1e9 << " mm3 of "
                                  << gpBefore.Mass() * 1e9 << ", r=" << r * 1e3
                                  << " grow=" << grow * 1e3 << " nstep=" << neighbourStep * 1e3
                                  << " depth=" << depth * 1e3 << " uNeg=" << uNeg * 1e3
                                  << " uPos=" << uPos * 1e3 << " vHalf=" << vHalf * 1e3
                                  << " mm; " << diagnosis << std::endl;
                }
                return trimmed;
            }
            if (std::getenv("MVB_MITRE_DIAG")) {
                std::cerr << "[stages] raw cut valid=" << rawValid << " vol=" << rawVol * 1e9
                          << " mm3 -> after fix vol=" << gpAfter.Mass() * 1e9
                          << " mm3 (input " << gpBefore.Mass() * 1e9 << " mm3)" << std::endl;
            }
            std::ostringstream why;
            why << "ConductorBuilder: the mitre knife at " << what << " REMOVED "
                << removed * 1e9 << " mm3 against a " << stubMax * 1e9
                << " mm3 bound for a corner stub"
                << " (before " << gpBefore.Mass() * 1e9 << " mm3, after "
                << gpAfter.Mass() * 1e9 << " mm3)"
                << "; knife at (" << P.X() * 1e3 << "," << P.Y() * 1e3 << "," << P.Z() * 1e3
                << ") mm, keep dir (" << keepDir.X() << "," << keepDir.Y() << "," << keepDir.Z()
                << "), wire r=" << r * 1e3 << " mm, growth=" << grow * 1e3 << " mm, neighbour step="
                << neighbourStep * 1e3 << " mm, " << diagnosis << "."
                << " The knife reached material away from the joint, so the trim is refused and"
                << " the piece would stand proud of the corner by its growth.";
            return refuseTrim(why.str(), solid);
        }
    } catch (const Standard_Failure& f) {
        return refuseTrim("ConductorBuilder: the mitre knife at " + what + " threw in OCC (" +
                              std::string(f.GetMessageString() ? f.GetMessageString() : "(null)") +
                              "); knife at (" + std::to_string(P.X() * 1e3) + "," +
                              std::to_string(P.Y() * 1e3) + "," + std::to_string(P.Z() * 1e3) +
                              ") mm, wire r=" + std::to_string(r * 1e3) + " mm, growth=" +
                              std::to_string(grow * 1e3) + " mm",
                          solid);
    }
    return refuseTrim("ConductorBuilder: the mitre cut at " + what + " did not complete (OCC "
                      "reported not-done or an empty result), so the piece would stand proud of "
                      "the corner by its " + std::to_string(grow * 1e3) + " mm growth",
                      solid);
}

TopoDS_Shape assembleWire(const std::vector<const Primitive*>& ptrs, double wireRadius,
                          int segments, CornerStyle corners,
                          std::vector<size_t>* primIndexPerSolid) {
    // WHEN IS A JUNCTION A CORNER? ABT #685 (Alf, 2026-08-18). Not "below 3 degrees", which was
    // another chosen number. A junction that is NOT mitred is BRIDGED: the earlier piece grows
    // flush past the joint until it fills the wedge the direction change opens on the outer side
    // of the bend, r*tan(theta), and the union swallows it. The only cost of bridging is that the
    // bridge stub pokes out of the neighbouring tube, by r*tan(theta)*sin(theta).
    // So mitre exactly when that poke-out is something this model can even REPRESENT -- i.e. when
    // it exceeds the sweep's own chordal sag, the error the faceted tube already carries
    // (samplingSag). Below that the mitre is beneath the model's resolution, and bridging is
    // strictly better: no boolean to fail, and no wedge gap either.
    // This matters because the boolean DOES fail. On 10_emi a 3.47 deg riser/wrap joint (poke-out
    // 0.4 um against a 2.2 um sag) was mitred, and OCC could not cut that spiral pipe surface at
    // ANY slice thickness -- 143 um or 3.3 um both came back INVALID, with the volume collapsing
    // from 0.411 mm3 to 0.0008 mm3 and a face LOST. The trim was refused, the piece stood proud,
    // and that was the reported "mitre overshoot".
    const double sag = samplingSag(wireRadius);
    // Growth that fills the wedge at a bridged (un-mitred) joint. Clamped well inside the range
    // where bridging is ever chosen, so tan() cannot run away.
    auto bridgeGrow = [&](double ang) { return wireRadius * std::tan(std::min(ang, 0.5)); };
    auto worthMitring = [&](double ang) { return bridgeGrow(ang) * std::sin(ang) > sag; };
    // |unit + unit| = 2 cos(angle/2): 0.2 admits joints up to ~168 degrees and rejects the folds
    // beyond, whose bisector carries no usable direction (ABT #685).
    const double kMinBisector = 0.2;
    // GROWTH PER MITRED END -- derived, not chosen. ABT #685 (Alf, 2026-08-18: "why 1.3? that
    // looks a magic number" -- it was). A bisector plane cuts a tube of radius r whose axis meets
    // the plane NORMAL at half the joint angle, so the ellipse it carves reaches exactly
    // r*tan(theta/2) past the joint along the axis. Grow that much and the plane is guaranteed to
    // find material everywhere it cuts, with nothing left over to slice away.
    // The old fixed 1.3*wireRadius is tan(52 deg): sized for a ~105 deg worst case and then
    // applied to every junction. It was wrong at BOTH ends -- on 10_emi's 3.47 deg riser/wrap
    // joint the mitre needs 3.3 um and it grew 143 um, a 43x stub that OCC then could not slice
    // off a spiral pipe surface (the cut came back INVALID, so the trim was refused and the piece
    // stood proud: the reported "mitre overshoot"); and past ~105 deg it grew LESS than the plane
    // cuts, leaving the mitre face notched.
    // The cap is not a new constant -- it is the SAME conditioning limit as kMinBisector:
    // |a+b| = 2cos(theta/2) > 0.2 admits theta < 168.5 deg, so the growth is bounded by ~10r.
    const double kMaxMitreAngle = 2.0 * std::acos(0.5 * kMinBisector);
    auto mitreGrow = [&](double ang) {
        return wireRadius * std::tan(0.5 * std::min(ang, kMaxMitreAngle));
    };
    const bool diag = std::getenv("MVB_MITRE_DIAG") != nullptr;
    BRep_Builder builder;
    TopoDS_Compound compound;
    builder.MakeCompound(compound);
    size_t n = ptrs.size();
    if (n == 0) return compound;
    std::vector<gp_Dir> fs(n), fe(n);
    for (size_t i = 0; i < n; ++i) {
        fs[i] = primFwdStart(*ptrs[i], wireRadius);
        fe[i] = primFwdEnd(*ptrs[i], wireRadius);
    }
    int nCut = 0, nRepaired = 0, nInvalid = 0;
    std::vector<TopoDS_Shape> built;
    // MITRED-JOINT ABUTMENT CHECK (ABT #685, Alf 2026-08-20: "there is no real bisection, just
    // overlapping solids ... learn how to detect those and add that to the checkings"). A
    // mitred corner slices BOTH sides on one shared bisector plane, so the two solids must
    // ABUT: their boolean COMMON is empty. Any residual volume is copper counted twice -- not
    // a valid FEM assembly and not a real winding. Checked for every CORNER joint as it is
    // built; bridged (near-tangent) joints are exempt, their wedge fill overlaps by design.
    // The floor is a (10 um)^3 cube in the mm frame: orders of magnitude above the boolean's
    // face-coincidence noise and orders below any real un-bisected corner.
    TopoDS_Shape prevBuilt;
    std::vector<std::string> overlappingJoints;
    // ABT #685: a silent check is worth nothing -- "no overlap reported" must not be able to mean
    // "the boolean never answered". Count the joints that got a real verdict against the ones that
    // did not, so a clean corpus can be told apart from a blind one (MVB_MITRE_DIAG prints both).
    size_t abutChecked = 0, abutNoVerdict = 0;
    // MEASURE ONLY WHAT THE KERNEL CAN VALIDLY PRODUCE, AND TAKE THE WORST OF IT. ABT #860
    // (2026-08-23): this gate used to read the volume of whatever BRepAlgoAPI_Common handed back,
    // in the millimetre frame, without asking whether the result was a valid solid -- the one thing
    // the rest of this file already knows never to do ("the volume of an invalid shell is
    // meaningless"). Measured on 13_current_sense at --segments 0: the Common of two VALID solids
    // came back INVALID reading 0.000973882 mm3 -- 87% of the 0.001118 mm3 dragback stub, i.e. the
    // stub reported as very nearly buried inside its neighbour -- and the design was refused for a
    // corner that does not overlap at all. The SAME two solids intersect in exactly 0, valid, when
    // the boolean runs in metres (and after a BRep round trip, and at any fuzzy value >= 1e-5 mm);
    // no point of the stub classifies inside the neighbour; and the faceted build of that very
    // design reports the joint clean. So: run the boolean in BOTH frames, discard any answer OCC
    // reports as invalid, and take the LARGEST of the trustworthy ones. That is conservative in the
    // direction a gate must be conservative -- a valid answer that shows overlap is never
    // discarded, and the frame that goes blind (also measured, see toMitreFrame) cannot argue a
    // corner clean on its own. If NEITHER frame gives a valid answer the joint is counted as
    // unchecked, loudly, and never as clean.
    auto checkMitredAbutment = [&](const TopoDS_Shape& a, const TopoDS_Shape& b,
                                   const std::string& joint) {
        if (a.IsNull() || b.IsNull()) { ++abutNoVerdict; return; }
        double worst = -1.0;                 // mm^3, over the frames that answered validly
        Bnd_Box worstBox;
        // scale: how many mm^3 one cubic unit of that frame's result is worth
        const std::array<double, 2> frameScale{1.0, 1e9};
        for (int frame = 0; frame < 2; ++frame) {
            try {
                const TopoDS_Shape x = frame == 0 ? toMitreFrame(a) : a;
                const TopoDS_Shape y = frame == 0 ? toMitreFrame(b) : b;
                BRepAlgoAPI_Common common(x, y);
                if (!common.IsDone() || common.Shape().IsNull()) continue;
                if (!BRepCheck_Analyzer(common.Shape()).IsValid()) continue;   // debris, not a volume
                GProp_GProps gcp;
                BRepGProp::VolumeProperties(common.Shape(), gcp);
                const double mm3 = gcp.Mass() * frameScale[frame];
                if (mm3 > worst) {
                    worst = mm3;
                    Bnd_Box cb;
                    BRepBndLib::Add(common.Shape(), cb);
                    if (frame == 1 && !cb.IsVoid()) {   // report every box in mm
                        double bx0, by0, bz0, bx1, by1, bz1;
                        cb.Get(bx0, by0, bz0, bx1, by1, bz1);
                        Bnd_Box scaled;
                        scaled.Update(bx0 * 1e3, by0 * 1e3, bz0 * 1e3,
                                      bx1 * 1e3, by1 * 1e3, bz1 * 1e3);
                        cb = scaled;
                    }
                    worstBox = cb;
                }
            } catch (const Standard_Failure&) {
                // The boolean itself failing is a different defect; the invalid-cut paths
                // already throw where the trim happens. Let the other frame answer.
            }
        }
        if (worst < 0.0) { ++abutNoVerdict; return; }   // no frame produced a valid answer
        ++abutChecked;
        if (worst > 1e-6) {
            std::ostringstream where;
            where << joint << " overlaps by " << worst << " mm^3";
            if (!worstBox.IsVoid()) {
                double cx0, cy0, cz0, cx1, cy1, cz1;
                worstBox.Get(cx0, cy0, cz0, cx1, cy1, cz1);
                where << ", common bbox mm [" << cx0 << "," << cy0 << "," << cz0 << "]..["
                      << cx1 << "," << cy1 << "," << cz1 << "]";
            }
            overlappingJoints.push_back(where.str());
        }
    };
    for (size_t i = 0; i < n; ++i) {
        // Bend at each end: a tangent joint (wrap SEG<->ARC) needs no growth and no boolean; only a
        // real corner (the entrance/exit leads) is grown and sliced on its angle-bisector plane.
        // Junctions may carry an ENDPOINT MISMATCH (layouts hand out lead ends offset from the wrap
        // start by up to a wire radius -- measured 0.36 mm on PQ33): the mitre plane goes through
        // the MIDPOINT of the two endpoints (same plane for both sides) and the growth is widened
        // by the mismatch so the overhang always reaches past the shared plane.
        const double dpS = (i > 0)
            ? primEndpoints(*ptrs[i - 1]).second.Distance(primEndpoints(*ptrs[i]).first) : 0.0;
        const double dpE = (i + 1 < n)
            ? primEndpoints(*ptrs[i]).second.Distance(primEndpoints(*ptrs[i + 1]).first) : 0.0;
        bool bentS = false;
        gp_Dir nS = fs[i];
        if (diag && i > 0) {
            static const char* kn[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            std::cerr << "[slope]   " << ptrs[i - 1]->label << "  dy/ds=" << fe[i - 1].Y()
                      << "   |   " << ptrs[i]->label << "  dy/ds=" << fs[i].Y() << "\n";
            std::cerr << "[mitre]   junction " << i - 1 << "->" << i << " " << kn[ptrs[i - 1]->kind]
                      << "->" << kn[ptrs[i]->kind] << " angle=" << fe[i - 1].Angle(fs[i]) * 180.0 / kPi
                      << " deg" << (worthMitring(fe[i - 1].Angle(fs[i])) ? " CORNER" : " bridged")
                      << " ['" << ptrs[i - 1]->label << "' -> '" << ptrs[i]->label << "']\n";
        }
        // ABT #685 (Alf, 2026-08-14): a lead/turn corner is mitred exactly like the elbow between
        // the two lead runs themselves (the vertical bridge into the radial run out) — one bisector
        // plane, both sides grown to it and sliced, so on a 90 degree corner the cut is at 45 and
        // the two pieces share the face. Tried and rejected on the way here: slicing on the LEAD's
        // face, and leaving both flush with a wire-radius overshoot — the first is a long slanted
        // ellipse, the second leaves the corner's outer side notched.
        if (i > 0 && worthMitring(fe[i - 1].Angle(fs[i]))) {
            gp_Vec s(fe[i - 1].XYZ());
            s += gp_Vec(fs[i].XYZ());
            // ABT #685 (Alf, 2026-08-18): the bisector must be WELL CONDITIONED. Two unit
            // directions meeting at nearly 180 degrees sum to nearly nothing, and normalising
            // that gives an essentially arbitrary plane -- the knife then keeps the wrong side
            // and eats the piece. Measured on 10_emi: a riser/wrap joint whose keep direction
            // came out (0, 0.03, -0.9995), pointing along the wrap rather than across it, and the
            // cut removed 0.415 mm3 of a 0.421 mm3 half-revolution. |sum| = 2 cos(angle/2), so
            // requiring 0.2 stops mitring anything past ~168 degrees; such a joint is a fold, not
            // a corner, and is left tangent -- ungrown and uncut, which cannot overshoot.
            if (s.Magnitude() > kMinBisector) { nS = gp_Dir(s); bentS = true; }
        }
        bool bentE = false;
        gp_Dir nE = fe[i];
        if (i + 1 < n && worthMitring(fe[i].Angle(fs[i + 1]))) {
            gp_Vec s(fe[i].XYZ());
            s += gp_Vec(fs[i + 1].XYZ());
            if (s.Magnitude() > kMinBisector) { nE = gp_Dir(s); bentE = true; }   // see above
        }
        // Growth per end: corners grow by over+mismatch (then trimmed on the mitre plane);
        // a TANGENT junction with an endpoint gap is bridged by growing the EARLIER prim's end
        // flush forward (no cut) -- one side only, so the bridge is never doubled.
        const double angS = (i > 0) ? fe[i - 1].Angle(fs[i]) : 0.0;
        const double angE = (i + 1 < n) ? fe[i].Angle(fs[i + 1]) : 0.0;
        const double overS = bentS ? mitreGrow(angS) + dpS : 0.0;
        // A bridged end grows by the WEDGE the bend opens (r*tan(theta)) as well as any endpoint
        // mismatch -- growing only by dpE left a wedge gap of up to r*tan(theta) on the outer side
        // of every near-tangent joint, open copper-to-copper. One side only, so it is never doubled.
        const double overE = bentE ? mitreGrow(angE) + dpE
                                   : dpE + (angE > 1e-12 ? bridgeGrow(angE) : 0.0);
        // SPLIT A CLOSED REVOLUTION. ABT #685 (Alf, 2026-08-18): a SPIRAL that closes a full
        // turn has its own other end one pitch away ALONG THE COLUMN AXIS -- and a coil wound
        // with touching turns puts that end's copper exactly one wire radius from this one's
        // surface (measured: PQ33 pitch 0.959133 against 2r 0.959000, 0.13 um to spare; the
        // single-switch forward is flush to 0.06 um). A mitre knife must span the whole section
        // to sever it, so on those coils NO box can cut the corner without biting the piece
        // itself: too wide and it gouges the turn below and shears off a fragment (the flyback's
        // "dent below the mitre" and its stray [2/2] solid), exactly wide enough and its face
        // lies ON the cylinder, which is the grazing degeneracy that makes the boolean return
        // debris.
        // Sweeping the revolution as two HALVES dissolves the conflict: each half's other end is
        // 180 degrees away, far outside any knife, and the pass one pitch along the axis now
        // belongs to a DIFFERENT solid -- which a knife may overlap freely, since it only ever
        // subtracts from its own. The halves overlap at a tangent seam and are fused back into
        // ONE solid after both corners are trimmed, so the part list is unchanged.
        auto revolutionHalves = [&](const Primitive& pr) {
            std::vector<Primitive> out;
            const Spiral& sp = pr.spiral;
            const bool closed = pr.kind == Primitive::SPIRAL && !sp.blend &&
                                std::abs(sp.r1 - sp.r0) < 1e-12 &&
                                std::abs(sp.az1 - sp.az0) > 0.97 * kTwoPi;
            if (!closed) {
                out.push_back(pr);
                return out;
            }
            const double azM = 0.5 * (sp.az0 + sp.az1);
            const double yM = 0.5 * (sp.y0 + sp.y1);
            // The halves ABUT EXACTLY at the mid-azimuth: same centreline point, same tangent,
            // meeting on one flush perpendicular disc -- the assembler's own junction doctrine.
            // The previous seam OVERLAPPED the halves along the same analytic curve, which made
            // their tube surfaces exactly coincident over the window: the one input class OCC's
            // fuse cannot handle. On the pushpull's steep wrap it returned "done", one valid
            // solid -- half A alone, half B silently gone (A=147.242, B=147.294, fused=147.240
            // mm3). Flush abutment plus glue mode below is the coincidence OCC does support.
            Primitive a = pr, b = pr;
            a.spiral.az1 = azM;
            a.spiral.y1 = yM;
            b.spiral.az0 = azM;
            b.spiral.y0 = yM;
            out.push_back(a);
            out.push_back(b);
            return out;
        };
        std::vector<Primitive> pieces = revolutionHalves(*ptrs[i]);
        auto sweepPieces = [&](double a0, double b0) {
            std::vector<TopoDS_Shape> got;
            for (size_t q = 0; q < pieces.size(); ++q) {
                TopoDS_Shape sp = rawGrownSolid(pieces[q], wireRadius,
                                                q == 0 ? a0 : 0.0,
                                                q + 1 == pieces.size() ? b0 : 0.0, segments);
                if (sp.IsNull()) return std::vector<TopoDS_Shape>{};
                got.push_back(sp);
            }
            return got;
        };
        std::vector<TopoDS_Shape> parts = sweepPieces(overS, overE);
        if (parts.empty()) {  // ARC clamp (near-full revolve): fall back to a flush, uncut tube
            parts = sweepPieces(0.0, 0.0);
            bentS = bentE = false;
        }
        TopoDS_Shape solid = parts.empty() ? TopoDS_Shape() : parts.front();
        if (solid.IsNull()) {
            // Never skip silently: a missing primitive is missing copper, and its neighbours
            // may still touch each other so even the connectivity contract can stay satisfied
            // (measured: faceted blend sweeps returned null and every rect-column turn quietly
            // detached from the next).
            static const char* kindName[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            const auto ends = primEndpoints(*ptrs[i]);
            throw std::runtime_error(
                "ConductorBuilder: the conformal assembler could not sweep primitive " +
                std::to_string(i) + " (" + kindName[ptrs[i]->kind] + ") '" + ptrs[i]->label +
                "' (" + std::to_string(ends.first.X()) + "," + std::to_string(ends.first.Y()) +
                "," + std::to_string(ends.first.Z()) + ")->(" + std::to_string(ends.second.X()) +
                "," + std::to_string(ends.second.Y()) + "," + std::to_string(ends.second.Z()) +
                "). Refusing to emit a conductor with missing copper.");
        }
        // How close the SAME wire comes to this joint on a neighbouring pass: for a spiral that
        // is its axial advance per revolution. The knife may not reach that far, or it cuts the
        // wire's other pass instead of the corner. Zero = unbounded (no neighbouring pass).
        auto junctionDiag = [](const gp_Dir& a, const gp_Dir& b, const Primitive& pr) {
            static const char* kn[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
            std::ostringstream d;
            d << "joint angle " << a.Angle(b) * 180.0 / kPi << " deg (in ("
              << a.X() << "," << a.Y() << "," << a.Z() << ") out (" << b.X() << "," << b.Y()
              << "," << b.Z() << ")), cut piece is " << kn[pr.kind];
            if (pr.kind == Primitive::SPIRAL) {
                d << " az span " << (pr.spiral.az1 - pr.spiral.az0) * 180.0 / kPi << " deg r "
                  << pr.spiral.r0 * 1e3 << " y " << pr.spiral.y0 * 1e3 << "->"
                  << pr.spiral.y1 * 1e3 << " mm";
            }
            return d.str();
        };
        auto reachLimitOf = [&](const Primitive& pr) {
            if (pr.kind != Primitive::SPIRAL) return 0.0;
            const double sweep = pr.spiral.az1 - pr.spiral.az0;
            if (std::abs(sweep) < 1e-9) return 0.0;
            // ONLY A CLOSED REVOLUTION HAS COPPER OF ITS OWN BESIDE THE CORNER. ABT #685 (Alf,
            // 2026-08-18): the knife subtracts from THIS solid and nothing else, so a neighbouring
            // pass that belongs to a different primitive is out of its reach by construction --
            // bounding against it only shrank the knife below the tube it had to sever. Applying
            // the bound to every spiral refused ordinary windings: PQ33's half-revolution wrap
            // 8->9 was rejected for a pass 0.959133 mm away that is not part of it at all, and
            // that design had been building for weeks. A full turn IS its own neighbour -- and is
            // split into halves above precisely so that neither half is.
            if (std::abs(sweep) < 0.97 * kTwoPi) return 0.0;
            const double pitch = std::abs((pr.spiral.y1 - pr.spiral.y0) * kTwoPi / sweep);
            if (pitch < 1e-12) return 0.0;
            return pitch;   // the neighbouring pass sits exactly one pitch along the axis
        };
        const double reachHere = reachLimitOf(*ptrs[i]);
        // Once split, a half's own other end is half a revolution away, so nothing of this piece
        // sits near either knife: the neighbour bound only applies to a still-closed revolution.
        const double neighbourHere = pieces.size() > 1 ? 0.0 : reachHere;
        if (bentS) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i - 1]).second.XYZ() +
                                  primEndpoints(*ptrs[i]).first.XYZ()));
            parts.front() =
                localMitreTrim(parts.front(), J, nS, gp_Dir(fs[i].XYZ() * -1.0),
                               primEndpoints(*ptrs[i]).first, wireRadius, overS,
                               "'" + ptrs[i - 1]->label + "' -> '" + ptrs[i]->label + "'",
                               neighbourHere, junctionDiag(fe[i - 1], fs[i], *ptrs[i]));
            ++nCut;
        }
        if (bentE) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i]).second.XYZ() +
                                  primEndpoints(*ptrs[i + 1]).first.XYZ()));
            parts.back() =
                localMitreTrim(parts.back(), J, gp_Dir(nE.XYZ() * -1.0), fe[i],
                               primEndpoints(*ptrs[i]).second, wireRadius, overE,
                               "'" + ptrs[i]->label + "' -> '" + ptrs[i + 1]->label + "'",
                               neighbourHere, junctionDiag(fe[i], fs[i + 1], *ptrs[i]));
            ++nCut;
        }
        solid = parts.empty() ? TopoDS_Shape() : parts.front();
        if (parts.size() > 1) {
            // Fuse the halves back: both corners are cut, so the closed revolution is safe again.
            gp_Trsf up, down;
            up.SetScale(gp_Pnt(0, 0, 0), 1000.0);
            down.SetScale(gp_Pnt(0, 0, 0), 1.0 / 1000.0);
            try {
                // Glue mode: the halves meet on one exactly-shared flush disc, and
                // BOPAlgo_GlueFull is OCC's declared contract for coincident faces -- the
                // general boolean treats them as an intersection problem and loses a half.
                BRepAlgoAPI_Fuse fu;
                TopTools_ListOfShape args, tools;
                args.Append(BRepBuilderAPI_Transform(parts[0], up, Standard_True).Shape());
                tools.Append(BRepBuilderAPI_Transform(parts[1], up, Standard_True).Shape());
                fu.SetArguments(args);
                fu.SetTools(tools);
                fu.SetGlue(BOPAlgo_GlueShift);
                fu.Build();
                if (!fu.IsDone() || fu.Shape().IsNull())
                    throw std::runtime_error("fuse of the two half-revolutions did not complete");
                solid = BRepBuilderAPI_Transform(fu.Shape(), down, Standard_True).Shape();
                // THE FUSE MUST ACCOUNT FOR EVERY CUBIC MICRON. ABT #685 (Alf, 2026-08-19,
                // "why is turn 2 -> turn 2_ending not connecting with the layer link?"): on the
                // pushpull's steep single-turn wrap, BRepAlgoAPI_Fuse reported IsDone, returned
                // one VALID solid -- and that solid was half A alone (A=147.242, B=147.294,
                // fused=147.240 mm3). Half B, the half that reaches the layer link, vanished
                // without a diagnostic; the conductor shipped visibly disconnected. IsDone and
                // IsValid cannot see this, only conservation can: the halves coincide exactly
                // over the seam (same analytic curve), so the union's volume must be
                // A + B - seamOverlap, where seamOverlap is the tube over the shared arc --
                // computable, not estimated. If the fuse comes up short, keep BOTH halves as
                // flush-overlapping solids (the weld-lens pattern the rect path already uses;
                // the consumer's fragment welds them) and say so.
                {
                    GProp_GProps gA, gB, gF;
                    BRepGProp::VolumeProperties(parts[0], gA);
                    BRepGProp::VolumeProperties(parts[1], gB);
                    BRepGProp::VolumeProperties(solid, gF);
                    // Flush abutment: no shared volume, the union must carry every cubic
                    // micron of both halves. The failure mode loses an entire half (~50%); one
                    // percent is pure measurement headroom, not a configuration in between.
                    const double expected = gA.Mass() + gB.Mass();
                    if (gF.Mass() < expected * 0.99) {
                        std::cerr << "[ConductorBuilder] fuse of '" << ptrs[i]->label
                                  << "' lost copper (A=" << gA.Mass() * 1e9
                                  << " B=" << gB.Mass() * 1e9 << " fused=" << gF.Mass() * 1e9
                                  << " mm3, expected >= " << expected * 0.99 * 1e9
                                  << "); keeping both halves as flush-abutting solids for the "
                                     "consumer's weld."
                                  << std::endl;
                        TopoDS_Compound both;
                        BRep_Builder bb2;
                        bb2.MakeCompound(both);
                        bb2.Add(both, parts[0]);
                        bb2.Add(both, parts[1]);
                        solid = both;
                    }
                    if (std::getenv("MVB_FUSE_DIAG")) {
                        int ns = 0;
                        for (TopExp_Explorer e(solid, TopAbs_SOLID); e.More(); e.Next()) ++ns;
                        std::cerr << "[fuse] '" << ptrs[i]->label << "' A=" << gA.Mass() * 1e9
                                  << " B=" << gB.Mass() * 1e9 << " fused=" << gF.Mass() * 1e9
                                  << " mm3, solids=" << ns << std::endl;
                        int si = 0;
                        for (TopExp_Explorer e(solid, TopAbs_SOLID); e.More(); e.Next(), ++si) {
                            Bnd_Box bx;
                            BRepBndLib::Add(e.Current(), bx);
                            double x0, y0z, z0, x1, y1z, z1;
                            bx.Get(x0, y0z, z0, x1, y1z, z1);
                            std::cerr << "[fuse]   solid " << si << " bbox x[" << x0 * 1e3 << ","
                                      << x1 * 1e3 << "] y[" << y0z * 1e3 << "," << y1z * 1e3
                                      << "] z[" << z0 * 1e3 << "," << z1 * 1e3 << "] mm"
                                      << std::endl;
                        }
                    }
                }
            } catch (const Standard_Failure& f) {
                throw std::runtime_error(
                    "ConductorBuilder: could not rejoin the two halves of the closed revolution '" +
                    ptrs[i]->label + "' (" +
                    std::string(f.GetMessageString() ? f.GetMessageString() : "(null)") + ")");
            }
        }
        // STEP round-trip robustness: the boolean's cut curve on a torus can carry an edge tolerance
        // that a STEP export degrades into an invalid solid (valid in memory, invalid on reload).
        // ShapeFix tightens edges/tolerances in place so the written solid survives the round-trip.
        if (!solid.IsNull() && (bentS || bentE)) {
            try {
                ShapeFix_Shape fix(solid);
                fix.Perform();
                if (!fix.Shape().IsNull()) {
                    if (std::getenv("MVB_FUSE_DIAG")) {
                        GProp_GProps g0, g1;
                        BRepGProp::VolumeProperties(solid, g0);
                        BRepGProp::VolumeProperties(fix.Shape(), g1);
                        if (std::abs(g1.Mass() - g0.Mass()) > 0.01 * std::abs(g0.Mass()))
                            std::cerr << "[shapefix] '" << ptrs[i]->label << "' volume "
                                      << g0.Mass() * 1e9 << " -> " << g1.Mass() * 1e9 << " mm3"
                                      << std::endl;
                    }
                    solid = fix.Shape();
                }
            } catch (const Standard_Failure&) {
            }
        }
        // A slice that self-intersected or emptied the tube (very short/sharp primitive) is caught
        // here and rebuilt as a plain flush-capped tube -- always valid, at worst a hair of overlap
        // at that single joint. Never emit an invalid solid.
        if (solid.IsNull() || !BRepCheck_Analyzer(solid).IsValid()) {
            TopoDS_Shape flush = rawGrownSolid(*ptrs[i], wireRadius, 0.0, 0.0, segments);
            if (!flush.IsNull() && BRepCheck_Analyzer(flush).IsValid()) { solid = flush; ++nRepaired; }
            else {
                // NEVER drop a primitive silently: the copper it carries simply disappears from
                // the conductor, and downstream checks can miss it (a short corner's neighbours
                // may still touch, so even the connectivity contract stays satisfied -- that is
                // exactly how a missing dragback corner reached a STEP). Fail loudly instead,
                // naming the primitive so the geometry that produced it can be fixed.
                static const char* kindName[] = {"SEG", "ARC3", "SPIRAL", "BLEND"};
                const auto ends = primEndpoints(*ptrs[i]);
                std::ostringstream detail;
                detail << "ConductorBuilder: the conformal assembler cannot build a valid solid for "
                       << "primitive " << i << " (" << kindName[ptrs[i]->kind] << ") '"
                       << ptrs[i]->label << "' (" << ends.first.X() << "," << ends.first.Y() << ","
                       << ends.first.Z() << ")->(" << ends.second.X() << "," << ends.second.Y()
                       << "," << ends.second.Z() << ")";
                if (ptrs[i]->kind == Primitive::ARC3) {
                    detail << " sweep=" << ptrs[i]->arc.sweep
                           << " centrelineRadius=" << ptrs[i]->arc.v0.Modulus()
                           << " wireRadius=" << wireRadius;
                    if (ptrs[i]->arc.v0.Modulus() <= wireRadius * (1.0 + 1e-6)) {
                        detail << " -- the centreline radius does not exceed the wire radius, so the "
                                  "swept tube is a HORN TORUS touching its own axis of revolution "
                                  "(no valid solid exists; give the corner a larger bend radius)";
                    }
                }
                detail << ". Refusing to emit a conductor with missing copper.";
                throw std::runtime_error(detail.str());
            }
        }
        // ABT #685: add SOLIDS, never a sub-compound. A mitre trim can fragment one primitive
        // into two solids, and adding that compound as one shape made the assembly's component
        // count differ from its solid count -- 21 solids, 20 components on the pushpull's
        // Primary 1. Every consumer that walks solids and then indexes components (the STEP
        // exporter's per-solid naming, above all) is off by one from that point on, which is how
        // 'turn 15 -> turn 16' ended up naming the wrong piece of copper. One solid, one
        // component, one name.
        if (bentS && i > 0) {
            checkMitredAbutment(prevBuilt, solid,
                                "'" + ptrs[i - 1]->label + "' -> '" + ptrs[i]->label + "'");
        }
        prevBuilt = solid;
        for (TopExp_Explorer solidExplorer(solid, TopAbs_SOLID); solidExplorer.More();
             solidExplorer.Next()) {
            builder.Add(compound, solidExplorer.Current());
            if (primIndexPerSolid != nullptr) {
                primIndexPerSolid->push_back(i);
            }
        }
        if (diag) {
            built.push_back(solid);
            GProp_GProps gpr;
            BRepGProp::VolumeProperties(solid, gpr);
            Bnd_Box bb;
            BRepBndLib::Add(solid, bb);
            double x0, y0, z0, x1, y1, z1;
            bb.Get(x0, y0, z0, x1, y1, z1);
            const double len = primEndpoints(*ptrs[i]).first.Distance(primEndpoints(*ptrs[i]).second);
            const double nominal = kPi * wireRadius * wireRadius * len;
            std::cerr << "[prim] " << i << " '" << ptrs[i]->label << "' vol=" << gpr.Mass()
                      << " nominal=" << nominal
                      << " ratio=" << (nominal > 0 ? gpr.Mass() / nominal : 0.0)
                      << " y=[" << y0 << "," << y1 << "]\n";
        }
    }
    if (abutNoVerdict > 0 || diag) {
        std::cerr << "[abut] mitred joints with a boolean verdict: " << abutChecked
                  << ", without one: " << abutNoVerdict
                  << (abutNoVerdict > 0 ? "  <-- those corners are UNCHECKED, not clean\n" : "\n");
    }
    if (!overlappingJoints.empty()) {
        std::string msg =
            "ConductorBuilder: " + std::to_string(overlappingJoints.size()) +
            " mitred corner(s) did not abut -- the two sides overlap instead of sharing the "
            "bisector face:";
        for (const auto& j : overlappingJoints) msg += "\n  " + j;
        // ABT #685: the same escape the certified enamel gate carries (MVB_ALLOW_ENAMEL) -- OFF by
        // default, so anything that ships is a design whose corners abut. Set it only to EXPORT a
        // known-faulty design for visual review; the fault is still printed in full.
        if (std::getenv("MVB_ALLOW_ABUTMENT")) {
            std::cerr << "[corners] " << msg
                      << "\n  (MVB_ALLOW_ABUTMENT set: reported, not refused)\n";
        } else {
            throw std::runtime_error(msg);
        }
    }
    // Per-solid analysis (union-find connectivity, volumes, interpenetration) is O(n^2) EXACT
    // BRepExtrema/classification -- ~15k queries on a 170-prim toroid, tens of minutes. Keep it
    // behind its own flag so plain MVB_MITRE_DIAG (junction angles, DROPPED prims, the counters)
    // stays usable during normal iteration.
    if (diag && std::getenv("MVB_MITRE_INTERPEN")) {
        // component structure over ALL pairs (mirrors the test's union-find), with volumes
        std::vector<int> uf(built.size());
        for (size_t i = 0; i < uf.size(); ++i) uf[i] = (int)i;
        std::function<int(int)> find = [&](int a) {
            while (uf[a] != a) a = uf[a] = uf[uf[a]];
            return a;
        };
        for (size_t i = 0; i < built.size(); ++i)
            for (size_t j = i + 1; j < built.size(); ++j) {
                try {
                    BRepExtrema_DistShapeShape d(built[i], built[j]);
                    if (d.IsDone() && d.Value() < 1e-6) uf[find((int)i)] = find((int)j);
                } catch (const Standard_Failure&) {
                }
            }
        for (size_t i = 0; i < built.size(); ++i) {
            GProp_GProps gp;
            BRepGProp::VolumeProperties(built[i], gp);
            int nsol = 0;
            std::string sub;
            for (TopExp_Explorer ex(built[i], TopAbs_SOLID); ex.More(); ex.Next()) {
                ++nsol;
                GProp_GProps gs;
                BRepGProp::VolumeProperties(ex.Current(), gs);
                char vb[32]; std::snprintf(vb, sizeof(vb), " %.4g", gs.Mass()); sub += vb;
            }
            std::cerr << "[mitre]   solid " << i << " comp=" << find((int)i) << " vol=" << gp.Mass()
                      << " nsolids=" << nsol << (nsol > 1 ? " subvols:" + sub : "") << "\n";
        }
        for (size_t i = 1; i < built.size(); ++i) {
            const double dp = primEndpoints(*ptrs[i - 1]).second.Distance(primEndpoints(*ptrs[i]).first);
            double ds = -1.0;
            try {
                BRepExtrema_DistShapeShape ext(built[i - 1], built[i]);
                if (ext.IsDone()) ds = ext.Value();
            } catch (const Standard_Failure&) {
            }
            if (dp > 1e-9 || ds > 1e-9)
                std::cerr << "[mitre]   junction " << i - 1 << "->" << i << " endpointGap=" << dp
                          << " solidGap=" << ds << "\n";
        }
        // Interpenetration sweep over ALL bbox-touching pairs (not just chain neighbours), the
        // same junction-grid classifier as the test battery, with LABELS -- post-prune solid
        // indices in test output cannot be mapped back to primitives. O(n^2) solid
        // classification dominates the run (tens of minutes on a 170-prim toroid), so it needs
        // its OWN flag: plain MVB_MITRE_DIAG stays fast enough for routine use.
        {
            std::vector<Bnd_Box> bxs(built.size());
            for (size_t i = 0; i < built.size(); ++i) BRepBndLib::Add(built[i], bxs[i]);
            for (size_t i = 0; i < built.size(); ++i)
                for (size_t j = i + 1; j < built.size(); ++j) {
                    if (bxs[i].Distance(bxs[j]) > 1e-9) continue;
                    double x0, y0, z0, x1, y1, z1, u0, v0, w0, u1, v1, w1;
                    bxs[i].Get(x0, y0, z0, x1, y1, z1);
                    bxs[j].Get(u0, v0, w0, u1, v1, w1);
                    const double ax = std::max(x0, u0), bx = std::min(x1, u1);
                    const double ay = std::max(y0, v0), by = std::min(y1, v1);
                    const double az = std::max(z0, w0), bz = std::min(z1, w1);
                    if (ax >= bx || ay >= by || az >= bz) continue;
                    // classify per SOLID: a built[] entry can be a multi-solid compound (a trim
                    // that fragmented), and BRepClass3d_SolidClassifier on a compound silently
                    // misclassifies -- exploding is what makes this sweep agree with the test
                    // battery's flat-solid enumeration.
                    auto insideAny = [](const TopoDS_Shape& sh, const gp_Pnt& pp) {
                        for (TopExp_Explorer ex(sh, TopAbs_SOLID); ex.More(); ex.Next()) {
                            BRepClass3d_SolidClassifier c(TopoDS::Solid(ex.Current()), pp, 1e-9);
                            if (c.State() == TopAbs_IN) return true;
                        }
                        return false;
                    };
                    int hits = 0;
                    gp_Pnt hit;
                    constexpr int N = 6;
                    for (int gx = 0; gx < N && hits == 0; ++gx)
                        for (int gy = 0; gy < N && hits == 0; ++gy)
                            for (int gz = 0; gz < N && hits == 0; ++gz) {
                                gp_Pnt pp(ax + (bx - ax) * (gx + 0.5) / N,
                                          ay + (by - ay) * (gy + 0.5) / N,
                                          az + (bz - az) * (gz + 0.5) / N);
                                if (insideAny(built[i], pp) && insideAny(built[j], pp)) {
                                    ++hits;
                                    hit = pp;
                                }
                            }
                    if (hits > 0)
                        std::cerr << "[mitre]   INTERPENETRATION prims " << i << "<->" << j
                                  << " ['" << ptrs[i]->label << "' vs '" << ptrs[j]->label
                                  << "'] at (" << hit.X() << "," << hit.Y() << "," << hit.Z()
                                  << ")\n";
                }
        }
    }
    if (diag) {
        std::cerr << "[mitre] prims=" << n << " boolean-cuts=" << nCut << " repaired=" << nRepaired
                  << " dropped-invalid=" << nInvalid << "\n";
    }
    return compound;
}

} // namespace mvb
