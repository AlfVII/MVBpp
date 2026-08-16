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
    int n = curveSampleCount(std::max(sp.r0, sp.r1), std::abs(sp.az1 - sp.az0), wireRadius);
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
TopoDS_Edge primEdge(const Primitive& pr, double wireRadius) {
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
                TopoDS_Edge e;
                if (!sp.blend) {
                    Handle(Geom2d_TrimmedCurve) seg2d =
                        GCE2d_MakeSegment(gp_Pnt2d(sp.az0, v0),
                                          gp_Pnt2d(sp.az1, v0 + v1)).Value();
                    e = BRepBuilderAPI_MakeEdge(seg2d, surf).Edge();
                } else {
                    // Cosine blend of V over U — end tangents purely azimuthal, so the
                    // junctions into the neighbouring on-station geometry are exact.
                    int n = curveSampleCount(std::max(sp.r0, sp.r1),
                                             std::abs(sp.az1 - sp.az0), wireRadius);
                    Handle(TColgp_HArray1OfPnt2d) arr =
                        new TColgp_HArray1OfPnt2d(1, n);
                    for (int i = 0; i < n; ++i) {
                        double t = static_cast<double>(i) / (n - 1);
                        double f = 0.5 * (1.0 - std::cos(kPi * t));
                        arr->SetValue(i + 1, gp_Pnt2d(sp.az0 + (sp.az1 - sp.az0) * t,
                                                      v0 + v1 * f));
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
        if (sp.blend) {
            // Flat blended spiral (radius varies at constant height): a 2D spiral in
            // the horizontal plane's parameter space, interpolated with exact azimuthal
            // end tangents.
            try {
                gp_Ax3 frame(gp_Pnt(sp.cx, sp.y0, sp.cz), gp_Dir(0, 1, 0),
                             gp_Dir(1, 0, 0));
                Handle(Geom_Plane) plane = new Geom_Plane(frame);
                // Plane P(U,V) = O + U XDir + V YDir with YDir = -Z, so the azPointC
                // trace maps to (U, V) = (r cos az, r sin az).
                int n = curveSampleCount(std::max(sp.r0, sp.r1),
                                         std::abs(sp.az1 - sp.az0), wireRadius);
                Handle(TColgp_HArray1OfPnt2d) arr = new TColgp_HArray1OfPnt2d(1, n);
                for (int i = 0; i < n; ++i) {
                    double t = static_cast<double>(i) / (n - 1);
                    double f = 0.5 * (1.0 - std::cos(kPi * t));
                    double az = sp.az0 + (sp.az1 - sp.az0) * t;
                    double r = sp.r0 + dr * f;
                    arr->SetValue(i + 1, gp_Pnt2d(r * std::cos(az), r * std::sin(az)));
                }
                gp_Vec2d ta(-std::sin(sp.az0), std::cos(sp.az0));
                gp_Vec2d tb(-std::sin(sp.az1), std::cos(sp.az1));
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
        // Flat LINEAR spiral: falls through to the sampled BSpline below.
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
    // Flat varying-radius spirals at exactly constant height: approximated BSpline
    // through the sampled centreline.
    auto pts = samplePrim(pr, wireRadius);
    if (pts.size() < 2) return TopoDS_Edge();
    try {
        TColgp_Array1OfPnt arr(1, static_cast<Standard_Integer>(pts.size()));
        for (size_t i = 0; i < pts.size(); ++i)
            arr.SetValue(static_cast<Standard_Integer>(i + 1), pts[i]);
        Handle(Geom_BSplineCurve) bs = GeomAPI_PointsToBSpline(arr).Curve();
        if (bs.IsNull()) return TopoDS_Edge();
        return BRepBuilderAPI_MakeEdge(bs).Edge();
    } catch (const Standard_Failure&) {
        return TopoDS_Edge();
    }
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
    TopoDS_Edge e = primEdge(pr, r);
    if (e.IsNull()) return {};
    try {
        const auto ends = primEndpoints(pr);
        const gp_Dir tA = primFwdStart(pr, r), tB = primFwdEnd(pr, r);
        BRepBuilderAPI_MakeWire mw;
        gp_Pnt spineStart = ends.first;
        if (overA > 0.0) {
            spineStart = gp_Pnt(ends.first.XYZ() - tA.XYZ() * overA);
            mw.Add(BRepBuilderAPI_MakeEdge(spineStart, ends.first).Edge());
        }
        mw.Add(e);
        if (overB > 0.0)
            mw.Add(BRepBuilderAPI_MakeEdge(ends.second,
                                           gp_Pnt(ends.second.XYZ() + tB.XYZ() * overB)).Edge());
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
        TopoDS_Wire prof = wireProfileWire(spineStart, tA, r, segments);
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
static TopoDS_Shape localMitreTrim(const TopoDS_Shape& solid, const gp_Pnt& P, const gp_Dir& keepDir,
                                   const gp_Dir& stubDir, double r, double grow) {
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
    const double depth = grow + r;
    gp_XYZ u = stubDir.XYZ() - w * stubDir.XYZ().Dot(w);
    if (u.Modulus() > 1e-9) u.Normalize();
    else u = gp_Ax2(P, gp_Dir(w)).XDirection().XYZ();   // stub ~normal to plane: no drift
    const gp_XYZ v = w.Crossed(u);
    const double vHalf = 1.05 * r;
    gp_Pnt corner(P.XYZ() - u * ell - v * vHalf);
    TopoDS_Shape box = BRepPrimAPI_MakeBox(gp_Ax2(corner, gp_Dir(w), gp_Dir(u)),
                                           drift + 2.0 * ell, 2.0 * vHalf, depth).Shape();
    // BOUND THE KNIFE TO THE JUNCTION. A box has straight sides, so on a CURVED piece the tube
    // bends back into its far region and the box bites a rectangular notch out of material that
    // should have stayed -- the square voids beside every dragback. Intersecting the box with a
    // sphere centred on the joint keeps the cut strictly local: nothing further from the corner
    // than the overhang itself can be reached, whatever the piece does downstream.
    try {
        TopoDS_Shape ball =
            BRepPrimAPI_MakeSphere(P, grow + 1.5 * r).Shape();
        BRepAlgoAPI_Common local(box, ball);
        if (local.IsDone() && !local.Shape().IsNull()) box = local.Shape();
    } catch (const Standard_Failure&) {
    }
    try {
        GProp_GProps gpBefore, gpAfter;
        BRepGProp::VolumeProperties(solid, gpBefore);
        BRepAlgoAPI_Cut cut(solid, box);
        if (cut.IsDone() && !cut.Shape().IsNull()) {
            BRepGProp::VolumeProperties(cut.Shape(), gpAfter);
            const double removed = gpBefore.Mass() - gpAfter.Mass();
            // stub bound: the overhang is at most a full-radius tube of length grow+2r, doubled
            // for the tilted-ellipse wedge. More than that = the knife ate distant material.
            const double stubMax = 2.0 * kPi * r * r * (grow + 2.0 * r);
            if (removed <= stubMax && gpAfter.Mass() > 0.0) return cut.Shape();
        }
    } catch (const Standard_Failure&) {
    }
    return solid;
}

TopoDS_Shape assembleWire(const std::vector<const Primitive*>& ptrs, double wireRadius,
                          int segments, CornerStyle corners,
                          std::vector<size_t>* primIndexPerSolid) {
    const double over = 1.3 * wireRadius;   // covers a bisector tilt up to ~50 deg half-angle
    const double tanThresh = 0.05;          // rad (~3 deg): below this a junction is tangent -> no cut
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
                      << " deg" << (fe[i - 1].Angle(fs[i]) > tanThresh ? " CORNER" : "")
                      << " ['" << ptrs[i - 1]->label << "' -> '" << ptrs[i]->label << "']\n";
        }
        // ABT #685 (Alf, 2026-08-14): a lead/turn corner is mitred exactly like the elbow between
        // the two lead runs themselves (the vertical bridge into the radial run out) — one bisector
        // plane, both sides grown to it and sliced, so on a 90 degree corner the cut is at 45 and
        // the two pieces share the face. Tried and rejected on the way here: slicing on the LEAD's
        // face, and leaving both flush with a wire-radius overshoot — the first is a long slanted
        // ellipse, the second leaves the corner's outer side notched.
        if (i > 0 && fe[i - 1].Angle(fs[i]) > tanThresh) {
            gp_Vec s(fe[i - 1].XYZ());
            s += gp_Vec(fs[i].XYZ());
            if (s.Magnitude() > 1e-9) { nS = gp_Dir(s); bentS = true; }
        }
        bool bentE = false;
        gp_Dir nE = fe[i];
        if (i + 1 < n && fe[i].Angle(fs[i + 1]) > tanThresh) {
            gp_Vec s(fe[i].XYZ());
            s += gp_Vec(fs[i + 1].XYZ());
            if (s.Magnitude() > 1e-9) { nE = gp_Dir(s); bentE = true; }
        }
        // Growth per end: corners grow by over+mismatch (then trimmed on the mitre plane);
        // a TANGENT junction with an endpoint gap is bridged by growing the EARLIER prim's end
        // flush forward (no cut) -- one side only, so the bridge is never doubled.
        const double overS = bentS ? over + dpS : 0.0;
        const double overE = bentE ? over + dpE : (dpE > 1e-9 ? dpE : 0.0);
        TopoDS_Shape solid = rawGrownSolid(*ptrs[i], wireRadius, overS, overE, segments);
        if (solid.IsNull()) {  // ARC clamp (near-full revolve): fall back to a flush, uncut tube
            solid = rawGrownSolid(*ptrs[i], wireRadius, 0.0, 0.0, segments);
            bentS = bentE = false;
        }
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
        if (bentS) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i - 1]).second.XYZ() +
                                  primEndpoints(*ptrs[i]).first.XYZ()));
            solid = localMitreTrim(solid, J, nS, gp_Dir(fs[i].XYZ() * -1.0), wireRadius, overS);
            ++nCut;
        }
        if (bentE) {
            const gp_Pnt J(0.5 * (primEndpoints(*ptrs[i]).second.XYZ() +
                                  primEndpoints(*ptrs[i + 1]).first.XYZ()));
            solid = localMitreTrim(solid, J, gp_Dir(nE.XYZ() * -1.0), fe[i], wireRadius, overE);
            ++nCut;
        }
        // STEP round-trip robustness: the boolean's cut curve on a torus can carry an edge tolerance
        // that a STEP export degrades into an invalid solid (valid in memory, invalid on reload).
        // ShapeFix tightens edges/tolerances in place so the written solid survives the round-trip.
        if (!solid.IsNull() && (bentS || bentE)) {
            try {
                ShapeFix_Shape fix(solid);
                fix.Perform();
                if (!fix.Shape().IsNull()) solid = fix.Shape();
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
