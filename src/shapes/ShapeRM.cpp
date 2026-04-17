#include "mvb/shapes/ShapeRM.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <gp_Pnt.hxx>
#include <cmath>

namespace mvb {
namespace shapes {

// Mirrors MVB.js RMShape.getShapeBase + PShape.getNegativeWindingWindow.
// Profile is the RM polygonal cross-section (no central hole), extruded by B.
// Then a donut-shaped winding window (E/2 outer, F/2 inner, height D) is cut
// at the top. Optional central H hole through the full height.
TopoDS_Shape ShapeRM::buildPiece(const MAS::CoreShape& shapeData) const {
    auto dimsOpt = shapeData.get_dimensions();
    if (!dimsOpt) return TopoDS_Shape();
    auto dims = flatten_dimensions(*dimsOpt);

    double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0, h = 0.0, j = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second / 2.0;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second / 2.0;
    it = dims.find("F"); if (it != dims.end()) f = it->second / 2.0;
    it = dims.find("G"); if (it != dims.end()) g = it->second / 2.0;
    it = dims.find("H"); if (it != dims.end()) h = it->second / 2.0;
    it = dims.find("J"); if (it != dims.end()) j = it->second;
    if (b == 0.0) return TopoDS_Shape();

    std::string familySubtype = shapeData.get_family_subtype().value_or("1");

    double p = std::sqrt(2.0) * j - 2.0 * a;

    double zv = 0.0;
    if (e > 0.0 && (g * 2.0) <= (e * 2.0)) {
        double alpha = std::asin((g * 2.0) / (e * 2.0));
        zv = e * std::cos(alpha);
    }

    double t = 0.0, n = 1.0, r = 0.0, s = 0.0;
    if (familySubtype == "1") {
        t = 0.0;
        n = (zv - c) / g;
    } else if (familySubtype == "2") {
        if (f > 0.0 && std::abs(c / f) <= 1.0) {
            t = f * std::sin(std::acos(c / f));
        }
        n = (zv - c) / g;
    } else if (familySubtype == "3") {
        if (e > 0.0 && std::abs((g * 2.0) / (e * 2.0)) <= 1.0) {
            t = c - e * std::cos(std::asin((g * 2.0) / (e * 2.0))) + g;
        }
        n = (zv - c) / g;
    } else { // subtype 4
        t = 0.0; n = 1.0;
    }
    r = (a + p / 2.0 - c + n * t) / (n + 1.0);
    s = n * r + c;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(a, -p / 2.0, 0.0));
    poly.Add(gp_Pnt(a, p / 2.0, 0.0));
    poly.Add(gp_Pnt(r, s, 0.0));
    if (familySubtype != "4") {
        poly.Add(gp_Pnt(t, c, 0.0));
        poly.Add(gp_Pnt(-t, c, 0.0));
    }
    poly.Add(gp_Pnt(-r, s, 0.0));
    poly.Add(gp_Pnt(-a, p / 2.0, 0.0));
    poly.Add(gp_Pnt(-a, -p / 2.0, 0.0));
    poly.Add(gp_Pnt(-r, -s, 0.0));
    if (familySubtype != "4") {
        poly.Add(gp_Pnt(-t, -c, 0.0));
        poly.Add(gp_Pnt(t, -c, 0.0));
    }
    poly.Add(gp_Pnt(r, -s, 0.0));
    poly.Close();

    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire()).Face();
    TopoDS_Shape piece = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, b)).Shape();

    // Donut winding window: outer cylinder E/2, inner column F/2, height D, top at z=B
    if (d > 0.0 && e > 0.0 && f > 0.0 && f < e) {
        TopoDS_Shape outerCyl = BRepPrimAPI_MakeCylinder(e, d).Shape();
        outerCyl = translate_shape(outerCyl, 0.0, 0.0, b - d);
        TopoDS_Shape innerCyl = BRepPrimAPI_MakeCylinder(f, d).Shape();
        innerCyl = translate_shape(innerCyl, 0.0, 0.0, b - d);
        BRepAlgoAPI_Cut donut(outerCyl, innerCyl);
        if (donut.IsDone()) {
            BRepAlgoAPI_Cut cutter(piece, donut.Shape());
            if (cutter.IsDone()) piece = cutter.Shape();
        }
    }

    // Optional central H hole through the full height
    if (h > 0.0) {
        TopoDS_Shape hole = BRepPrimAPI_MakeCylinder(h, b).Shape();
        BRepAlgoAPI_Cut cutter(piece, hole);
        if (cutter.IsDone()) piece = cutter.Shape();
    }

    piece = translate_shape(piece, 0.0, 0.0, -b);
    piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    return piece;
}

} // namespace shapes
} // namespace mvb
