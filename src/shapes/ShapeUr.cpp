#include "mvb/shapes/ShapeUr.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace mvb {
namespace shapes {

// Mirrors MVB.js URShape: 4 subtypes with rectangular profile + cylindrical
// winding/lateral columns.
TopoDS_Shape ShapeUr::buildPiece(const MAS::CoreShape& shapeData) const {
    auto dimsOpt = shapeData.get_dimensions();
    if (!dimsOpt) return TopoDS_Shape();
    auto dims = flatten_dimensions(*dimsOpt);

    double aFull = 0.0, b = 0.0, cFull = 0.0, d = 0.0, fFull = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) aFull = it->second;
    it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) cFull = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("F"); if (it != dims.end()) fFull = it->second;

    if (b == 0.0) return TopoDS_Shape();

    std::string subtype = shapeData.get_family_subtype().value_or("1");
    double c = cFull / 2.0;

    double windingColumnWidth = 0.0;
    bool isHalfRect = (subtype == "1" || subtype == "3");

    if (subtype == "1" || subtype == "2") windingColumnWidth = cFull;
    else /* "3" or "4" */                windingColumnWidth = fFull;

    double leftA;
    if (isHalfRect) leftA = aFull - windingColumnWidth / 2.0;
    else            leftA = aFull - windingColumnWidth;

    BRepBuilderAPI_MakePolygon poly;
    if (isHalfRect) {
        poly.Add(gp_Pnt(0.0, c, 0.0));
        poly.Add(gp_Pnt(-leftA, c, 0.0));
        poly.Add(gp_Pnt(-leftA, -c, 0.0));
        poly.Add(gp_Pnt(0.0, -c, 0.0));
        poly.Add(gp_Pnt(windingColumnWidth / 2.0, -c, 0.0));
        poly.Add(gp_Pnt(windingColumnWidth / 2.0, c, 0.0));
    } else {
        poly.Add(gp_Pnt(0.0, c, 0.0));
        poly.Add(gp_Pnt(-leftA, c, 0.0));
        poly.Add(gp_Pnt(-leftA, -c, 0.0));
        poly.Add(gp_Pnt(0.0, -c, 0.0));
    }
    poly.Close();

    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire()).Face();
    TopoDS_Shape piece = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, b)).Shape();

    // Winding window cutout: box(2A, 2C, D) at top
    if (d > 0.0 && aFull > 0.0 && cFull > 0.0) {
        gp_Pnt corner(-aFull, -cFull, b - d);
        TopoDS_Shape ww = BRepPrimAPI_MakeBox(corner, 2.0 * aFull, 2.0 * cFull, d).Shape();
        BRepAlgoAPI_Cut cut(piece, ww);
        if (cut.IsDone()) piece = cut.Shape();
    }

    // Add winding & lateral columns
    auto fuse = [&](TopoDS_Shape& base, const TopoDS_Shape& add) {
        if (add.IsNull()) return;
        BRepAlgoAPI_Fuse f(base, add);
        if (f.IsDone()) base = f.Shape();
    };

    if (subtype == "1") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(cFull / 2.0, d).Shape();
        wcCol = translate_shape(wcCol, 0.0, 0.0, b - d);
        fuse(piece, wcCol);
        // lateral box H × C × D at -(A - wcw/2 - H/2)
        double h = 0.0;
        auto itH = dims.find("H"); if (itH != dims.end()) h = itH->second;
        if (h > 0.0) {
            double xOff = -(aFull - cFull / 2.0 - h / 2.0);
            gp_Pnt cn(xOff - h / 2.0, -cFull / 2.0, b - d);
            TopoDS_Shape latBox = BRepPrimAPI_MakeBox(cn, h, cFull, d).Shape();
            fuse(piece, latBox);
        }
    } else if (subtype == "2") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(cFull / 2.0, b).Shape();
        fuse(piece, wcCol);
        TopoDS_Shape latCol = BRepPrimAPI_MakeCylinder(cFull / 2.0, b).Shape();
        latCol = translate_shape(latCol, -(aFull - cFull), 0.0, 0.0);
        fuse(piece, latCol);
    } else if (subtype == "3") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(fFull / 2.0, d).Shape();
        wcCol = translate_shape(wcCol, 0.0, 0.0, b - d);
        fuse(piece, wcCol);
        double h = 0.0;
        auto itH = dims.find("H"); if (itH != dims.end()) h = itH->second;
        if (h > 0.0) {
            double xOff = -(aFull - fFull / 2.0 - h / 2.0);
            gp_Pnt cn(xOff - h / 2.0, -cFull / 2.0, b - d);
            TopoDS_Shape latBox = BRepPrimAPI_MakeBox(cn, h, cFull, d).Shape();
            fuse(piece, latBox);
        }
    } else if (subtype == "4") {
        TopoDS_Shape wcCol = BRepPrimAPI_MakeCylinder(fFull / 2.0, b).Shape();
        fuse(piece, wcCol);
        TopoDS_Shape latCol = BRepPrimAPI_MakeCylinder(fFull / 2.0, b).Shape();
        latCol = translate_shape(latCol, -(aFull - fFull), 0.0, 0.0);
        fuse(piece, latCol);
    }

    piece = translate_shape(piece, 0.0, 0.0, -b);
    piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    return piece;
}

} // namespace shapes
} // namespace mvb
