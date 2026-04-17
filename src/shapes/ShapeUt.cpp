#include "mvb/shapes/ShapeUt.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace mvb {
namespace shapes {

// Mirrors MVB Python CadQueryBuilder.Ut: rectangle A×C profile, box WW,
// then two asymmetric columns added on top.
TopoDS_Shape ShapeUt::buildPiece(const MAS::CoreShape& shapeData) const {
    auto dimsOpt = shapeData.get_dimensions();
    if (!dimsOpt) return TopoDS_Shape();
    auto dims = flatten_dimensions(*dimsOpt);

    double a = 0.0, b = 0.0, c = 0.0, d = 0.0, e = 0.0, fdim = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second;
    it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) fdim = it->second;
    if (b == 0.0) return TopoDS_Shape();

    double ha = a / 2.0, hc = c / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-ha,  hc, 0.0));
    poly.Add(gp_Pnt( ha,  hc, 0.0));
    poly.Add(gp_Pnt( ha, -hc, 0.0));
    poly.Add(gp_Pnt(-ha, -hc, 0.0));
    poly.Close();
    TopoDS_Face face = BRepBuilderAPI_MakeFace(poly.Wire()).Face();
    TopoDS_Shape piece = BRepPrimAPI_MakePrism(face, gp_Vec(0, 0, b)).Shape();

    // Cut WW: box(2A, 2C, D) translated to z = B/2 (centered at z=B/2)
    if (d > 0.0) {
        gp_Pnt cn(-a, -c, b / 2.0 - d / 2.0);
        TopoDS_Shape ww = BRepPrimAPI_MakeBox(cn, 2.0 * a, 2.0 * c, d).Shape();
        BRepAlgoAPI_Cut cut(piece, ww);
        if (cut.IsDone()) piece = cut.Shape();
    }

    // Top column: F × C × D at x = -A/2 + F/2, z = B/2
    if (fdim > 0.0 && d > 0.0) {
        double xc = -a / 2.0 + fdim / 2.0;
        gp_Pnt cn(xc - fdim / 2.0, -hc, b / 2.0 - d / 2.0);
        TopoDS_Shape topCol = BRepPrimAPI_MakeBox(cn, fdim, c, d).Shape();
        BRepAlgoAPI_Fuse fuser(piece, topCol);
        if (fuser.IsDone()) piece = fuser.Shape();
    }

    // Bottom column: (A-E-F) × C × D at x = A/2 - bcw/2, z = B/2
    double bcw = a - e - fdim;
    if (bcw > 0.0 && d > 0.0) {
        double xc = a / 2.0 - bcw / 2.0;
        gp_Pnt cn(xc - bcw / 2.0, -hc, b / 2.0 - d / 2.0);
        TopoDS_Shape botCol = BRepPrimAPI_MakeBox(cn, bcw, c, d).Shape();
        BRepAlgoAPI_Fuse fuser(piece, botCol);
        if (fuser.IsDone()) piece = fuser.Shape();
    }

    piece = translate_shape(piece, 0.0, 0.0, -b);
    piece = rotate_shape(piece, -std::numbers::pi / 2.0, 0.0, 0.0);
    return piece;
}

} // namespace shapes
} // namespace mvb
