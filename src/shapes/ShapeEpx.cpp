#include "mvb/shapes/ShapeEpx.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// EPX (EP variant with elongated/stadium central column).
// Mirrors MVB.js EPXShape.
TopoDS_Face ShapeEpx::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0, c = 0.0, k = 0.0, f = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("K"); if (it != dims.end()) k = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;

    double columnLength = k + f / 2.0;
    double topC = c - columnLength / 2.0;
    double bottomC = columnLength / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, topC, 0.0));
    poly.Add(gp_Pnt( a, topC, 0.0));
    poly.Add(gp_Pnt( a, -bottomC, 0.0));
    poly.Add(gp_Pnt(-a, -bottomC, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

TopoDS_Shape ShapeEpx::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double b = 0.0, c = 0.0, d = 0.0, e = 0.0, f = 0.0, g = 0.0, k = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("G"); if (it != dims.end()) g = it->second;
    it = dims.find("K"); if (it != dims.end()) k = it->second;
    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    double zBase = b - d;
    double rectangularPartWidth = k - f / 2.0;
    double columnWidth = k + f / 2.0;

    // Stadium-shaped central column: rectangle F × rectangularPartWidth + cylinders at ±rpw/2 in Y
    TopoDS_Shape centralCol;
    {
        gp_Pnt cn(-f / 2.0, -rectangularPartWidth / 2.0, zBase);
        TopoDS_Shape rect = BRepPrimAPI_MakeBox(cn, f, rectangularPartWidth, d).Shape();
        TopoDS_Shape topCyl = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
        topCyl = translate_shape(topCyl, 0.0, rectangularPartWidth / 2.0, zBase);
        TopoDS_Shape botCyl = BRepPrimAPI_MakeCylinder(f / 2.0, d).Shape();
        botCyl = translate_shape(botCyl, 0.0, -rectangularPartWidth / 2.0, zBase);
        BRepAlgoAPI_Fuse f1(rect, topCyl);
        centralCol = f1.IsDone() ? f1.Shape() : rect;
        BRepAlgoAPI_Fuse f2(centralCol, botCyl);
        if (f2.IsDone()) centralCol = f2.Shape();
    }

    // Outer winding-window cylinder, offset by +rpw/2 in Y
    TopoDS_Shape ww = BRepPrimAPI_MakeCylinder(e / 2.0, d).Shape();
    ww = translate_shape(ww, 0.0, rectangularPartWidth / 2.0, zBase);
    {
        BRepAlgoAPI_Cut cut(ww, centralCol);
        if (cut.IsDone()) ww = cut.Shape();
    }

    // Top G-slot cutout if applicable
    if (g > 0.0) {
        gp_Pnt cn(-g / 2.0, c / 2.0 + columnWidth / 2.0 - c / 2.0, zBase);
        // MVB.js: translate [0, width/2 + columnWidth/2, ...] where width=C
        // so Y center = C/2 + columnWidth/2; box centered Y, so Y range = (C/2+columnWidth/2 - C/2, C/2+columnWidth/2 + C/2) = (columnWidth/2, C+columnWidth/2)
        // Actually box spans Y = (yc - C/2, yc + C/2) where yc = C/2 + columnWidth/2.
        gp_Pnt cn2(-g / 2.0, columnWidth / 2.0, zBase);
        TopoDS_Shape topBox = BRepPrimAPI_MakeBox(cn2, g, c, d).Shape();
        BRepAlgoAPI_Fuse fuser(ww, topBox);
        if (fuser.IsDone()) ww = fuser.Shape();
    }

    // Bottom rectangle E × C × D, offset by -C/2 + rpw/2 in Y, minus central column
    {
        // MVB.js: translate [0, -width/2 + rpw/2, ...] where width=C
        // Y center = -C/2 + rpw/2; box spans Y = (-C + rpw/2, rpw/2)
        gp_Pnt cn(-e / 2.0, -c + rectangularPartWidth / 2.0, zBase);
        TopoDS_Shape bottomBox = BRepPrimAPI_MakeBox(cn, e, c, d).Shape();
        BRepAlgoAPI_Cut cut(bottomBox, centralCol);
        if (cut.IsDone()) bottomBox = cut.Shape();
        BRepAlgoAPI_Fuse fuser(ww, bottomBox);
        if (fuser.IsDone()) ww = fuser.Shape();
    }

    return ww;
}

} // namespace shapes
} // namespace mvb
