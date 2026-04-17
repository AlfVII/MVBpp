#include "mvb/shapes/ShapeE.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeE::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0, cFull = 0.0, f = 0.0, f2 = 0.0, q = 0.0, k = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("C"); if (it != dims.end()) cFull = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("F2"); if (it != dims.end()) f2 = it->second;
    it = dims.find("q"); if (it != dims.end()) q = it->second;
    it = dims.find("K"); if (it != dims.end()) k = it->second;

    // EFD variant: when q > 0 and F2 > 0, the profile is asymmetric with a
    // triangular dent notch at the top-center (and an optional mini-indent
    // at the bottom-center when K > 0). Matches MVB Python Efd.get_shape_base.
    bool isEfd = (q > 0.0 && f2 > 0.0 && f > 0.0);
    if (isEfd) {
        double topC = cFull - k - f2 / 2.0;
        double bottomC = k + f2 / 2.0;
        double dentHeight = cFull * 2.0 / 5.0;
        double dentTopW = f / 2.0;
        double dentBottomW = f / 2.0 - q;

        BRepBuilderAPI_MakePolygon poly;
        poly.Add(gp_Pnt(-a, topC, 0.0));
        poly.Add(gp_Pnt(-dentTopW, topC, 0.0));
        poly.Add(gp_Pnt(-dentBottomW, topC - dentHeight, 0.0));
        poly.Add(gp_Pnt(dentBottomW, topC - dentHeight, 0.0));
        poly.Add(gp_Pnt(dentTopW, topC, 0.0));
        poly.Add(gp_Pnt(a, topC, 0.0));
        poly.Add(gp_Pnt(a, -bottomC, 0.0));
        if (k > 0.0) {
            double miniW = f / 2.0 - q;
            if (miniW < 0.0) miniW = 0.0;
            poly.Add(gp_Pnt(miniW, -bottomC, 0.0));
            poly.Add(gp_Pnt(miniW, -bottomC + k, 0.0));
            poly.Add(gp_Pnt(-miniW, -bottomC + k, 0.0));
            poly.Add(gp_Pnt(-miniW, -bottomC, 0.0));
        }
        poly.Add(gp_Pnt(-a, -bottomC, 0.0));
        poly.Close();
        BRepBuilderAPI_MakeWire wire(poly.Wire());
        BRepBuilderAPI_MakeFace face(wire.Wire());
        return face.Face();
    }

    double c = cFull / 2.0;
    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a, c, 0.0));
    poly.Add(gp_Pnt(a, c, 0.0));
    poly.Add(gp_Pnt(a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Close();

    BRepBuilderAPI_MakeWire wire(poly.Wire());
    BRepBuilderAPI_MakeFace face(wire.Wire());
    return face.Face();
}

TopoDS_Shape ShapeE::buildWindingWindow(const std::map<std::string, double>& dims) const {
    double b = 0.0, d = 0.0, e = 0.0, f = 0.0, c = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("D"); if (it != dims.end()) d = it->second;
    it = dims.find("E"); if (it != dims.end()) e = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("C"); if (it != dims.end()) c = it->second;

    if (b == 0.0 || d == 0.0 || e == 0.0) return TopoDS_Shape();

    // Winding window box: E × 2C × D, Y-centered. The 2C height (vs C)
    // ensures the cut fully covers asymmetric profiles (EFD with K>0 has
    // top_c ≠ bottom_c). For symmetric profiles the extra Y outside the
    // profile is clipped away by the profile intersection, so it's safe.
    gp_Pnt wwCorner(-e / 2.0, -c, b - d);
    TopoDS_Shape wwBox = BRepPrimAPI_MakeBox(wwCorner, e, 2.0 * c, d).Shape();

    double f2 = f;
    it = dims.find("F2");
    if (it != dims.end()) f2 = it->second;
    bool oblong = std::abs(f2 - f) > 0.0001;

    TopoDS_Shape column;
    if (oblong && f2 > f) {
        // Stadium-shaped central column for EL cores
        double half_width = f / 2.0;      // X direction, semicircle radius
        double half_depth = f2 / 2.0;     // Y direction
        double rect_half_length = half_depth - half_width;
        double zCenter = b - d / 2.0;

        if (rect_half_length <= 0.0) {
            // Actually round
            column = BRepPrimAPI_MakeCylinder(half_width, d).Shape();
            // MakeCylinder is base-at-zero; center it at zCenter
            column = translate_shape(column, 0.0, 0.0, zCenter - d / 2.0);
        } else {
            gp_Pnt rectCorner(-half_width, -rect_half_length, b - d);
            TopoDS_Shape centerRect = BRepPrimAPI_MakeBox(rectCorner, f, rect_half_length * 2.0, d).Shape();

            TopoDS_Shape topCyl = BRepPrimAPI_MakeCylinder(half_width, d).Shape();
            topCyl = translate_shape(topCyl, 0.0, rect_half_length, zCenter - d / 2.0);

            TopoDS_Shape bottomCyl = BRepPrimAPI_MakeCylinder(half_width, d).Shape();
            bottomCyl = translate_shape(bottomCyl, 0.0, -rect_half_length, zCenter - d / 2.0);

            BRepAlgoAPI_Fuse f1(centerRect, topCyl);
            if (!f1.IsDone()) return TopoDS_Shape();
            TopoDS_Shape fused = f1.Shape();
            BRepAlgoAPI_Fuse f2(fused, bottomCyl);
            if (f2.IsDone()) fused = f2.Shape();
            column = fused;
        }
    } else if (oblong && f2 < f) {
        // EFD: leave the WW cut as a plain E × 2C × D box (no column carved
        // out). The column is added back separately in applyExtras so it
        // survives the profile dent intact at full F × F2 × B size.
        return wwBox;
    } else {
        // Rectangular central column: F x C x D, centered at same Z
        gp_Pnt ccCorner(-f / 2.0, -c / 2.0, b - d);
        column = BRepPrimAPI_MakeBox(ccCorner, f, c, d).Shape();
    }

    BRepAlgoAPI_Cut cutter(wwBox, column);
    return cutter.IsDone() ? cutter.Shape() : TopoDS_Shape();
}

TopoDS_Shape ShapeE::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    double b = 0.0, f = 0.0, f2 = 0.0, q = 0.0;
    auto it = dims.find("B"); if (it != dims.end()) b = it->second;
    it = dims.find("F"); if (it != dims.end()) f = it->second;
    it = dims.find("F2"); if (it != dims.end()) f2 = it->second;
    it = dims.find("q"); if (it != dims.end()) q = it->second;

    TopoDS_Shape result = piece;
    // EFD (F2 < F, both > 0): add central column with chamfered (octagonal)
    // corners — F × F2 with q-sized 45° chamfers at each corner. Matches
    // MVB.js EFDShape.getShapeExtras.
    if (f > 0.0 && f2 > 0.0 && f2 + 1e-9 < f) {
        double hf = f / 2.0;
        double hf2 = f2 / 2.0;
        TopoDS_Shape column;
        if (q > 0.0) {
            BRepBuilderAPI_MakePolygon poly;
            poly.Add(gp_Pnt( hf - q,  hf2,    0.0));
            poly.Add(gp_Pnt(-hf + q,  hf2,    0.0));
            poly.Add(gp_Pnt(-hf,      hf2 - q, 0.0));
            poly.Add(gp_Pnt(-hf,     -hf2 + q, 0.0));
            poly.Add(gp_Pnt(-hf + q, -hf2,    0.0));
            poly.Add(gp_Pnt( hf - q, -hf2,    0.0));
            poly.Add(gp_Pnt( hf,     -hf2 + q, 0.0));
            poly.Add(gp_Pnt( hf,      hf2 - q, 0.0));
            poly.Close();
            BRepBuilderAPI_MakeFace face(poly.Wire());
            column = BRepPrimAPI_MakePrism(face.Face(), gp_Vec(0, 0, b)).Shape();
        } else {
            gp_Pnt cc(-hf, -hf2, 0.0);
            column = BRepPrimAPI_MakeBox(cc, f, f2, b).Shape();
        }
        BRepAlgoAPI_Fuse fuser(result, column);
        if (fuser.IsDone() && !fuser.Shape().IsNull()) {
            result = fuser.Shape();
        }
    }
    return translate_shape(result, 0.0, 0.0, -b);
}

} // namespace shapes
} // namespace mvb
