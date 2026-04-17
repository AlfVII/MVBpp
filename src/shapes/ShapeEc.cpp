#include "mvb/shapes/ShapeEc.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// EC: E-shape with two side dents (between ±(t+s) and ±a, height ±s).
// Mirrors MVB.js ECShape.getShapeBase.
TopoDS_Face ShapeEc::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0, c = 0.0, t = 0.0, sv = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("C"); if (it != dims.end()) c = it->second / 2.0;
    it = dims.find("T"); if (it != dims.end()) t = it->second / 2.0;
    it = dims.find("s"); if (it != dims.end()) sv = it->second / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a,  c, 0.0));
    poly.Add(gp_Pnt( a,  c, 0.0));
    poly.Add(gp_Pnt( a,  sv, 0.0));
    poly.Add(gp_Pnt( t + sv,  sv, 0.0));
    poly.Add(gp_Pnt( t + sv, -sv, 0.0));
    poly.Add(gp_Pnt( a, -sv, 0.0));
    poly.Add(gp_Pnt( a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -sv, 0.0));
    poly.Add(gp_Pnt(-(t + sv), -sv, 0.0));
    poly.Add(gp_Pnt(-(t + sv),  sv, 0.0));
    poly.Add(gp_Pnt(-a,  sv, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

} // namespace shapes
} // namespace mvb
