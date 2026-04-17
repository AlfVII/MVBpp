#include "mvb/shapes/ShapeLp.h"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <gp_Pnt.hxx>

namespace mvb {
namespace shapes {

// LP shares EP's winding window topology (donut + top G-slot + bottom E-rect),
// but uses a centered A × C rectangular profile (ER-like).
TopoDS_Face ShapeLp::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0, c = 0.0;
    auto it = dims.find("A"); if (it != dims.end()) a = it->second / 2.0;
    it = dims.find("C"); if (it != dims.end()) c = it->second / 2.0;

    BRepBuilderAPI_MakePolygon poly;
    poly.Add(gp_Pnt(-a,  c, 0.0));
    poly.Add(gp_Pnt( a,  c, 0.0));
    poly.Add(gp_Pnt( a, -c, 0.0));
    poly.Add(gp_Pnt(-a, -c, 0.0));
    poly.Close();
    return BRepBuilderAPI_MakeFace(poly.Wire()).Face();
}

} // namespace shapes
} // namespace mvb
