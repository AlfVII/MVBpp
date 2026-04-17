#include "mvb/shapes/ShapeT.h"
#include "mvb/Utils.h"
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>

namespace mvb {
namespace shapes {

TopoDS_Face ShapeT::buildProfile(const std::map<std::string, double>& dims) const {
    double a = 0.0; // outer diameter
    double b = 0.0; // inner diameter
    auto itA = dims.find("A");
    auto itB = dims.find("B");
    if (itA != dims.end()) a = itA->second / 2.0;
    if (itB != dims.end()) b = itB->second / 2.0;

    TopoDS_Wire outer = build_polygon_circle(a, DEFAULT_CORE_POLYGON_SEGMENTS);
    TopoDS_Wire inner = build_polygon_circle(b, DEFAULT_CORE_POLYGON_SEGMENTS);

    BRepBuilderAPI_MakeFace face(outer);
    face.Add(inner);
    return face.Face();
}

TopoDS_Shape ShapeT::applyExtras(const std::map<std::string, double>& dims,
                                 const TopoDS_Shape& piece) const {
    double c = 0.0;
    auto it = dims.find("C");
    if (it != dims.end()) c = it->second;
    return translate_shape(piece, 0.0, 0.0, -c / 2.0);
}

} // namespace shapes
} // namespace mvb
