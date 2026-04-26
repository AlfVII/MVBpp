#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <map>
#include <string>
#include <memory>

namespace mvb {
namespace shapes {

class ShapeBuilder {
public:
    virtual ~ShapeBuilder() = default;

    // Polygon segment count used by derived builders when approximating
    // circles/cylinders. Set by the factory so subclasses can consult it
    // via member access without changing their virtual signatures.
    void setCorePolygonSegments(int segments) { m_corePolygonSegments = segments; }
    int  corePolygonSegments() const { return m_corePolygonSegments; }

    // Build a single core piece (includes intrinsic family rotation)
    virtual TopoDS_Shape buildPiece(const MAS::CoreShape& shapeData) const;

    // Apply machining (gap) to an already-rotated piece
    virtual TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                        const MAS::Machining& machining,
                                        const std::map<std::string, double>& dims) const;

protected:
    // Subclasses implement the 2D profile in the XY plane
    virtual TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const = 0;

    // Subclasses implement the winding window cutout (default: none)
    virtual TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const;

    // Optional post-processing before intrinsic rotation (default: identity)
    virtual TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                                     const TopoDS_Shape& piece) const;

    // True for toroidal families (no -90° X rotation)
    virtual bool isToroidal() const { return false; }

    // Helper: extrude a face along Z by height
    static TopoDS_Shape extrude(const TopoDS_Face& face, double height);

    // Helper: create a rectangular box centered at origin
    static TopoDS_Shape makeBox(double x, double y, double z);

    // Helper: create a cylinder along Z axis, centered at origin
    static TopoDS_Shape makeCylinder(double height, double radius, int segments);

    int m_corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS;
};

// Factory — `corePolygonSegments` is applied to the returned builder.
std::unique_ptr<ShapeBuilder> createShapeBuilder(MAS::CoreShapeFamily family,
                                                  const std::string& subtype = "",
                                                  int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS);

} // namespace shapes
} // namespace mvb
