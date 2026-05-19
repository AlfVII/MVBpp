#pragma once

#include "mvb/NamedShape.h"
#include <TopoDS_Shape.hxx>
#include <gp_Pln.hxx>
#include <string>
#include <vector>

namespace mvb {

enum class SectionPlane {
    XY,  // z = const
    XZ,  // y = const
    YZ,  // x = const
};

// Parse "XY" / "XZ" / "YZ" (case-insensitive). Throws on anything else.
SectionPlane parseSectionPlane(const std::string& s);

class SectionBuilder {
public:
    // Cuts `solid` with the given plane and returns a compound of the
    // resulting section edges (1D wires lying in the plane).
    static TopoDS_Shape sectionCore(const TopoDS_Shape& solid, SectionPlane plane);

    // Cuts each input named shape with `plane` shifted by `offset` (metres
    // along the plane normal) and returns a list of (face-or-wire-compound,
    // name) pairs — one per input shape that produced any section. The
    // section result is a compound of edges; closed wires are upgraded to
    // planar TopoDS_Faces where possible (suitable for ElmerFEM/Gmsh
    // surface meshing). Returns NamedShape entries with the same name as
    // the input.
    static std::vector<NamedShape> cut2DFaces(const std::vector<NamedShape>& shapes,
                                              SectionPlane plane,
                                              double offset = 0.0);

    // Renders the section edges as an SVG document. The 3D edges are
    // projected onto the section plane (Z dropped for XY, Y dropped for XZ,
    // X dropped for YZ). Returns the SVG XML as a string.
    static std::string edgesToSvg(const TopoDS_Shape& edges, SectionPlane plane,
                                   double width_px = 800.0,
                                   double margin_px = 40.0,
                                   double stroke_width = 0.5,
                                   const std::string& stroke_color = "#000000");

    // Convenience: section the core and write SVG to disk.
    static void writeSectionSvg(const TopoDS_Shape& solid, SectionPlane plane,
                                  const std::string& outputPath);
};

} // namespace mvb
