#pragma once

#include <TopoDS_Shape.hxx>
#include <gp_Pln.hxx>
#include <string>
#include <vector>

namespace mvb {

enum class SectionPlane {
    XY,  // z = 0
    XZ,  // y = 0
    YZ,  // x = 0
};

class SectionBuilder {
public:
    // Cuts `solid` with the given plane and returns a compound of the
    // resulting section edges (1D wires lying in the plane).
    static TopoDS_Shape sectionCore(const TopoDS_Shape& solid, SectionPlane plane);

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
