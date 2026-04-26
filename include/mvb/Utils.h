#pragma once

#include "MAS.hpp"
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_Circ.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <numbers>

// Forward declaration for MKF Magnetic
namespace OpenMagnetics { class Magnetic; }

namespace mvb {

constexpr int DEFAULT_CORE_POLYGON_SEGMENTS = 16;
constexpr int DEFAULT_WIRE_POLYGON_SEGMENTS = 16;
// Azimuthal segmentation of revolved wire turns. Even with a polygonal
// cross-section, MakeRevol produces analytic surfaces of revolution that
// STEP serialises as CYLINDRICAL/TOROIDAL_SURFACE. Setting this > 0 replaces
// the revolve with a ThruSections loft over N angular slices → fully planar.
// Set to 0 for the legacy MakeRevol behaviour.
constexpr int DEFAULT_WIRE_REVOLUTION_SEGMENTS = 12;

// Extract nominal double from MAS Dimension variant
double flatten_dimension(const MAS::Dimension& dim);

// Extract all nominal dimensions from a shape dimension map
std::map<std::string, double> flatten_dimensions(const std::map<std::string, MAS::Dimension>& dims);

// Build a polygon-approximated circle wire in the XY plane centered at origin
// segments = 0 yields a perfect BRep circle edge
TopoDS_Wire build_polygon_circle(double radius, int segments);

// Build a polygon-approximated cylinder solid along Z axis
// segments = 0 yields exact revolved circle
TopoDS_Shape build_polygon_cylinder(double height, double radius, int segments);

// Build a ring (torus approximation) as a solid by lofting a polygonal
// circular cross-section (cross_segments sides) at `revolution_segments`
// azimuthal stations around the Y axis. Produces only PLANAR faces — STEP
// export contains no CYLINDRICAL / TOROIDAL / SURFACE_OF_REVOLUTION.
// turn_radius = major radius, wire_radius = minor, y = ring plane height.
// revolution_segments <= 0 falls back to a classical MakeRevol (analytic).
TopoDS_Shape build_polygon_ring(double turn_radius, double wire_radius,
                                 double y, int cross_segments,
                                 int revolution_segments);

// Apply a 3D rotation (radians) to a shape about X, Y, Z axes in order
TopoDS_Shape rotate_shape(const TopoDS_Shape& shape, double rx, double ry, double rz);

// Translate a shape
gp_Trsf translation_trsf(double x, double y, double z);
TopoDS_Shape translate_shape(const TopoDS_Shape& shape, double x, double y, double z);

// Get family string from enum
std::string core_shape_family_to_string(MAS::CoreShapeFamily family);

// Enumerate every core shape family that has a builder registered in the
// factory (mirrors MVB.js getSupportedFamilies()).
std::vector<std::string> get_supported_families();

// Preprocess JSON to add missing "nominal" fields to dimension objects
// (MVB test data often has only min/max)
void patch_dimension_nominals(nlohmann::json& j);

// Safely enrich a raw MAS::Magnetic using MKF's magnetic_autocomplete.
// This avoids object-slicing and Coil::wind() crashes for raw MAS files
// (e.g. with "Basic" bobbins) by constructing OpenMagnetics::Coil with
// windInConstructor=false before calling MKF enrichment.
OpenMagnetics::Magnetic magnetic_autocomplete_safe(const MAS::Magnetic& magnetic);
OpenMagnetics::Magnetic magnetic_autocomplete_safe(const nlohmann::json& magneticJson);

// Helpers for cutting bobbin with cores/turns to match Python MVB behavior
bool is_shape_usable(const TopoDS_Shape& shape);
TopoDS_Shape cut_bobbin(const TopoDS_Shape& bobbin, const std::vector<TopoDS_Shape>& cutters);

} // namespace mvb
