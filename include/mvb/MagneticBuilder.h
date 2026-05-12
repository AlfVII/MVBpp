#pragma once

#include "MAS.hpp"
#include "Utils.h"
#include "mvb/NamedShape.h"
#include "mvb/Symmetry.h"
#include <TopoDS_Shape.hxx>
#include <vector>
#include <string>
#include <map>

// Forward declarations for MKF-enriched magnetic overloads
namespace OpenMagnetics { class Magnetic; class Coil; }

namespace mvb {

// Unified configuration for drawMagnetic / drawMagneticToBytes.
// All bindings (Python, WASM) accept this struct so the API surface
// is identical across languages.
struct DrawConfig {
    std::string format                = "step";  // "step" or "stl"
    bool        includeBobbin         = true;
    double      scale                 = 1.0;     // 1000 for mm export
    int         symmetryPlanes        = 0;       // 0=full, 1=half, 2=quarter
    int         wirePolygonSegments    = DEFAULT_WIRE_POLYGON_SEGMENTS;
    int         corePolygonSegments    = DEFAULT_CORE_POLYGON_SEGMENTS;
};

class MagneticBuilder {
public:
    // Build geometry and export to STEP ("step") or STL ("stl").
    //
    // Parameters
    // ----------
    //   outputPath            : directory where magnetic.step / magnetic.stl is written
    //   format                : "step" or "stl"
    //   includeBobbin         : include bobbin geometry
    //   scale                 : uniform scale factor (use 1000 for mm export)
    //   symmetryPlanes        : 0=full, 1=half, 2=quarter domain
    //   wirePolygonSegments   : <=0 = exact torus, >0 = faceted polygon (wire cross-section)
    //   corePolygonSegments   : polygon segments for core cylinders/circles
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const std::string& format = "step",
                             bool includeBobbin = true,
                             double scale = 1.0,
                             int symmetryPlanes = 0,
                             int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                             int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Config-based overload — preferred for bindings.
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& cfg) const;

    // Overload that accepts an already-enriched OpenMagnetics::Magnetic
    // to avoid object-slicing issues with MAS::Magnetic
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const std::string& format = "step",
                             bool includeBobbin = true,
                             double scale = 1.0,
                             int symmetryPlanes = 0,
                             int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                             int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Config-based overload — preferred for bindings.
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& cfg) const;

    std::vector<TopoDS_Shape> buildCore(const MAS::MagneticCore& core,
                                         int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    TopoDS_Shape buildBobbin(const MAS::Coil& coil, const MAS::MagneticCore& core) const;
    TopoDS_Shape buildBobbin(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core) const;
    std::vector<TopoDS_Shape> buildTurns(const MAS::Coil& coil, const MAS::MagneticCore& core,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) const;
    std::vector<TopoDS_Shape> buildTurns(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) const;
    std::vector<TopoDS_Shape> buildTurns(const MAS::Coil& coil, const MAS::MagneticCore& core,
                                          std::vector<std::string>& names,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) const;
    std::vector<TopoDS_Shape> buildTurns(const OpenMagnetics::Coil& coil, const MAS::MagneticCore& core,
                                          std::vector<std::string>& names,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) const;

    // Named-shape overloads. Each returned element carries the logical
    // name (core name / "Turn_<i>" or Turn::get_name / bobbin name) so the
    // identity survives downstream operations (symmetry cut, STEP export
    // with XCAF labels, mesh tagging).
    std::vector<NamedShape> buildCoreNamed(const MAS::MagneticCore& core,
                                            int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    std::vector<NamedShape> buildTurnsNamed(const MAS::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) const;
    std::vector<NamedShape> buildTurnsNamed(const OpenMagnetics::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS) const;
    NamedShape buildBobbinNamed(const MAS::Coil& coil,
                                const MAS::MagneticCore& core) const;
    NamedShape buildBobbinNamed(const OpenMagnetics::Coil& coil,
                                const MAS::MagneticCore& core) const;

    // Assemble all geometry as named shapes (no export). Optionally applies
    // symmetry cuts according to symmetryPlanes (0=full, 1=half, 2=quarter).
    std::vector<NamedShape> buildAllNamed(const MAS::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    std::vector<NamedShape> buildAllNamed(const OpenMagnetics::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
};

} // namespace mvb
