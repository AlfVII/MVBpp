#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <array>
#include <limits>
#include <vector>

namespace OpenMagnetics { class Coil; }

namespace mvb {

// Real-winding geometry: ONE continuous copper body per (winding, parallel), replacing
// TurnBuilder's independent closed loops when DrawConfig::useRealWindingGeometry is on.
//
// Contract (the hard constraint):
//   - Every MKF turn keeps its exact MAS position: each turn is an on-station arc swept at
//     exactly (coordinates[0], coordinates[1]); MVB++ adds ONLY connecting geometry, and
//     only inside the seam sector / terminal routes MKF leaves unspecified.
//   - Electrical order comes from turnsDescription vector order filtered by
//     (turn.winding, turn.parallel) — MKF's documented contract.
//   - The path planner collision-checks every primitive (capsule model) against all
//     conductors and the winding-window bounds, and THROWS on any overlap beyond contact;
//     nothing is ever moved to "make it fit".
//
// Scope: concentric ROUND, RECTANGULAR and OBLONG columns and TOROIDS, with round/litz
// wire. Rectangular/planar/foil wire throws (the lead cross-section orientation through
// the exit bends is undefined until specified).
class ConductorBuilder {
public:
    struct Options {
        Options() {}
        int  wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS;
        // true  -> conductor swept at the OUTER (insulation) footprint (visualisation)
        // false -> swept at the CONDUCTING (copper) footprint (FEM)
        bool paintCoating = true;
        // false (OM / drawing) -> fast per-run compound: overlapping swept pipes, perfect for the
        //   3D viewer, seconds to build.
        // true  (FEM export)   -> the slow one-piece machinery: single continuous body per parallel
        //   where a sweep closes (columns, sparse toroids), a conformal mitre-jointed compound for
        //   dense toroids, fuse for single-turns. Meshable (no self-overlap) but far slower.
        bool femReady = false;
        // Azimuth (radians, azPointC convention: az -> (x = r cos az, z = -r sin az)) the
        // terminal leads should EXIT through -- the core's winding-window OPENING. NAN =
        // keep the default (az = pi/2, the -Z direction). Cores whose lateral legs cover
        // -Z (PQ, EP, ...) MUST redirect the leads or they terminate inside the core plate
        // (measured on 03_buck_pq3230: the lead tip landed exactly on the oblique plate
        // face and tetgen could not recover the contact edge).
        double leadExitAzimuth = std::numeric_limits<double>::quiet_NaN();
        // The CORE solids, for aiming the leads at the true window opening. Column metadata
        // under-describes real cores (a PQ's plates wrap most of the perimeter, leaving only
        // narrow slots); classifying the actual solids around the lead-tip radius finds the
        // genuine free arc. Empty -> no aiming beyond leadExitAzimuth.
        std::vector<TopoDS_Shape> coreObstacles;
    };

    // REAL-PATH POLYLINES: the fully-assembled, collision-checked conductor centrelines
    // (wraps, links, Z end-runs, terminal leads -- everything emitConductor would sweep),
    // sampled per primitive, WITHOUT building any solid. No booleans anywhere: this is the
    // input for implicit (signed-distance / level-set) winding meshing, which sidesteps the
    // OCC weld/fragment failure classes entirely (ABT #490/#491).
    struct PathPolyline {
        std::string name;                 // "<winding> parallel <k>"
        double wireRadius = 0.0;          // round capsule radius (= half min dim for rect)
        bool isRectangular = false;
        double wireWidth = 0.0, wireHeight = 0.0;
        // Per-primitive sampled points (metres). Distance = min over prims of the
        // point-to-polyline capsule distance; per-prim grouping avoids phantom bridges
        // between non-contiguous runs.
        std::vector<std::vector<std::array<double, 3>>> prims;
        std::array<double, 3> end0{}, end1{};   // free ends (terminal port centres)
        std::array<double, 3> dir0{}, dir1{};   // OUTWARD end tangents (port normals)
    };
    static std::vector<PathPolyline> buildAllPaths(const OpenMagnetics::Coil& coil,
                                                   const MAS::CoreBobbinProcessedDescription& bobbin,
                                                   bool isToroidal,
                                                   const Options& opts = {});

    // Builds one solid per (winding, parallel), named "<winding> parallel <k>" (MAS style).
    // `coil` must be an MKF-wound coil (turnsDescription present, positions final);
    // `bobbin` is the resolved+patched processed description (same object MagneticBuilder
    // hands to TurnBuilder).
    static std::vector<NamedShape> buildAll(const MAS::Coil& coil,
                                            const MAS::CoreBobbinProcessedDescription& bobbin,
                                            bool isToroidal,
                                            const Options& opts = {});
    static std::vector<NamedShape> buildAll(const OpenMagnetics::Coil& coil,
                                            const MAS::CoreBobbinProcessedDescription& bobbin,
                                            bool isToroidal,
                                            const Options& opts = {});
};

} // namespace mvb
