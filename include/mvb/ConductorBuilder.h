#pragma once

#include "MAS.hpp"
#include "mvb/NamedShape.h"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
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
// Milestone 1 scope: concentric ROUND columns with round/litz wire. Rectangular/oblong
// columns and toroids throw (milestones 2 and 3).
class ConductorBuilder {
public:
    struct Options {
        Options() {}
        int  wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS;
        // true  -> conductor swept at the OUTER (insulation) footprint (visualisation)
        // false -> swept at the CONDUCTING (copper) footprint (FEM)
        bool paintCoating = true;
    };

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
