#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <utility>
#include <vector>

namespace mvb {

class TurnBuilder {
public:
    // GLOBAL setting: fuse each turn's sub-solid compound (straight segments + corner pieces) into
    // ONE solid. A round/rect turn is built as a COMPOUND of overlapping cylinders/boxes + corner
    // tori/swept-rects -- fast, and perfectly fine for STL/STEP VISUALISATION (the web frontend), so
    // the default is OFF. But the overlapping sub-solids make a tetrahedral MESH fail ("overlapping
    // facets"), so the FEM meshing path must enable fusing to get a single clean conductor per turn.
    // Fusing is a boolean per turn (slower) -- hence the toggle. Also honoured via env MVB_FUSE_TURNS.
    static void setFuseTurnParts(bool on);
    static bool fuseTurnParts();

    // paintCoating: true  → turn solid drawn at the OUTER (insulation) footprint
    //                       (default; unchanged behaviour for visualisation).
    //               false → turn solid drawn at the CONDUCTING (copper) footprint
    //                       — what FEM winding-loss extraction must mesh. For
    //                       LITZ this is the bare bundle treated as a single
    //                       solid conductor (diameter from MKF).
    static TopoDS_Shape buildTurn(const MAS::Turn& turn,
                                  const MAS::Wire& wire,
                                  const MAS::CoreBobbinProcessedDescription& bobbin,
                                  bool isToroidal,
                                  int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                  int wireRevolutionSegments = DEFAULT_WIRE_REVOLUTION_SEGMENTS,
                                  bool paintCoating = true);

    // Build a turn using ONLY the data on the Turn itself: coordinates,
    // dimensions, cross_sectional_shape, additional_coordinates and
    // rotation. No Wire or Bobbin lookup. Used by drawTurns/drawWinding
    // when the caller does not provide (or want) that surrounding context.
    //
    // Topology heuristic:
    //   - additional_coordinates present → toroidal racetrack turn
    //   - otherwise                       → concentric round-column turn
    //                                       (turn radius = |coordinates|)
    //
    // Throws `std::runtime_error` if any required field is missing.
    // paintCoating must stay true here: a standalone Turn carries only its
    // outer footprint, so requesting the conducting cross-section throws.
    static TopoDS_Shape buildFromTurnAlone(const MAS::Turn& turn,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = true);

    // Resolve the wire cross-section {width, height} without a Turn context (outer
    // footprint or copper per paintCoating). Used by ConductorBuilder, which sweeps whole
    // conductors instead of per-turn solids. Throws on missing wire data (no fallbacks).
    static std::pair<double, double> wireDimensions(const MAS::Wire& wire, bool paintCoating);

    // Fuse all SOLIDs of a compound into one — the same routine the per-turn fuse toggle
    // uses; exposed for ConductorBuilder's per-conductor compounds.
    static TopoDS_Shape fuseSolids(const TopoDS_Shape& compound);

    static void clearCache();
};

} // namespace mvb
