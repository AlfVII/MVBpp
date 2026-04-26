#pragma once

#include "MAS.hpp"
#include "mvb/Utils.h"
#include <TopoDS_Shape.hxx>
#include <vector>

namespace mvb {

class TurnBuilder {
public:
    static TopoDS_Shape buildTurn(const MAS::Turn& turn,
                                  const MAS::Wire& wire,
                                  const MAS::CoreBobbinProcessedDescription& bobbin,
                                  bool isToroidal,
                                  int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                  int wireRevolutionSegments = DEFAULT_WIRE_REVOLUTION_SEGMENTS);
    static void clearCache();
};

} // namespace mvb
