#pragma once

#include "MAS.hpp"
#include "Utils.h"
#include <TopoDS_Shape.hxx>

namespace mvb {

class BobbinBuilder {
public:
    // polygonSegments: number of segments for round-column cylinders.
    //   >0  → faceted polygonal prism (matches the core/turn tessellation,
    //          so subsequent boolean cuts with the turns are an order of
    //          magnitude faster than NURBS-vs-NURBS).
    //   <=0 → exact NURBS cylinder (legacy behaviour).
    // Has no effect on rectangular-column bobbins (already boxes).
    static TopoDS_Shape buildBobbin(const MAS::CoreBobbinProcessedDescription& bobbin,
                                    double flangeThickness,
                                    bool axisIsY = false,
                                    int polygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS);
};

} // namespace mvb
