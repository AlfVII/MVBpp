#pragma once

#include "MAS.hpp"
#include <TopoDS_Shape.hxx>

namespace mvb {

class BobbinBuilder {
public:
    static TopoDS_Shape buildBobbin(const MAS::CoreBobbinProcessedDescription& bobbin, double flangeThickness, bool axisIsY = false);
};

} // namespace mvb
