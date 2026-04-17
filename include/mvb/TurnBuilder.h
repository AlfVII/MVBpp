#pragma once

#include "MAS.hpp"
#include <TopoDS_Shape.hxx>
#include <vector>

namespace mvb {

class TurnBuilder {
public:
    static TopoDS_Shape buildTurn(const MAS::Turn& turn,
                                  const MAS::Wire& wire,
                                  const MAS::CoreBobbinProcessedDescription& bobbin,
                                  bool isToroidal);
};

} // namespace mvb
