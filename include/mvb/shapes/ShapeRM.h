#pragma once

#include "mvb/shapes/ShapeP.h"

namespace mvb {
namespace shapes {

class ShapeRM : public ShapeP {
public:
    TopoDS_Shape buildPiece(const MAS::CoreShape& shapeData) const override;
};

} // namespace shapes
} // namespace mvb
