#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

class ShapeUt : public ShapeBuilder {
public:
    TopoDS_Shape buildPiece(const MAS::CoreShape& shapeData) const override;
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>&) const override {
        return TopoDS_Face();  // unused — buildPiece is fully overridden
    }
};

} // namespace shapes
} // namespace mvb
