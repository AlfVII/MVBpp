#pragma once

#include "mvb/shapes/ShapeEr.h"

namespace mvb {
namespace shapes {

class ShapeEc : public ShapeEr {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
};

} // namespace shapes
} // namespace mvb
