#pragma once

#include "mvb/shapes/ShapeEp.h"

namespace mvb {
namespace shapes {

class ShapeEpx : public ShapeEp {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
};

} // namespace shapes
} // namespace mvb
