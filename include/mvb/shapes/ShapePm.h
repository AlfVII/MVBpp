#pragma once

#include "mvb/shapes/ShapeP.h"
#include <string>

namespace mvb {
namespace shapes {

class ShapePm : public ShapeP {
public:
    explicit ShapePm(const std::string& subtype) : familySubtype_(subtype) {}
protected:
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
private:
    std::string familySubtype_;
};

} // namespace shapes
} // namespace mvb
