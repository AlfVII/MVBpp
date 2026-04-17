#pragma once

#include "mvb/shapes/ShapeE.h"

namespace mvb {
namespace shapes {

// Round-center-column E-family cores (ER, ETD)
class ShapeEr : public ShapeE {
protected:
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                const MAS::Machining& machining,
                                const std::map<std::string, double>& dims) const override;
};

} // namespace shapes
} // namespace mvb
