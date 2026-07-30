#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

// DRUM (bobbin) core, ABT #331: round post between two flange discs, wound in the groove.
// Dimension convention (MAS drum records): A flange OD, A2 optional second flange OD
// (asymmetric drums, TDG DRC convention), B total height, C post OD, D top flange
// thickness, E groove height, F bottom flange thickness, H optional bore.
class ShapeDrum : public ShapeBuilder {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
};

} // namespace shapes
} // namespace mvb
