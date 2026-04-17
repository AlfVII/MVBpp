#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

class ShapeP : public ShapeBuilder {
public:
    TopoDS_Shape applyMachining(const TopoDS_Shape& piece,
                                const MAS::Machining& machining,
                                const std::map<std::string, double>& dims) const override;
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
};

} // namespace shapes
} // namespace mvb
