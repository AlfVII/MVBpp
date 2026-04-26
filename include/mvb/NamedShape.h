#pragma once

#include <TopoDS_Shape.hxx>
#include <string>
#include <utility>
#include <vector>

namespace mvb {

// A geometric shape paired with a human-readable identifier. Names are
// produced at build time from MAS metadata (`MagneticCore::get_name`,
// `Turn::get_name`, bobbin name, etc.) with fallbacks, then carried
// through every downstream operation — boolean cuts, STEP export, mesh
// tagging — so the original logical identity survives the pipeline.
struct NamedShape {
    TopoDS_Shape shape;
    std::string  name;

    NamedShape() = default;
    NamedShape(TopoDS_Shape s, std::string n)
        : shape(std::move(s)), name(std::move(n)) {}
};

} // namespace mvb
