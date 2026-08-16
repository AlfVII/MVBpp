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
    // ABT #685: names for the INDIVIDUAL SOLIDS of a multi-solid shape, in the order
    // TopExp_Explorer visits them. A round-wire conductor is a compound of per-primitive solids
    // (the conformal mitre assembly), and with only the compound named a viewer falls back to
    // numbering them itself — "Primary parallel 001", "…012" — a convention that lives in the
    // viewer, not in the file, so a solid cannot be named unambiguously in conversation. Filling
    // this makes exportSTEP write the compound as an assembly with a real name on every part.
    // Empty (the default) keeps the previous single-product behaviour exactly.
    std::vector<std::string> partNames;

    NamedShape() = default;
    NamedShape(TopoDS_Shape s, std::string n)
        : shape(std::move(s)), name(std::move(n)) {}
    NamedShape(TopoDS_Shape s, std::string n, std::vector<std::string> parts)
        : shape(std::move(s)), name(std::move(n)), partNames(std::move(parts)) {}
};

} // namespace mvb
