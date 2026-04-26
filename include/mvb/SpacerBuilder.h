#pragma once

#include "MAS.hpp"
#include <TopoDS_Shape.hxx>
#include <vector>

namespace mvb {

class SpacerBuilder {
public:
    // Build spacer boxes from MAS geometrical description entries with
    // type == SPACER. Each spacer is a box of (dimensions[0], dimensions[2],
    // dimensions[1]) placed at coordinates (X, Y, Z) in MVB++'s column-along-Y
    // frame, centred on its Z (height) coordinate.
    //
    // Returns an empty vector if no spacers are found.
    static std::vector<TopoDS_Shape> buildSpacers(
        const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription);

    // Convenience: fuse all spacer shapes into a single compound (or the
    // single spacer itself if only one exists). Returns a null shape when
    // the list is empty.
    static TopoDS_Shape buildSpacersCompound(
        const std::vector<MAS::CoreGeometricalDescriptionElement>& geometricalDescription);
};

} // namespace mvb
