#pragma once

#include "MAS.hpp"
#include <string>

namespace OpenMagnetics { class Magnetic; }

namespace mvb {

// Renders a dimensioned technical-drawing-style SVG of the given core's
// front-elevation cross-section (section plane XY in MVB++'s Y-axis-column
// convention; this is the view Python MVB calls "FrontView"). Annotates:
//   - Core-shape B (height) and D (window height).
//   - Every gap in functionalDescription.gapping as a vertical arrow
//     labeled with the gap length in mm or μm.
// Requires the Magnetic to have been enriched by magnetic_autocomplete so
// geometricalDescription, gapping coordinates and processedDescription
// columns are populated.
class SectionDrawing {
public:
    // FrontView — X horizontal, Y vertical (column axis). Annotates B, D
    // (core shape) plus every gap length and the chunk sizes between gaps.
    // For toroidal cores, C is drawn instead of B.
    static std::string drawDimensionedFrontView(
        const OpenMagnetics::Magnetic& magnetic,
        double width_px = 800.0,
        double label_font_px = 12.0,
        const std::string& projection_color = "#000000",
        const std::string& dimension_color = "#1976d2");

    // TopView — X horizontal, Z vertical (depth). Annotates the full per-
    // family set of dimensions (A, E, F, C, G, H, J, L, K, F2 as applicable).
    // Dispatches on the core-shape family: E-family (default), UR, UT, T.
    static std::string drawDimensionedTopView(
        const OpenMagnetics::Magnetic& magnetic,
        double width_px = 800.0,
        double label_font_px = 12.0,
        const std::string& projection_color = "#000000",
        const std::string& dimension_color = "#1976d2");

    // FrontView showing only gap dimensions (no core shape labels).
    // Used by the gapping customizer panel.
    static std::string drawCoreGappingTechnicalDrawing(
        const OpenMagnetics::Magnetic& magnetic,
        double width_px = 800.0,
        double label_font_px = 12.0,
        const std::string& projection_color = "#000000",
        const std::string& dimension_color = "#1976d2");

    // Convenience: build the view and write to disk.
    static void writeDimensionedFrontView(
        const OpenMagnetics::Magnetic& magnetic,
        const std::string& outputPath);
    static void writeDimensionedTopView(
        const OpenMagnetics::Magnetic& magnetic,
        const std::string& outputPath);

    // Unified pictorial-view dispatcher used by drawView in the bindings.
    //   - plane="XY" : front elevation (delegates to drawDimensionedFrontView
    //                  if dimensions=true, otherwise plain section SVG).
    //   - plane="XZ" : top view (delegates to drawDimensionedTopView when
    //                  dimensions=true).
    //   - plane="YZ" : not currently dimensioned — falls back to a plain
    //                  section SVG cut at x=offset.
    // `offset` (metres) shifts the cut plane along its normal. For
    // dimensioned views only offset=0.0 is supported (throws otherwise).
    static std::string drawView(
        const OpenMagnetics::Magnetic& magnetic,
        const std::string& plane = "XY",
        double offset = 0.0,
        bool dimensions = true,
        double width_px = 800.0,
        double label_font_px = 12.0,
        const std::string& projection_color = "#000000",
        const std::string& dimension_color = "#1976d2");
};

} // namespace mvb
