#pragma once
// THE TERMINAL CORNER IS A FILLET, NOT A MITRE (ABT #961 / #839 follow-up, 2026-09-02).
//
// The terminal stub used to be the helix's OSCULATING CIRCLE so that its corner with the lead
// could be mitred analytically (a spiral pipe end cannot be cut reliably by OCC). Every circle
// leaves the helix by the torsion cubic (-kappa*tau/6 * s^3), ~100 nm over a mm-scale stub, and on
// boost_inductor_complete that ate a sibling lead's +60 nm margin (47 nm certified overlap).
// Now the stub IS the helix (an exact SPIRAL primitive), and the corner between it and the lead
// is what a bent wire actually is: a tangent BIARC (two circular arcs, radii >= the bend radius),
// consuming equal tangent lengths from the helix's end and the lead's start. Every joint is
// tangent, so nothing is mitred and nothing is bridged.
#include <mvb/WireAssembler.h>

#include <string>
#include <vector>

namespace mvb {

// Fillets every terminal corner in `prims` -- an exact-helix "(terminal stub)" SPIRAL adjacent to a
// "lead seg" SEG, in either order -- in place. Returns how many corners were filleted. Throws,
// naming the corner, when a fillet with radius >= minBend does not fit.
size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                             const std::string& who);
// Bend radius = kRoundCornerBendFactor * wire radius (the corner rule everywhere else).
inline size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend,
                                    const std::string& who) {
    return filletTerminalCorners(prims, minBend, minBend / kRoundCornerBendFactor, who);
}

}  // namespace mvb
