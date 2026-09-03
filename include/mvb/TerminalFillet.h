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

#include <functional>
#include <string>
#include <vector>

namespace mvb {

// Fillets every terminal corner in `prims` -- an exact-helix "(terminal stub)" SPIRAL adjacent to a
// "lead seg" SEG, in either order -- in place. Returns how many corners were filleted. Throws,
// naming the corner, when a fillet with radius >= minBend does not fit.
// entranceVariant / exitVariant pick WHICH admissible fillet to use at that corner: the search
// walks the helix-side bend angle from 1 degree up and 0 takes the first that closes, 1 the next,
// and so on (clamped to the last available). The corner alone cannot tell which is best -- the
// arcs leave the helix by a few nanometres to a few microns, and whether that matters depends on
// where the NEIGHBOURING conductor is, which only the fan knows. So the fan tries variants when a
// member does not clear, and hands the winner to the emitter (leadFilletIn / leadFilletOut), the
// same way it hands over the attach leg. Measured on 14_dab: variant 0 leaves 4.16 nm and 2.30 nm
// of interpenetration on two members; other variants clear them.
size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend, double wireRadius,
                             const std::string& who, int entranceVariant = 0,
                             int exitVariant = 0);

// Replaces each straight radial LAYER LINK by a biarc tangent to the two wraps it joins, taking
// the azimuth it needs from each (ABT #969, Alf's option 2). `clears` decides between the
// candidate transitions: the roundest one is the kindest to the wire, but a long transition
// sweeps azimuth and can reach a neighbouring parallel (14_dab: 11.7 um), so the caller -- which
// is the only place that can see the other conductors -- vetoes. Candidates run from the roundest
// down to the tightest that still respects the bend radius; without a predicate the roundest wins.
size_t smoothLayerLinks(std::vector<Primitive>& prims, double minBend, const std::string& who,
                        const std::function<bool(const Primitive&, const Primitive&)>& clears = {});
// Bend radius = kRoundCornerBendFactor * wire radius (the corner rule everywhere else).
inline size_t filletTerminalCorners(std::vector<Primitive>& prims, double minBend,
                                    const std::string& who, int entranceVariant = 0,
                                    int exitVariant = 0) {
    return filletTerminalCorners(prims, minBend, minBend / kRoundCornerBendFactor, who,
                                 entranceVariant, exitVariant);
}

}  // namespace mvb
