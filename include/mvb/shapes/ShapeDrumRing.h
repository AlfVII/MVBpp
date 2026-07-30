#pragma once

#include "mvb/shapes/ShapeDrum.h"

namespace mvb {
namespace shapes {

// DRUM RING (shielded drum) core, ABT #366: a drum closed by a concentric shield ring, sold
// as a matched pair and described by ONE shape record. Drum letters as ShapeDrum (A flange
// OD, A2 optional second flange OD, B total height, C post OD, D top flange thickness,
// E groove height, F bottom flange thickness, H optional bore) plus the ring: J ring OD,
// K ring ID, L ring height.
//
// The ring does NOT touch the drum: K > A by the annular clearance that closes the magnetic
// circuit through air (MKF synthesizes those clearances as structural residual gaps). The
// assembly is therefore emitted as the drum solid plus a disjoint ring solid.
class ShapeDrumRing : public ShapeDrum {
protected:
    TopoDS_Shape applyExtras(const std::map<std::string, double>& dims,
                             const TopoDS_Shape& piece) const override;
};

} // namespace shapes
} // namespace mvb
