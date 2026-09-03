#pragma once

#include "mvb/shapes/ShapeDrum.h"

namespace mvb {
namespace shapes {

// SEMI-SHIELDED DRUM (ABT #362): a wound drum whose winding is overcoated with MAGNETIC
// EPOXY — a polymer binder loaded with magnetic powder, mu ~3-15 — forming a low-permeability
// return shell. The WE-LQS class.
//
// Letters: the drum convention (A..H, see ShapeDrum) plus the shell OUTER ENVELOPE
// J width, K depth, L height. Square LQS-class bodies have J = K. The shell fills from the
// winding and flange rims out to that envelope.
//
// This builder deliberately produces the DRUM ONLY. The shell is not a core piece: MAS
// models it as a coating (coatingType 'magneticEpoxy', whose material resolves against CORE
// materials for its permeability), and it is delivered separately by drawCoreShell so a
// viewer can render it translucent over an opaque drum. Fusing the two would make the
// winding — the thing a semi-shielded part exists to show — permanently invisible.
class ShapeDrumSemishielded : public ShapeDrum {
    // Geometry is the drum's; the shell lives in buildSemishieldedShell().
};

// The magnetic-epoxy shell as its own solid: the J x K x L envelope minus the drum it
// encases. Returns a null shape when the record states no envelope, which is a drum that
// happens to be tagged semi-shielded rather than an error.
TopoDS_Shape buildSemishieldedShell(const std::map<std::string, double>& dims,
                                    const TopoDS_Shape& drum,
                                    int corePolygonSegments);

} // namespace shapes
} // namespace mvb
