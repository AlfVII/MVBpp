#pragma once

#include "mvb/shapes/ShapeBuilder.h"

namespace mvb {
namespace shapes {

// MOLDED core (ABT #357): a coil compression-moulded inside a homogeneous soft-magnetic
// composite block — the WE-MAPI / IHLP / XAL class. There is no air gap and no separate
// bobbin: the composite IS the body, and the winding sits in a cavity inside it.
//
// Dimension convention, matching MKF's CorePieceMolded (which documents it as the pot-core
// P/PM convention): A body width, B body height along the coil axis, C body depth,
// D coil-cavity height, E coil-cavity outer diameter, F coil-cavity inner diameter — the
// composite post under the coil bore.
//
// NOTE the letters D/E/F are the ones MKF implements and the ones the WE ingest emits. The
// MAS commit that introduced the family (ABT #357) described D as the cavity inner diameter
// and F as its height, which is NOT what shipped; MKF's own header comment and its
// validation ("cavity outer diameter E must exceed the cavity inner diameter F") are
// authoritative here.
//
// The body is drawn as a closed block: the cavity is INTERNAL, so from outside a moulded
// inductor is a plain rectangular body. It is rendered translucent so the winding inside is
// visible, which is the only way to see anything of the construction.
class ShapeMolded : public ShapeBuilder {
protected:
    TopoDS_Face buildProfile(const std::map<std::string, double>& dims) const override;
    TopoDS_Shape buildWindingWindow(const std::map<std::string, double>& dims) const override;
};

} // namespace shapes
} // namespace mvb
