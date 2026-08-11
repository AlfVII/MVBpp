#pragma once

#include "MAS.hpp"
#include "Utils.h"
#include "mvb/ConductorBuilder.h"
#include "mvb/NamedShape.h"
#include "mvb/Symmetry.h"
#include <TopoDS_Shape.hxx>
#include <vector>
#include <string>
#include <map>

// Forward declarations for MKF-enriched magnetic overloads
namespace OpenMagnetics { class Magnetic; class Coil; }

namespace mvb {

// Unified configuration for drawMagnetic / drawMagneticToBytes.
// All bindings (Python, WASM) accept this struct so the API surface
// is identical across languages.
struct DrawConfig {
    std::string format                = "step";  // "step" or "stl"
    bool        includeBobbin         = true;
    double      scale                 = 1.0;     // extra USER scale; exporters emit mm natively (ABT #317) — keep 1.0
    int         symmetryPlanes        = 0;       // 0=full, 1=half, 2=quarter
    int         wirePolygonSegments    = DEFAULT_WIRE_POLYGON_SEGMENTS;
    int         corePolygonSegments    = DEFAULT_CORE_POLYGON_SEGMENTS;
    // true  → turns drawn at the OUTER (insulation) diameter (visualisation).
    // false → turns drawn at the CONDUCTING (copper) diameter — required for
    //         FEM winding-loss meshing (LITZ → bare bundle as a solid).
    bool        paintCoating           = true;
    // true → real winding: ONE continuous conductor per (winding, parallel) instead of
    //        independent per-turn loops; MKF enriches with real-winding turn blocking on.
    bool        useRealWindingGeometry = false;
    // Real-winding only. false (OM / drawing) → fast per-run compound conductors for the 3D viewer.
    // true (FEM export) → the slow one-piece/conformal machinery: single body per parallel where a
    // sweep closes, a mitre-jointed conformal compound for dense toroids. Meshable but far slower.
    bool        femReady               = false;
};

class MagneticBuilder {
public:
    // Build geometry and export to STEP ("step") or STL ("stl").
    //
    // Parameters
    // ----------
    //   outputPath            : directory where magnetic.step / magnetic.stl is written
    //   format                : "step" or "stl"
    //   includeBobbin         : include bobbin geometry
    //   scale                 : extra uniform USER scale on top of the exporters'
    //                           native millimetre output (ABT #317). Keep 1.0 for a
    //                           correct mm file; 1000 now double-scales (1e6x too big).
    //   symmetryPlanes        : 0=full, 1=half, 2=quarter domain
    //   wirePolygonSegments   : <=0 = exact torus, >0 = faceted polygon (wire cross-section)
    //   corePolygonSegments   : polygon segments for core cylinders/circles
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const std::string& format = "step",
                             bool includeBobbin = true,
                             double scale = 1.0,
                             int symmetryPlanes = 0,
                             int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                             int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Config-based overload — preferred for bindings.
    std::string drawMagnetic(const MAS::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& cfg) const;

    // Overload that accepts an already-enriched OpenMagnetics::Magnetic
    // to avoid object-slicing issues with MAS::Magnetic
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const std::string& format = "step",
                             bool includeBobbin = true,
                             double scale = 1.0,
                             int symmetryPlanes = 0,
                             int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                             int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Config-based overload — preferred for bindings.
    std::string drawMagnetic(const OpenMagnetics::Magnetic& magnetic,
                             const std::string& outputPath,
                             const DrawConfig& cfg) const;

    // Named-shape overloads. Each returned element carries the logical
    // name (core name / "Turn_<i>" or Turn::get_name / bobbin name) so the
    // identity survives downstream operations (symmetry cut, STEP export
    // with XCAF labels, mesh tagging). These are the only public builders;
    // the prior unnamed variants (buildCore / buildBobbin / buildTurns) have
    // been removed in favour of this Named API.
    std::vector<NamedShape> buildCoreNamed(const MAS::MagneticCore& core,
                                            int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    std::vector<NamedShape> buildTurnsNamed(const MAS::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = true) const;
    std::vector<NamedShape> buildTurnsNamed(const OpenMagnetics::Coil& coil,
                                            const MAS::MagneticCore& core,
                                            int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                            bool paintCoating = true) const;
    NamedShape buildBobbinNamed(const MAS::Coil& coil,
                                const MAS::MagneticCore& core,
                                int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;
    NamedShape buildBobbinNamed(const OpenMagnetics::Coil& coil,
                                const MAS::MagneticCore& core,
                                int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Assemble all geometry as named shapes (no export). Optionally applies
    // symmetry cuts according to symmetryPlanes (0=full, 1=half, 2=quarter).
    // emitCoatingShells: emit BOTH the bare-copper turn ("<turn>") and the outer insulated
    // footprint ("<turn> coating") per turn, so a thermal mesh can resolve the low-k wire enamel
    // (the mesher fragments the overlap into copper core + coating annulus). Implies copper turns.
    // includeInsulation: also emit INSULATION-layer solids ("insulation_layer_<i>") for the
    // inter-layer/inter-section tape, when they carry real thickness (zero-thickness placeholders
    // are skipped). For thermal FEA -- a low-k conduction barrier between windings.
    // useRealWindingGeometry: replace the per-turn closed loops with ONE continuous copper
    // body per (winding, parallel) ("<winding> parallel <k>", ConductorBuilder), enriching
    // through MKF with its real-winding turn blocking on. A magnetic that already carries a
    // geometricalDescription THROWS with the flag on — MKF must (re-)wind; caller-provided
    // turn positions are never silently re-interpreted.
    // Real-winding conductor CENTRELINES (sampled, collision-checked, seam-aimed), for
    // implicit/level-set winding meshing. No solids are built. See
    // ConductorBuilder::PathPolyline.
    std::vector<ConductorBuilder::PathPolyline> buildRealWindingPaths(
        const OpenMagnetics::Magnetic& magnetic) const;

    std::vector<NamedShape> buildAllNamed(const MAS::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
                                          bool paintCoating = true,
                                          bool emitCoatingShells = false,
                                          bool includeInsulation = false,
                                          double coreCoatingThickness = 0.0,
                                          bool useRealWindingGeometry = false,
                                          bool femReady = false) const;
    std::vector<NamedShape> buildAllNamed(const OpenMagnetics::Magnetic& magnetic,
                                          bool includeBobbin = true,
                                          int symmetryPlanes = 0,
                                          int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
                                          int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
                                          bool paintCoating = true,
                                          bool emitCoatingShells = false,
                                          bool includeInsulation = false,
                                          // >0: build the core's insulating coating as a conformal
                                          // shell solid ("<core> coating") of this thickness [m].
                                          double coreCoatingThickness = 0.0,
                                          bool useRealWindingGeometry = false,
                                          bool femReady = false) const;

    // ---- Standalone builders for the unified bindings API -----------------
    //
    // Build a single core piece from a MAS::CoreShape. Validates the shape
    // through OpenMagnetics::CorePiece::factory (which also fills in derived
    // parameters), then dispatches to mvb's ShapeBuilder for geometry.
    // Returns a single NamedShape with name = shape.name (or family code).
    NamedShape buildCorePieceNamed(const MAS::CoreShape& shape,
                                    int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Build a bobbin from a fully-populated MAS::Bobbin. Throws if
    // `processed_description` is not present. `axisIsY=true` orients the
    // bobbin tube along Y (matches concentric-core convention); pass false
    // for toroidal bobbins (kept along Z).
    NamedShape buildBobbinNamedFromBobbin(const MAS::Bobbin& bobbin,
                                          bool axisIsY = true,
                                          int polygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS) const;

    // Build a list of turns where every Turn carries its own dimensions and
    // cross_sectional_shape (no Wire/bobbin lookup). Throws if any turn is
    // missing required fields. Toroidal layout is auto-detected from the
    // presence of `additional_coordinates`; otherwise concentric round
    // column is assumed.
    // paintCoating must stay true: standalone turns carry only their outer
    // footprint, so the conducting cross-section cannot be recovered here and
    // paintCoating=false throws (use a full Magnetic JSON via drawMagnetic).
    std::vector<NamedShape> buildTurnsNamedFromTurns(
        const std::vector<MAS::Turn>& turns,
        int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
        bool paintCoating = true) const;

    // Turns ONLY, in real-winding form: one continuous copper body per (winding, parallel),
    // with the real leads, pitch and dragbacks. The turns-only counterpart of
    // buildAllNamed(..., useRealWindingGeometry=true), for a consumer that draws core,
    // bobbin and turns as SEPARATE meshes (the web 3D viewer, which colours and toggles
    // them independently) and so cannot use the single-call assembly.
    //
    // The magnetic must already be enriched through MKF's real-winding autocomplete —
    // same precondition as buildAllNamed with the flag on. The core is built internally
    // (not returned) because the conductor builder aims the terminal leads at the real
    // window opening; see buildRealWindingConductorsNamed.
    //
    // diagnosticSkipCollisionCheck builds the conductors even when two of them overlap, so
    // the overlap can be LOOKED AT instead of only read about in an exception. The result is
    // interpenetrating copper: a picture, never a part. Tools and tests investigating a
    // specific refusal only — nothing in the library sets it.
    std::vector<NamedShape> buildRealWindingTurnsNamed(
        const OpenMagnetics::Magnetic& magnetic,
        int wirePolygonSegments = DEFAULT_WIRE_POLYGON_SEGMENTS,
        int corePolygonSegments = DEFAULT_CORE_POLYGON_SEGMENTS,
        bool paintCoating = true,
        bool femReady = false,
        bool diagnosticSkipCollisionCheck = false) const;

  private:
    // Single implementation of real-winding conductor emission, shared by buildAllNamed
    // and buildRealWindingTurnsNamed so the assembly and the viewer cannot drift apart.
    // coreShapes are handed to ConductorBuilder as lead-aiming obstacles (ignored for
    // toroids); they are not part of the result.
    std::vector<NamedShape> buildRealWindingConductorsNamed(
        const OpenMagnetics::Magnetic& magnetic,
        const std::vector<NamedShape>& coreShapes,
        int wirePolygonSegments,
        bool paintCoating,
        bool emitCoatingShells,
        bool femReady,
        bool diagnosticSkipCollisionCheck = false) const;
};

} // namespace mvb
