# Handoff Summary: MVB++ Bobbin & Turn Geometry Fixes

## Project Context
**MVB++** is a C++ magnetics geometry builder using OpenCASCADE (OCCT 7.9), located at `/home/alf/OpenMagnetics/MVB++/`. It must match the behavior of:
- **Python MVB** (`/home/alf/OpenMagnetics/MVB/`) - CadQuery-based
- **MVB.js** (`/home/alf/OpenMagnetics/MVB.js/`) - Replicad-based

## Goal
Fix MVB++ STEP output to match Python MVB for bobbins and turns, specifically for examples `01_simple_inductor_etd34_n87` and `03_buck_inductor_pq3230_n95`.

---

## Completed Work

### 1. Fixed `column_width` / `column_depth` interpretation
**Discovery:** MKF sets `column_width` in bobbin `processedDescription` to **outer radius/half-width** (`core_column_width/2 + wall_thickness`), not diameter. MVB++ was incorrectly dividing by 2.0.

**Files:**
- `src/MagneticBuilder.cpp:54-77` - `patchBobbinDimensions` now uses `centralCol.get_width() / 2.0 + wallThickness`
- `src/TurnBuilder.cpp:438-439` - removed `/2.0` from `half_col_width` and `half_col_depth`

### 2. Fixed rectangular turns on round columns
**Problem:** `BRepOffsetAPI_MakePipe` with a rectangular profile on a polygon-circle path produced degenerate/empty solids in OCCT 7.9.

**Fix:** Replaced with `BRepPrimAPI_MakeRevol` of a rectangle face (360° revolve around Y axis), matching MVB.js approach.

**File:** `src/TurnBuilder.cpp:117-158` (`build_concentric_round_column_turn`)

### 3. Fixed round wire turns on rectangular columns
**Problem:** After removing `/2.0`, `BRepOffsetAPI_MakePipe` for round wire on rectangular columns threw `BRepAdaptor_Curve::No geometry` in OCCT 7.9.

**Fix:** Implemented manual construction using 4 cylinders (straight segments) + 4 quarter tori (corners), then fused together.

**File:** `src/TurnBuilder.cpp:235-289` (`build_concentric_rect_column_turn`)

### 4. Added bobbin cutting to match Python MVB
**Problem:** Python MVB cuts the bobbin with cores and turns before export. MVB++ did not.

**Fix:** Added `cut_bobbin()` utility in `src/Utils.cpp:210-224` and called it in both `MagneticBuilder::drawMagnetic` overloads.

**Files:** `src/Utils.cpp`, `src/Utils.h`, `src/MagneticBuilder.cpp:127-145` and `:184-199`

### 5. Fixed `build_circle_profile` transformation bug
**Problem:** The profile builder was placing profiles in the wrong plane due to incorrect `gp_Ax2` construction.

**Fix:** Corrected axis orientation in `build_circle_profile`.

**File:** `src/TurnBuilder.cpp:66-93`

### 6. Reconstructed accidentally truncated `MagneticBuilder.cpp`
**Problem:** An `edit` command truncated `src/MagneticBuilder.cpp` from ~375 lines to 161 lines. Reconstructed from memory of earlier reads.

**File:** `src/MagneticBuilder.cpp` - now complete with `buildCore`, `buildBobbin`, `buildTurns`, `patchBobbinDimensions`, `extractBobbinDims`

### 7. Fixed critical ETD/ER core crash
**Problem:** `etd49_5t` test and ETD34 example crashed with SIGSEGV in `ShapeEr::applyMachining`.

**Root Cause:** A **dangling reference bug** in `MagneticBuilder::buildCore`:
```cpp
// BUGGY - temporary optional destroyed too early
for (const auto& mach : *piece.get_machining()) { ... }
```
`piece.get_machining()` returns `std::optional<std::vector<Machining>>` **by value**. In this GCC 14 / C++23 build, the temporary optional's lifetime was NOT properly extended through the range-for, causing `mach` to reference garbage memory. `machining.get_coordinates()` returned a corrupted vector with size `809564666238718298`, and accessing `coords[0]` caused the segfault.

**Fix:** Store the optional in a local variable first:
```cpp
auto machiningOpt = piece.get_machining();
if (machiningOpt) {
    for (const auto& mach : *machiningOpt) { ... }
}
```

**Files:** `src/MagneticBuilder.cpp:258-264`

Also cleaned up `ShapeEr::applyMachining` to use `BRepPrimAPI_MakeCylinder` with correct axis directly instead of creating a Z-cylinder then rotating it.

**File:** `src/shapes/ShapeEr.cpp`

Also fixed `ShapeEr::buildWindingWindow` to use exact cylinders (`segments=0`) instead of polygon-approximated ones to avoid OCCT coincident-face boolean issues.

---

## Current Test Status
All tests pass:
```bash
cd /home/alf/OpenMagnetics/MVB++ && cmake --build build -j$(nproc)
LD_LIBRARY_PATH=build/occt-install/lib ./build/mvb_tests
```
**Result:** `All tests passed (82 assertions in 8 test cases)`

The step generator no longer crashes on ETD34:
```bash
LD_LIBRARY_PATH=build/occt-install/lib ./build/mvbpp_step_generator \
  -o /tmp/mvb_examples_steps/01_simple_inductor_etd34_n87.step \
  /home/alf/OpenMagnetics/MAS/examples/01_simple_inductor_etd34_n87.json
```

---

## Remaining Issues

### 1. ETD34 bounding box mismatch
**Comparison result:**
```
MVB++ vol=10664.1 solids=26 | MVB vol=10250.4 solids=26
BBOX mismatch: MVB++ dims=[16.14, 34.2, 34.6] vs MVB dims=[23.0, 34.2, 34.6]
```

The 34.2 and 34.6 (core width/height) match. The 16.14 vs 23.0 discrepancy is the **third sorted dimension**.

**Analysis:**
- MVB++'s 16.14 mm corresponds to the **turn ring diameter** in the XZ plane: `2 * (radial_pos + wire_radius) = 2 * (7.75 + 0.4275) ≈ 16.36 mm`
- Python MVB reports 23.0 mm, which suggests Python's limiting dimension is **not** the turns but the **bobbin flange width after cutting**

This implies Python MVB's turns might be built differently (smaller effective ring), OR Python's bobbin flange is being cut differently by OCCT vs CadQuery's boolean operations, OR the `flangeExtension` default in Python is different.

**Key question:** Does Python MVB actually use `flangeExtension = 0.002` (the MVB++ default when missing), or does it extract it from somewhere else? For ETD34 "Basic" bobbin, there is no `functionalDescription` with dimensions.

**Action needed:** Trace exactly how Python MVB computes the ETD34 bobbin flange width and why the post-cut width is ~23.0 mm instead of ~16.14 mm.

### 2. PQ3230 rectangular wire behavior
**Problem:** The PQ3230 example uses rectangular wire but the wire data lacks `crossSectionalShape`.

**Python MVB bug:** Python defaults to round wire with `outer_diameter = dims[0]` when `crossSectionalShape` is missing. This is a known Python bug.

**MVB++ behavior:** MVB++ correctly detects rectangular wire from `wire.type == RECTANGULAR` and builds rectangular turns. This produces different geometry (rectangular wire rings vs round wire tori).

**Decision needed:** Should MVB++ replicate Python's buggy behavior to pass comparison, or keep the correct behavior and accept the mismatch?

---

## Architecture Notes

### Key Files
| File | Purpose |
|------|---------|
| `src/MagneticBuilder.cpp` | Orchestrates core/bobbin/turn building, applies cutting, calls exporter |
| `src/BobbinBuilder.cpp` | Builds bobbin geometry (cylinders/boxes for body + flanges) |
| `src/TurnBuilder.cpp` | Builds individual turn solids (tori, revolved rings, manual cylinder+tori) |
| `src/Utils.cpp` | Utilities: `cut_bobbin`, `build_polygon_cylinder`, `magnetic_autocomplete_safe` |
| `src/shapes/ShapeE.cpp` | E-family core profile, winding window, machining |
| `src/shapes/ShapeEr.cpp` | ER/ETD round-column winding window and machining |
| `src/shapes/ShapeBuilder.cpp` | Base class: `buildPiece`, `applyMachining`, `extrude` |
| `src/shapes/Factory.cpp` | Maps `CoreShapeFamily` to builder subclass |
| `tests/test_step_comparison.cpp` | Tests comparing against Python MVB reference STEPs |

### Coordinate Systems
- **Cores:** Built in XY-extruded-Z, then translated by `-B`, then rotated `-90°` around X (so Z→-Y, Y→Z). Final orientation: Y is the column axis.
- **Bobbins:** `BobbinBuilder::buildBobbin` takes `axisIsY`. For non-toroidal cores, it rotates the bobbin `-90°` around X to align with the core.
- **Turns:** Built with Y as the revolution/sweep axis for concentric cores.

### MKF Data Flow
1. Raw MAS JSON → `patch_dimension_nominals` (adds missing `nominal` to min/max dimensions)
2. → `MAS::Magnetic` parsing
3. → `magnetic_autocomplete_safe(magnetic)` (enriches core geo desc + turns desc via MKF)
4. → `MagneticBuilder::drawMagnetic(enriched, ...)`

### Comparison Method
Use `/tmp/compare_example.py` (copied from `tests/debug_bobbin.py`) to generate Python MVB reference and compare volumes, solid counts, and bounding boxes.

---

## Recommended Next Steps
1. **Investigate ETD34 bobbin flange width:** Check Python MVB's exact bobbin dimensions after `cut()` with E-core. Compare `CadQueryBuilder.StandardBobbin.get_bobbin_flanges` behavior for "Basic" bobbins. The 23.0 mm likely comes from a different `flangeExtension` or from how CadQuery handles the boolean cut.
2. **Decide on PQ3230:** Either match Python's round-wire bug or update the test to expect correct rectangular-wire geometry.
3. **Regenerate reference STEPs** once the new behavior is accepted as correct.
4. **Check for similar dangling-reference bugs** elsewhere in the codebase (search for `*.*get_machining()` in range-fors or similar patterns with returned-by-value optionals).
