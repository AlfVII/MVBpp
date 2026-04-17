# MVB++ Bobbin Progress Report

## Goal
Make MVB++ (C++/OCCT) produce bobbins that exactly match MVB (Python/CadQuery) output. This is a prerequisite step before moving on to full assembly comparisons (cores + bobbins + turns).

## What was discovered

### 1. MVB++ was building "dummy" bobbins
- `BobbinBuilder` originally produced solid blocks/cylinders covering the winding window area.
- Python MVB `StandardBobbin` builds hollow bobbins with walls, central holes, top/bottom flanges, and flange holes.
- These geometries have completely different volumes and bounding boxes, so they can never match.

### 2. Cylinder conventions differ between CadQuery and OCCT
- **CadQuery** `cq.Workplane("XY").cylinder(height, radius)` creates a cylinder **centered at the origin** (z range = [-height/2, height/2]).
- **OCCT** `BRepPrimAPI_MakeCylinder(radius, height)` creates a cylinder with its **base at z=0** (z range = [0, height]).
- This caused round bobbin hole cylinders to be positioned differently in MVB++ vs Python, changing the boolean-cut overlap with the flanges and producing different volumes.

### 3. Python custom comparison script had a bug
- `scripts/generate_mvb_bobbins.py` originally translated the round cylinder body down by `height/2`, making it misaligned with the flanges.
- Fixed by removing the extra translate so the cylinder stays centered, matching the real MVB `StandardBobbin` behavior.

### 4. `wallThickness` handling differed
- Python `StandardBobbin` honors `wallThickness = 0` from the JSON (produces a body with zero wall thickness, essentially just flanges + thin rings).
- MVB++ was forcing `wallThickness = 0.0005` whenever it was `<= 0`, which added a body tube that the Python reference did not have.
- Fix: only default `wallThickness` when it is `< 0` or `NaN`, preserving explicit zeros.

### 5. MKF enrichment vs Python defaults
- Raw MAS test JSONs do not contain `flangeThickness` or `flangeExtension` in the bobbin `functionalDescription`.
- Python `StandardBobbin` uses hardcoded defaults (`0.001`, `0.002`) when these are missing.
- MKF `magnetic_autocomplete_safe` enriches the bobbin and sometimes calculates very different flange dimensions (e.g., ETD49 gets `flangeExtension ≈ 0.014`).
- This causes assembly volume/bbox mismatches between MVB++ (using MKF-enriched dims) and the Python reference (using `StandardBobbin` defaults).

## Accomplished

### Library changes
- **Rewrote `BobbinBuilder`** (`include/mvb/BobbinBuilder.h`, `src/BobbinBuilder.cpp`)
  - Added `buildBobbin(bobbin, flangeThickness, flangeExtension, axisIsY)` overload.
  - Builds proper hollow rectangular bobbins with walls, central hole, top/bottom flanges, and flange holes.
  - Builds proper hollow round bobbins with inner/outer cylinders, top/bottom flanges with cylindrical holes.
  - Aligns to Y for concentric cores via a `-90°` rotation around X (matching Python's core orientation).
  - Uses strict validation: throws `std::invalid_argument` for missing/invalid flange dimensions.

- **Updated `MagneticBuilder`** (`src/MagneticBuilder.cpp`)
  - Added `getDimValue()` helper to extract doubles from `MAS::Dimension` variants.
  - Added `extractBobbinDims()` to read `flangeThickness` and `flangeExtension` from `MAS::Bobbin::functional_description::dimensions`.
  - `buildBobbin()` now passes extracted flange dimensions to `BobbinBuilder` instead of using dummy defaults.
  - Throws exceptions if required bobbin dimensions are missing (no silent fallbacks).

### Test / reference changes
- **Updated `tests/generate_reference_steps.py`**
  - Replaced `_build_dummy_bobbins()` with `_build_proper_bobbin()` that uses Python MVB's `StandardBobbin` class.
  - Rotates the bobbin to Y for concentric cores so the reference matches MVB++ orientation.
  - Regenerated all reference STEPs in `tests/reference_steps/`.

- **Updated `tests/test_step_comparison.cpp`**
  - Added `loadEnrichedMagneticFromJson()` helper that calls `mvb::magnetic_autocomplete_safe()` and then patches bobbin flange dimensions to `0.001` / `0.002` when missing, ensuring test data matches the Python reference defaults.
  - Added `compareFullAssemblyAgainstReference` overload for `const OpenMagnetics::Magnetic&` to avoid object-slicing issues.
  - Assembly tests now use enriched magnetic data so `BobbinBuilder` receives valid flange dimensions.

### Standalone comparison tool
- Created `tools/generate_mvbpp_bobbins.cpp` to generate rectangular and round bobbins in MVB++ using the same dimensions as the Python script.
- Verified volume and bounding-box equality:
  - **Rectangular**: 319.498 mm³, bbox `[-4.5,4.5] × [-4.5,4.5] × [-6,6]` mm
  - **Round**: 640.575 mm³, bbox `[-10,10] × [-5,5] × [-6,6]` mm

## Current status

- `BobbinBuilder` produces proper bobbins that match Python MVB `StandardBobbin` geometry.
- `rect_one_turn` assembly test **passes** after fixing `wallThickness` and flange dimension alignment.
- `etd49_5t` assembly test is blocked by the MKF-enrichment flange-dimension mismatch (being fixed by post-enrichment patching in the test loader).
- There is a **compilation error** in `tests/test_step_comparison.cpp`: `OpenMagnetics::Coil::get_mutable_bobbin()` returns `std::variant<OpenMagnetics::Bobbin, std::string>`, not `std::variant<MAS::Bobbin, std::string>`. The test code uses `MAS::Bobbin` in the variant access, which fails to compile.

## Next steps
1. Fix the compilation error in `tests/test_step_comparison.cpp` by using `OpenMagnetics::Bobbin` instead of `MAS::Bobbin` when accessing the variant from `OpenMagnetics::Coil`.
2. Rebuild and run all tests to confirm both assembly comparisons pass.
3. Clean up any temporary debug prints in `src/Utils.cpp`.
4. Report success on bobbin matching and ask whether to proceed to turn comparison or full assembly batch generation.

## Relevant files
- `include/mvb/BobbinBuilder.h`
- `src/BobbinBuilder.cpp`
- `src/MagneticBuilder.cpp`
- `include/mvb/Utils.h`
- `src/Utils.cpp`
- `tests/test_step_comparison.cpp`
- `tests/generate_reference_steps.py`
- `tools/generate_mvbpp_bobbins.cpp`
- `scripts/generate_mvb_bobbins.py`
