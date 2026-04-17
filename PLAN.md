# MVB++ Overall Project Plan

## Vision
Build **MVB++**, a C++ equivalent of the Python/CadQuery-based `../MVB` library, using OpenCASCADE 7.9.3 directly to generate 3D geometry (cores, bobbins, turns) from MAS JSON. The output must match Python MVB reference STEPs, including full assemblies with bobbin and turns.

## High-level goals
1. **Self-contained build** — fetch and build OCCT 7.9.3 from source automatically.
2. **Direct MAS.hpp usage** — use generated `MAS` C++ types everywhere, no lightweight wrappers.
3. **Meters everywhere** — no mm scaling internally; scale only at export time (Python MVB exports in mm).
4. **MKF integration** — safely use `OpenMagnetics::magnetic_autocomplete` to enrich raw MAS JSONs (adds `geometricalDescription`, `turnsDescription`, bobbin details) without object slicing or Coil constructor crashes.
5. **Full assembly export** — generate STEP files containing core pieces, bobbin, and turns that match Python MVB reference outputs.

## Coordinate conventions
- **Concentric cores** (E, T, U, etc.): profile built in XY plane, extruded along Z, then rotated `-90°` around X so the column axis aligns with **Y**.
- **Toroidal cores** (T family): profile built in XY plane, extruded along Z, origin at toroid center.

## Architecture
```
include/mvb/         — public headers (Utils, MagneticBuilder, StepExporter, BobbinBuilder, TurnBuilder)
src/shapes/          — one builder per core shape family (ShapeE, ShapeT, ShapeEr, etc.)
src/                 — library implementation (BobbinBuilder, TurnBuilder, StepExporter, etc.)
tests/               — Catch2 tests + Python reference generator + OCCT-based comparison tests
tools/               — CLI utilities (step generator, bobbin generator)
```

## Phase 1: Infrastructure (done)
- [x] Set up CMake with OCCT 7.9.3 external project
- [x] Integrate MAS.hpp auto-generated from JSON schema
- [x] Integrate nlohmann/json and Catch2
- [x] Build system produces static `libmvb++.a` and test executable

## Phase 2: Core shapes (done)
- [x] Implement shape builders: `ShapeE`, `ShapeT`, `ShapeEr`, `ShapeP`, `ShapeEtd`, `ShapeU`, `ShapeToroidal`
- [x] Generic profile → extrude → cut winding window → apply machining → intrinsic rotation pipeline
- [x] Tests verify individual core STEPs against Python MVB references (E core, T core)

## Phase 3: Bobbin matching (in progress)
- [x] Analyze Python MVB `StandardBobbin` geometry (hollow body + flanges + holes)
- [x] Implement proper `BobbinBuilder` with rectangular and round support
- [x] Fix OCCT/CadQuery cylinder convention mismatch (base-at-zero vs centered)
- [x] Fix `wallThickness = 0` handling (Python keeps it zero, MVB++ was forcing 0.0005)
- [x] Update Python reference generator to use real `StandardBobbin` instead of dummy solids
- [x] Standalone bobbin comparison tool verifies rectangular and round volumes/bboxes match exactly
- [x] Wire flange dimensions through `MagneticBuilder`
- [ ] Fix remaining compilation error in `test_step_comparison.cpp` (`MAS::Bobbin` vs `OpenMagnetics::Bobbin` variant)
- [ ] Make both `rect_one_turn` and `etd49_5t` assembly tests pass

## Phase 4: Turns (blocked by bobbin phase)
- [ ] Build turns from `turnsDescription` with proper wire profiles (round, rectangular, foil)
- [ ] Fix or work around `BRepOffsetAPI_MakePipe` crashes/segfaults on certain turn geometries
- [ ] Match Python MVB turn positioning, orientation, and cross-sections
- [ ] Verify per-turn volumes and full assembly volumes match references

## Phase 5: Batch generation & validation (blocked by turns)
- [ ] Run batch STEP generation over all `../MAS/examples` JSONs
- [ ] Compare every generated STEP against Python MVB reference (volume, bbox, solid count)
- [ ] Investigate and fix remaining mismatches
- [ ] Document which examples pass/fail and why

## Phase 6: Cleanup & documentation
- [ ] Remove all temporary debug prints
- [ ] Run lint/typecheck if available
- [ ] Update README with build instructions, design notes, and known limitations
- [ ] Commit final state (only when explicitly requested)

## Known blockers / risks
1. **`BRepOffsetAPI_MakePipe` crashes** — OCCT throws `StdFail_NotDone` or segfaults for ~12/27 MAS examples during turn generation. May need to replace `MakePipe` with manual primitive construction (sweeping a profile along a polyline by fusing translated copies).
2. **MKF bobbin dimension divergence** — `magnetic_autocomplete` produces flange dimensions that differ from Python `StandardBobbin` defaults. Need a policy decision on whether MVB++ should follow Python defaults or enriched MKF values for missing inputs.
3. **Object slicing in MKF integration** — `drawMagnetic(const MAS::Magnetic&)` was slicing enriched `OpenMagnetics::Magnetic` objects. Already fixed with an overload for `OpenMagnetics::Magnetic`.

## File map
| File | Purpose |
|------|---------|
| `CMakeLists.txt` | Build configuration, OCCT fetch + link |
| `include/mvb/Utils.h` / `src/Utils.cpp` | Dimension flattening, OCCT helpers, safe MKF wrapper |
| `include/mvb/MagneticBuilder.h` / `src/MagneticBuilder.cpp` | Full assembly builder (core + bobbin + turns) |
| `include/mvb/BobbinBuilder.h` / `src/BobbinBuilder.cpp` | Bobbin geometry builder |
| `include/mvb/TurnBuilder.h` / `src/TurnBuilder.cpp` | Turn geometry builder |
| `include/mvb/StepExporter.h` / `src/StepExporter.cpp` | STEP export with mm scaling |
| `src/shapes/Shape*.cpp` | Per-shape core builders |
| `tests/test_step_comparison.cpp` | Catch2 tests comparing MVB++ vs Python reference STEPs |
| `tests/generate_reference_steps.py` | Python script to regenerate reference STEPs from MVB |
| `tools/mvbpp_step_generator.cpp` | CLI for single/batch STEP generation |
| `tools/generate_mvbpp_bobbins.cpp` | Standalone bobbin comparison tool |
| `scripts/generate_mvb_bobbins.py` | Python reference bobbin generator |
| `BOBBIN_PROGRESS.md` | Detailed bobbin phase notes |
| `PLAN.md` | This file |
