#!/usr/bin/env python3
"""Generate reference STEP files from MVB (Python/CadQuery) for MVB++ comparison tests."""

import sys
import os
import json
import copy
import math

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "MVB", "src", "OpenMagneticsVirtualBuilder"))

import cadquery as cq
from cadquery import exporters
from cadquery_builder import CadQueryBuilder, BobbinProcessedDescription, WireDescription, WireType, TurnDescription

# Force exact circular profile / torus geometry in the reference so MVB++ (which
# uses exact OCCT primitives) can be compared directly.
CadQueryBuilder.WIRE_POLYGON_SEGMENTS = 0

OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "reference_steps")
os.makedirs(OUTPUT_DIR, exist_ok=True)


def _build_proper_bobbin(bobbin_data, is_toroidal=False):
    """Build a proper bobbin matching MVB++ geometry using StandardBobbin."""
    if not isinstance(bobbin_data, dict):
        return None

    processed = bobbin_data.get("processedDescription", {})
    if not processed:
        return None

    std_bobbin = CadQueryBuilder().StandardBobbin()
    ww_list = processed.get("windingWindows", [{}])
    ww = ww_list[0] if ww_list else {}
    winding_window = {
        "width": ww.get("width", 0.01),
        "height": ww.get("height", 0.01),
        "radialHeight": ww.get("radialHeight", ww.get("width", 0.01)),
    }

    bobbin = std_bobbin.get_bobbin(bobbin_data, winding_window, name="Bobbin", save_files=False, export_files=False)
    if bobbin is None:
        return None

    if not is_toroidal:
        # Rotate to align with Y axis, matching MVB++ concentric core orientation
        bobbin = bobbin.rotate((1, 0, 0), (-1, 0, 0), -90)

    return bobbin


def _patch_cadquery_builder():
    """Monkey-patch get_magnetic to use dummy bobbins."""
    orig_get_magnetic = CadQueryBuilder.get_magnetic

    def patched_get_magnetic(self, magnetic_data, project_name="Magnetic", output_path=None, export_files=True, include_bobbin=True, include_residual_gaps=False):
        if output_path is None:
            output_path = f"{os.path.dirname(os.path.abspath(__file__))}/../../output/"

        core_pieces = []
        turn_pieces = []
        bobbin_geom = None
        is_toroidal = False

        core_data = magnetic_data.get("core", {})
        geometrical_description = core_data.get("geometricalDescription", [])
        for index, geometrical_part in enumerate(geometrical_description):
            if geometrical_part["type"] == "toroidal":
                is_toroidal = True
            if geometrical_part["type"] in ["half set", "toroidal"]:
                shape_data = geometrical_part["shape"]
                family = shape_data.get("family", "").lower()
                if family == "t":
                    is_toroidal = True
                part_builder = CadQueryBuilder().factory(shape_data)
                piece = part_builder.get_piece(data=copy.deepcopy(shape_data), name=f"Core_{index}", save_files=False, export_files=False)
                piece = piece.rotate((1, 0, 0), (-1, 0, 0), geometrical_part["rotation"][0] / math.pi * 180)
                piece = piece.rotate((0, 1, 0), (0, -1, 0), geometrical_part["rotation"][1] / math.pi * 180)
                piece = piece.rotate((0, 0, 1), (0, 0, -1), geometrical_part["rotation"][2] / math.pi * 180)
                if "machining" in geometrical_part and geometrical_part["machining"] is not None:
                    from utils import flatten_dimensions
                    for machining in geometrical_part["machining"]:
                        piece = part_builder.apply_machining(piece=piece, machining=machining, dimensions=flatten_dimensions(shape_data))
                piece = piece.translate((geometrical_part["coordinates"][0], geometrical_part["coordinates"][1], geometrical_part["coordinates"][2]))
                if include_residual_gaps and geometrical_part["type"] in ["half set"]:
                    residual_gap = 5e-6
                    if geometrical_part["rotation"][0] > 0:
                        piece = piece.translate((0, residual_gap / 2, 0))
                    else:
                        piece = piece.translate((0, -residual_gap / 2, 0))
                core_pieces.append(piece)

        coil_data = magnetic_data.get("coil", {})
        bobbin_data = coil_data.get("bobbin", {})
        bobbin_processed = BobbinProcessedDescription()
        if isinstance(bobbin_data, dict):
            bobbin_processed = BobbinProcessedDescription.from_dict(bobbin_data.get("processedDescription", {}))

        # Use proper bobbin geometry matching MVB++
        if include_bobbin:
            bobbin_geom = _build_proper_bobbin(bobbin_data, is_toroidal=is_toroidal)

        # Per-winding wire descriptors (multi-winding magnetics keep each
        # winding's wire separate and indexed by the turn's winding name).
        functional_desc = coil_data.get("functionalDescription", [])
        wires_by_winding = {}
        default_wire_desc = WireDescription(WireType.round)
        for fd in functional_desc:
            name = fd.get("name")
            wire_data = fd.get("wire", {})
            if isinstance(wire_data, dict):
                wd = WireDescription.from_dict(wire_data)
            else:
                wd = WireDescription(WireType.round)
            if name:
                wires_by_winding[name] = wd
            if default_wire_desc.wire_type == WireType.round and fd == functional_desc[0]:
                default_wire_desc = wd

        turns_data = coil_data.get("turnsDescription", [])
        for turn_data in turns_data:
            turn_desc = TurnDescription.from_dict(turn_data)
            # Pick the wire description for this turn's winding so the real
            # wire type (round/rectangular) is preserved. Falls back to the
            # first winding's wire if the lookup fails.
            winding_name = turn_data.get("winding")
            wire_desc = wires_by_winding.get(winding_name, default_wire_desc)
            turn_geom = self.get_turn(turn_desc, wire_desc, bobbin_processed, is_toroidal=is_toroidal)
            turn_pieces.append(turn_geom)

        if bobbin_geom is not None:
            for core_piece in core_pieces:
                try:
                    bobbin_geom = bobbin_geom.cut(core_piece)
                except Exception:
                    pass
            for turn in turn_pieces:
                try:
                    bobbin_geom = bobbin_geom.cut(turn)
                except Exception:
                    pass

        all_pieces = core_pieces[:]
        if bobbin_geom is not None:
            all_pieces.append(bobbin_geom)
        all_pieces.extend(turn_pieces)

        if export_files and all_pieces:
            scaled_pieces = []
            for piece in all_pieces:
                for o in piece.objects:
                    scaled_pieces.append(o.scale(1000))
            compound = cq.Compound.makeCompound(scaled_pieces)
            step_path = f"{output_path}/{project_name}.step"
            stl_path = f"{output_path}/{project_name}.stl"
            exporters.export(compound, step_path, "STEP")
            return step_path, None
        elif all_pieces:
            return all_pieces
        else:
            return None, None

    CadQueryBuilder.get_magnetic = patched_get_magnetic


def build_magnetic_reference(json_path, out_name):
    with open(json_path, "r") as f:
        data = json.load(f)

    magnetic_data = data.get("magnetic", {})
    builder = CadQueryBuilder()

    step_path = builder.get_magnetic(
        magnetic_data=magnetic_data,
        project_name=out_name,
        output_path=OUTPUT_DIR,
        export_files=True,
        include_bobbin=True,
        include_residual_gaps=False,
    )
    if step_path and step_path[0]:
        print(f"Reference assembly: {step_path[0]}")
        return step_path[0]
    else:
        print(f"Failed to generate reference for {out_name}")
        return None


def build_e_core_reference():
    shape_data = {
        "name": "E 19/8/5",
        "family": "e",
        "type": "standard",
        "dimensions": {
            "A": 0.019,
            "B": 0.008,
            "C": 0.005,
            "D": 0.005,
            "E": 0.012,
            "F": 0.006,
        },
    }
    geometrical_description = [
        {"type": "half set", "coordinates": [0.0, 0.0, 0.0], "rotation": [0.0, 0.0, 0.0], "shape": shape_data},
        {"type": "half set", "coordinates": [0.0, 0.0, 0.0], "rotation": [3.141592653589793, 0.0, 0.0], "shape": shape_data},
    ]
    builder = CadQueryBuilder()
    step_path, _ = builder.get_core("reference_e_core", geometrical_description, output_path=OUTPUT_DIR, save_files=False, export_files=True)
    print(f"E-core reference: {step_path}")


def build_t_core_reference():
    shape_data = {
        "name": "T 20/10/8",
        "family": "t",
        "type": "standard",
        "dimensions": {"A": 0.020, "B": 0.010, "C": 0.008},
    }
    geometrical_description = [
        {"type": "toroidal", "coordinates": [0.0, 0.0, 0.0], "rotation": [0.0, 0.0, 0.0], "shape": shape_data},
    ]
    builder = CadQueryBuilder()
    step_path, _ = builder.get_core("reference_t_core", geometrical_description, output_path=OUTPUT_DIR, save_files=False, export_files=True)
    print(f"T-core reference: {step_path}")


if __name__ == "__main__":
    _patch_cadquery_builder()
    build_e_core_reference()
    build_t_core_reference()

    test_data_dir = os.path.join(os.path.dirname(__file__), "..", "..", "MVB", "tests", "testData")
    build_magnetic_reference(os.path.join(test_data_dir, "concentric_rectangular_column_one_turn.json"), "reference_rect_one_turn")
    build_magnetic_reference(os.path.join(test_data_dir, "ETD49_N87_10uH_5T.json"), "reference_etd49_5t")
