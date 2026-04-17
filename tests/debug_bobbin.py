import sys, os, json, math, types
sys.path.insert(0, '../MVB/src/OpenMagneticsVirtualBuilder')
from cadquery_builder import CadQueryBuilder, BobbinProcessedDescription, WireDescription, WireType, TurnDescription

# Monkey-patch get_magnetic to pass full bobbin_data dict to get_bobbin
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
            piece = part_builder.get_piece(data=dict(shape_data), name=f"Core_{index}", save_files=False, export_files=False)
            piece = piece.rotate((1,0,0), (-1,0,0), geometrical_part['rotation'][0] / math.pi * 180)
            piece = piece.rotate((0,1,0), (0,-1,0), geometrical_part['rotation'][1] / math.pi * 180)
            piece = piece.rotate((0,0,1), (0,0,-1), geometrical_part['rotation'][2] / math.pi * 180)
            if 'machining' in geometrical_part and geometrical_part['machining'] is not None:
                import copy
                from utils import flatten_dimensions
                for machining in geometrical_part['machining']:
                    piece = part_builder.apply_machining(piece=piece, machining=machining, dimensions=flatten_dimensions(shape_data))
            piece = piece.translate((geometrical_part['coordinates'][0], geometrical_part['coordinates'][1], geometrical_part['coordinates'][2]))
            if include_residual_gaps and geometrical_part['type'] in ['half set']:
                residual_gap = 5e-6
                if geometrical_part['rotation'][0] > 0:
                    piece = piece.translate((0, residual_gap / 2, 0))
                else:
                    piece = piece.translate((0, -residual_gap / 2, 0))
            core_pieces.append(piece)
    
    coil_data = magnetic_data.get("coil", {})
    bobbin_data = coil_data.get("bobbin", {})
    if isinstance(bobbin_data, str):
        bobbin_processed = CadQueryBuilder.BobbinProcessedDescription()
    else:
        bobbin_processed_data = bobbin_data.get("processedDescription", {})
        bobbin_processed = BobbinProcessedDescription.from_dict(bobbin_processed_data)
    
    # FIX: pass full bobbin_data dict instead of bobbin_processed
    if not is_toroidal and include_bobbin:
        bobbin_geom = self.get_bobbin(bobbin_data, winding_window={}, export_files=False)
    
    wire_desc = WireDescription(WireType.round)
    functional_desc = coil_data.get("functionalDescription", [])
    if functional_desc:
        wire_data = functional_desc[0].get("wire", {})
        if wire_data and isinstance(wire_data, dict):
            wire_desc = WireDescription.from_dict(wire_data)
    
    turns_data = coil_data.get("turnsDescription", [])
    for turn_data in turns_data:
        turn_desc = TurnDescription.from_dict(turn_data)
        if turn_data.get("dimensions"):
            dims = turn_data["dimensions"]
            if len(dims) >= 2:
                wire_desc = WireDescription(
                    wire_type=WireType.round if turn_data.get("crossSectionalShape", "round") == "round" else WireType.rectangular,
                    outer_diameter=dims[0], conducting_diameter=dims[0]
                )
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
        import cadquery as cq
        from cadquery import exporters
        scaled_pieces = []
        for piece in all_pieces:
            for o in piece.objects:
                scaled_pieces.append(o.scale(1000))
        compound = cq.Compound.makeCompound(scaled_pieces)
        step_path = f"{output_path}/{project_name}.step"
        stl_path = f"{output_path}/{project_name}.stl"
        exporters.export(compound, step_path, "STEP")
        exporters.export(compound, stl_path, "STL")
        return step_path, stl_path
    elif all_pieces:
        return all_pieces
    else:
        return None, None

CadQueryBuilder.get_magnetic = patched_get_magnetic

# Now test
builder = CadQueryBuilder()
with open("../MVB/tests/testData/concentric_rectangular_column_one_turn.json") as f:
    data = json.load(f)
result = builder.get_magnetic(data.get("magnetic", {}), "test", export_files=False, include_bobbin=True, include_residual_gaps=False)
print(f"rect_one_turn pieces: {len(result)}")
for i, p in enumerate(result):
    bb = p.val().BoundingBox()
    print(f"  {i}: solids={len(list(p.val().Solids()))} bbox=[{bb.xmin:.4f}, {bb.ymin:.4f}, {bb.zmin:.4f}] - [{bb.xmax:.4f}, {bb.ymax:.4f}, {bb.zmax:.4f}]")

with open("../MVB/tests/testData/ETD49_N87_10uH_5T.json") as f:
    data = json.load(f)
result = builder.get_magnetic(data.get("magnetic", {}), "test", export_files=False, include_bobbin=True, include_residual_gaps=False)
print(f"etd49_5t pieces: {len(result)}")
for i, p in enumerate(result):
    bb = p.val().BoundingBox()
    print(f"  {i}: solids={len(list(p.val().Solids()))} bbox=[{bb.xmin:.4f}, {bb.ymin:.4f}, {bb.zmin:.4f}] - [{bb.xmax:.4f}, {bb.ymax:.4f}, {bb.zmax:.4f}]")
