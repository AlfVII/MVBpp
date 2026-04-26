"""
mvbpp Python binding tests.

Run with:  pytest bindings/python/tests/ -v
The mvbpp wheel must be installed first:
    pip install dist/mvbpp-*.whl
"""
import json
import os
import struct
import tempfile
from pathlib import Path

import pytest

import mvbpp

TESTS_DIR = Path(__file__).parent


# ── helpers ──────────────────────────────────────────────────────────────────

def load_magnetic_json(filename: str) -> str:
    """Return the 'magnetic' object from a test data file as a JSON string."""
    path = TESTS_DIR / filename
    data = json.loads(path.read_text())
    return json.dumps(data["magnetic"])


def is_step_bytes(data: bytes) -> bool:
    """Minimal check: STEP files start with 'ISO-10303-21'."""
    return data[:12] == b"ISO-10303-21"


# ── draw_magnetic_to_bytes ────────────────────────────────────────────────────

class TestDrawMagneticToBytes:

    def test_returns_bytes(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j)
        assert isinstance(result, bytes)

    def test_returns_nonempty(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j)
        assert len(result) > 1000

    def test_is_valid_step_header(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j)
        assert is_step_bytes(result), "Output does not start with ISO-10303-21"

    def test_step_contains_geometry_data(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j)
        text = result.decode("ascii", errors="replace")
        # STEP files with solid geometry contain CLOSED_SHELL or ADVANCED_BREP
        assert "CLOSED_SHELL" in text or "ADVANCED_BREP" in text

    def test_etd49_returns_valid_step(self):
        j = load_magnetic_json("ETD49_N87_10uH_5T.json")
        result = mvbpp.draw_magnetic_to_bytes(j)
        assert is_step_bytes(result)
        assert len(result) > 10_000

    def test_etd49_larger_than_basic(self):
        basic = mvbpp.draw_magnetic_to_bytes(load_magnetic_json("concentric_basic.json"))
        etd   = mvbpp.draw_magnetic_to_bytes(load_magnetic_json("ETD49_N87_10uH_5T.json"))
        # ETD49 with 5 turns should produce more geometry than the basic example
        assert len(etd) > len(basic)


# ── draw_magnetic (to file) ───────────────────────────────────────────────────

class TestDrawMagneticToFile:

    def test_writes_file(self):
        j = load_magnetic_json("concentric_basic.json")
        with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
            path = f.name
        try:
            mvbpp.draw_magnetic(j, path)
            assert os.path.exists(path)
            assert os.path.getsize(path) > 1000
        finally:
            os.unlink(path)

    def test_written_file_is_valid_step(self):
        j = load_magnetic_json("concentric_basic.json")
        with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
            path = f.name
        try:
            mvbpp.draw_magnetic(j, path)
            data = Path(path).read_bytes()
            assert is_step_bytes(data)
        finally:
            os.unlink(path)

    def test_file_matches_bytes_output(self):
        j = load_magnetic_json("concentric_basic.json")
        with tempfile.NamedTemporaryFile(suffix=".step", delete=False) as f:
            path = f.name
        try:
            mvbpp.draw_magnetic(j, path)
            file_data  = Path(path).read_bytes()
            bytes_data = mvbpp.draw_magnetic_to_bytes(j)
            # Both paths should produce STEP of similar size (within 1%)
            ratio = abs(len(file_data) - len(bytes_data)) / max(len(file_data), len(bytes_data))
            assert ratio < 0.01, f"File vs bytes size differs by {ratio:.1%}"
        finally:
            os.unlink(path)


# ── draw parameters (no config object) ────────────────────────────────────────

class TestDrawParameters:

    def test_default_parameters_produce_output(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j)
        assert is_step_bytes(result)

    def test_no_bobbin_produces_smaller_output(self):
        j = load_magnetic_json("concentric_basic.json")
        with_bobbin    = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=True)
        without_bobbin = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=False)
        assert len(with_bobbin) > len(without_bobbin)

    def test_scale_affects_output(self):
        """Scale=1000 (mm) should be accepted without error."""
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j, scale=1000.0)
        assert is_step_bytes(result)

    def test_wire_polygon_segments_accepted(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j, wire_polygon_segments=8)
        assert is_step_bytes(result)

    def test_core_polygon_segments_accepted(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j, core_polygon_segments=32)
        assert is_step_bytes(result)

    def test_no_draw_config_class(self):
        """DrawConfig must not be exposed — configuration goes via kwargs only."""
        assert not hasattr(mvbpp, "DrawConfig")


# ── 2D dimensioned drawings ───────────────────────────────────────────────────

class TestDimensionedDrawings:

    def test_front_view_returns_svg(self):
        j = load_magnetic_json("concentric_basic.json")
        svg = mvbpp.draw_dimensioned_front_view(j)
        assert isinstance(svg, str)
        assert svg.lstrip().startswith("<?xml") or svg.lstrip().startswith("<svg")

    def test_top_view_returns_svg(self):
        j = load_magnetic_json("concentric_basic.json")
        svg = mvbpp.draw_dimensioned_top_view(j)
        assert isinstance(svg, str)
        assert svg.lstrip().startswith("<?xml") or svg.lstrip().startswith("<svg")

    def test_front_view_custom_colors(self):
        j = load_magnetic_json("concentric_basic.json")
        svg = mvbpp.draw_dimensioned_front_view(
            j, width_px=600.0, label_font_px=10.0,
            projection_color="#222222", dimension_color="#ff0000")
        assert "#ff0000" in svg or "ff0000" in svg.lower()

    def test_write_front_view_to_file(self):
        j = load_magnetic_json("concentric_basic.json")
        with tempfile.NamedTemporaryFile(suffix=".svg", delete=False) as f:
            path = f.name
        try:
            mvbpp.write_dimensioned_front_view(j, path)
            data = Path(path).read_text()
            assert "<svg" in data
        finally:
            os.unlink(path)


# ── symmetry ──────────────────────────────────────────────────────────────────

class TestSymmetry:

    def test_get_symmetry_planes_returns_list(self):
        j = load_magnetic_json("concentric_basic.json")
        planes = mvbpp.get_symmetry_planes(j)
        assert isinstance(planes, list)

    def test_e_core_has_symmetry_planes(self):
        j = load_magnetic_json("concentric_basic.json")
        planes = mvbpp.get_symmetry_planes(j)
        assert len(planes) >= 1, f"Expected at least 1 symmetry plane, got {planes}"

    def test_symmetry_plane_names_are_valid(self):
        j = load_magnetic_json("concentric_basic.json")
        planes = mvbpp.get_symmetry_planes(j)
        valid = {"XY", "YZ", "XZ"}
        for p in planes:
            assert p in valid, f"Unexpected plane name: {p}"

    def test_one_plane_produces_smaller_step(self):
        j = load_magnetic_json("concentric_basic.json")
        full = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=False, symmetry_planes=0)
        half = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=False, symmetry_planes=1)
        assert len(half) < len(full), "Half-symmetry STEP should be smaller than full"

    def test_two_planes_smaller_than_one_plane(self):
        j = load_magnetic_json("concentric_basic.json")
        planes = mvbpp.get_symmetry_planes(j)
        if len(planes) < 2:
            pytest.skip("Magnetic has fewer than 2 symmetry planes")
        half    = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=False, symmetry_planes=1)
        quarter = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=False, symmetry_planes=2)
        assert len(quarter) < len(half), "Quarter-symmetry STEP should be smaller than half"

    def test_symmetry_output_is_valid_step(self):
        j = load_magnetic_json("concentric_basic.json")
        result = mvbpp.draw_magnetic_to_bytes(j, include_bobbin=False, symmetry_planes=1)
        assert is_step_bytes(result)


# ── error handling ────────────────────────────────────────────────────────────

class TestErrorHandling:

    def test_invalid_json_raises(self):
        with pytest.raises(Exception):
            mvbpp.draw_magnetic_to_bytes("not json at all")

    def test_empty_json_object_raises(self):
        with pytest.raises(Exception):
            mvbpp.draw_magnetic_to_bytes("{}")
