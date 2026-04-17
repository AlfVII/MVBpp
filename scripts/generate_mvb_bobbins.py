#!/usr/bin/env python3
"""Generate simple bobbin STEP files from Python MVB for comparison."""

import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "MVB", "src", "OpenMagneticsVirtualBuilder"))

import cadquery as cq
from cadquery import exporters

OUTPUT_DIR = "/home/alf/OpenMagnetics/MVB++/bobbin_comparison"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Rectangular bobbin parameters (meters)
col_width = 0.005
col_depth = 0.004
ww_width = 0.003
ww_depth = ww_width
ww_height = 0.010
wall_thick = 0.001

# Build rectangular bobbin like MVB does
outer_width = ww_width + 2 * wall_thick
outer_depth = ww_depth + 2 * wall_thick
height = ww_height

hole_width = ww_width * 0.8
hole_depth = ww_depth * 0.8

oCorner = (-outer_width / 2.0, -outer_depth / 2.0, -height / 2.0)
iCorner = (-ww_width / 2.0, -ww_depth / 2.0, -height / 2.0 - 0.0005)
cCorner = (-hole_width / 2.0, -hole_depth / 2.0, -height / 2.0 - 0.001)

outer = cq.Workplane("XY").box(outer_width, outer_depth, height, centered=False).translate(oCorner)
inner = cq.Workplane("XY").box(ww_width, ww_depth, height * 1.001, centered=False).translate(iCorner)
central = cq.Workplane("XY").box(hole_width, hole_depth, height * 1.002, centered=False).translate(cCorner)

body = outer.cut(inner).cut(central)

flangeThickness = 0.001
flangeExtension = 0.002
flangeWidth = outer_width + flangeExtension * 2.0
flangeDepth = outer_depth + flangeExtension * 2.0

tfCorner = (-flangeWidth / 2.0, -flangeDepth / 2.0, height / 2.0)
topFlange = cq.Workplane("XY").box(flangeWidth, flangeDepth, flangeThickness, centered=False).translate(tfCorner)
topHole = cq.Workplane("XY").box(hole_width, hole_depth, flangeThickness * 1.1, centered=False).translate((-hole_width/2, -hole_depth/2, height/2 - 0.0005))
topFlange = topFlange.cut(topHole)

bfCorner = (-flangeWidth / 2.0, -flangeDepth / 2.0, -(height / 2.0 + flangeThickness))
bottomFlange = cq.Workplane("XY").box(flangeWidth, flangeDepth, flangeThickness, centered=False).translate(bfCorner)
bottomHole = cq.Workplane("XY").box(hole_width, hole_depth, flangeThickness * 1.1, centered=False).translate((-hole_width/2, -hole_depth/2, -(height/2 + flangeThickness) - 0.0005))
bottomFlange = bottomFlange.cut(bottomHole)

rect_bobbin = body.union(topFlange).union(bottomFlange)

exporters.export(rect_bobbin, os.path.join(OUTPUT_DIR, "mvb_rectangular_bobbin.step"))
print(f"MVB rectangular bobbin exported")

# Round bobbin parameters
outer_r = col_width
inner_r = col_depth - wall_thick

outer_cyl = cq.Workplane("XY").cylinder(height, outer_r)
inner_cyl = cq.Workplane("XY").cylinder(height * 1.001, inner_r)
body = outer_cyl.cut(inner_cyl)

flangeOuterX = outer_r + ww_width + flangeExtension
flangeHalfY = col_depth / 2.0 if col_depth > 0 else outer_r

tfCorner = (-flangeOuterX, -flangeHalfY, height / 2.0)
topFlange = cq.Workplane("XY").box(flangeOuterX * 2.0, flangeHalfY * 2.0, flangeThickness, centered=False).translate(tfCorner)
topHole = cq.Workplane("XY").cylinder(flangeThickness * 1.1, inner_r).translate((0, 0, height / 2.0))
topFlange = topFlange.cut(topHole)

bfCorner = (-flangeOuterX, -flangeHalfY, -(height / 2.0 + flangeThickness))
bottomFlange = cq.Workplane("XY").box(flangeOuterX * 2.0, flangeHalfY * 2.0, flangeThickness, centered=False).translate(bfCorner)
bottomHole = cq.Workplane("XY").cylinder(flangeThickness * 1.1, inner_r).translate((0, 0, -(height / 2.0 + flangeThickness * 1.1)))
bottomFlange = bottomFlange.cut(bottomHole)

round_bobbin = body.union(topFlange).union(bottomFlange)

exporters.export(round_bobbin, os.path.join(OUTPUT_DIR, "mvb_round_bobbin.step"))
print(f"MVB round bobbin exported")

print(f"\nMVB Rectangular bobbin bbox:")
bb = rect_bobbin.val().BoundingBox()
print(f"  x: [{bb.xmin*1000:.3f}, {bb.xmax*1000:.3f}] mm")
print(f"  y: [{bb.ymin*1000:.3f}, {bb.ymax*1000:.3f}] mm")
print(f"  z: [{bb.zmin*1000:.3f}, {bb.zmax*1000:.3f}] mm")

print(f"\nMVB Round bobbin bbox:")
bb = round_bobbin.val().BoundingBox()
print(f"  x: [{bb.xmin*1000:.3f}, {bb.xmax*1000:.3f}] mm")
print(f"  y: [{bb.ymin*1000:.3f}, {bb.ymax*1000:.3f}] mm")
print(f"  z: [{bb.zmin*1000:.3f}, {bb.zmax*1000:.3f}] mm")
