# Custom-core MAS variants

Design variants that differ from the stock MAS example ONLY in the core shape, kept out of
`_deps/mas-src/examples` (a vendored checkout) and out of the sweep corpus so they change no
existing result.

## 13_current_sense_er95_widecore.json

`13_current_sense_er95_n87` on stock **ER 9.5/2.5/5** puts 1.157 mm3 of BARE COPPER inside both
core halves (99 pairs, worst 0.0263 mm3) — verified pre-existing, identical before and after the
2026-08-25 geometry work, and the reason that design lost -3.4% volume in the mesh fragment.
Copper reaches x = +-5.233 mm against a window wall at 3.8125 mm.

This variant widens the window until the audit is clean (`omfem_step_intersect`:
`VERDICT: NO OVERLAPS`). Per Alf, 2026-08-25: larger A and E, column and heights untouched.

| dim | stock | here | note |
|-----|-------|------|------|
| A   | 9.35  | 14.325 | outer width |
| E   | 7.625 | 12.60  | window width |
| G   | 7.2   | 10.20  | HAD to change: MKF rejects the shape otherwise ("Shape constants cannot be negative or 0") — G is a window dimension for `planarER` and must scale with E |
| B, C, D, F | 2.455 / 4.9 / 1.675 / 3.4 | unchanged | height, depth, window height, COLUMN diameter |

Took three iterations: the winding spreads as the window grows (copper reached 5.233 -> 5.461 mm),
so the wall has to clear it with margin rather than meet it.
