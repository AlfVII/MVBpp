#!/usr/bin/env bash
# Real-winding FEM STEP sweep over the whole MAS corpus (numbered examples + the
# `complete` fixtures). One process per design, capped memory, bounded parallelism.
#
#   scripts/sweep_real_steps.sh <outdir> [jobs]
#
# Every design produces a STEP to look at:
#   <name>.step             — built clean, collision gate satisfied, AND proven meshable.
#                             "FEM-grade" is a measured claim: see the MESH GATE below.
#   <name>.DIAGNOSTIC.step  — the gate REFUSED this geometry; exported anyway with
#                             MVB_LEAD_NO_VALIDATE so the fault can be seen. The copper
#                             in it interpenetrates. Never feed one of these to a mesher.
# RESULTS.tsv columns: name, build status, OVERLAP verdict, MESH verdict, seconds, reason.
# A design is FEM-grade only when all three read OK / CLEAN / MESHABLE.
#
# NB every design gets its OWN directory: mvbpp_step_generator always writes
# <outdir>/magnetic.step before renaming, so parallel jobs sharing one directory
# clobber each other's intermediate file.
#
# Env passed through: MVB_ALLOW_SQUISH, MVB_PATH_DUMP.
set -u

ROOT=/home/alf/OpenMagnetics/MVB++
GEN=$ROOT/build/mvbpp_step_generator
# The FEM product is FACETED (Alf, 2026-08-25: "for meshing it's alright to use segments=16").
# Exact round wire puts every conductor face on a cylinder/torus, and gmsh dispatches to its
# fragile periodic mesher off the SURFACE, not the face (OCCFace.cpp:142) -- 16 of the 23
# corpus mesh failures came from that path. The CORE stays exact regardless (the generator
# forces it under --fem) so faceting cannot bias its effective area.
SEGMENTS=${SEGMENTS:-0}
# Per-job virtual-memory cap. This box is SHARED with other agents' work, so a job that eats
# the machine OOM-kills somebody else (2026-08-25: a 36 GB generator of mine left 1.9 GB free and
# killed another agent's build). 16 GB is the standing budget and it is enough at BOTH segment
# settings now: the earlier need for 40 GB was the weld's growing fuse accumulator, not the
# faceted geometry, and gating the weld on real bridge growth removed it.
GEN_VMEM=${GEN_VMEM:-16000000}
STEPCHECK=${STEPCHECK:-/home/alf/OpenMagnetics/OMFEM/build/omfem_stepcheck}
INTERSECT=${INTERSECT:-/home/alf/OpenMagnetics/OMFEM/build/omfem_step_intersect}
OUT=${1:?usage: sweep_real_steps.sh <outdir> [jobs]}
# ONE GENERATOR AT A TIME (Alf, 2026-08-23: "never run multiple mvbpp_step_generator --real;
# if you launch another one, kill the previous ones"). This script used to default to THREE,
# which is how the box reached 0 GB available on 2026-08-23 -- one design
# (06_llc_xfmr_eq4128_3c97) took a single generator to 33 GB, and three at once is the exact
# condition that forced the 2026-07-27 OOM reboot. The per-job `ulimit -v 16G` below caps each
# process; it does not stop three of them adding up.
JOBS=${2:-1}
if [ "$JOBS" -gt 1 ] && [ -z "${MVB_ALLOW_PARALLEL_SWEEP:-}" ]; then
    echo "sweep_real_steps.sh: refusing jobs=$JOBS -- one generator at a time on this box." >&2
    echo "  A single design has been measured at 33 GB; concurrent sweeps took it to 0 GB free." >&2
    echo "  Set MVB_ALLOW_PARALLEL_SWEEP=1 to override deliberately." >&2
    exit 2
fi

mkdir -p "$OUT"
[ -n "${SWEEP_RESUME:-}" ] || : > "$OUT/RESULTS.tsv"

run_one() {
    local src="$1" name="$2" out="$3"
    local dir="$out/$name"
    # SWEEP_RESUME=1: a design already in RESULTS.tsv is done -- skip it. Long sweeps keep
    # getting killed by session teardown; without this every relaunch repays every completed
    # design (weld-all builds run ~10 min each).
    if [ -n "${SWEEP_RESUME:-}" ] && grep -q "^$name	" "$out/RESULTS.tsv" 2>/dev/null; then
        printf '[SKIP done] %s\n' "$name"
        return
    fi
    # Planar real windings are NOT implemented yet (Alf 2026-08-26: "09_planar should be
    # avoided, we are not implementing planar real windings YET") -- out of scope, not a fault.
    case "$name" in *planar*)
        printf '%s\tSKIPPED_PLANAR\t-\t-\t0s\tplanar real winding not implemented (ABT #492)\n' \
            "$name" >> "$out/RESULTS.tsv"
        printf '[SKIPPED_PLANAR] %s\n' "$name"
        return ;;
    esac
    # Start from an EMPTY directory. A design that used to build clean and now refuses (or the
    # reverse) would otherwise keep the previous run's artefact beside the current one, and the
    # pair gives no hint which is which except its timestamp -- 14_dab sat with a stale clean
    # .step next to a current .DIAGNOSTIC.step for a day.
    rm -rf "$dir"
    mkdir -p "$dir"
    local log="$dir/build.log"
    local step="$dir/$name.step"
    local t0=$SECONDS
    # --copper is NOT optional for the FEM product (Alf: MKF places by the OUTER envelope,
    # MVB++ paints only the CONDUCTING envelope, "so there will be a gap"). Without it the
    # paint fills the coated envelope and every MKF ride/rest (lead riser top on the bobbin
    # flange plane, first layer on the tube) becomes an EXACT-0 tangency that the fragment
    # turns into an unmeshable volume -- measured on 00_debug: riser top at z=-2.213 flush
    # with the bobbin bottom, "No elements in volume 7"; with --copper the enamel half-gap
    # (~22 um there) separates them.
    ( ulimit -v $GEN_VMEM; timeout 3600 "$GEN" --real --fem --copper --segments $SEGMENTS -o "$step" "$src" ) \
        > "$log" 2>&1
    local rc=$?
    local status reason=""
    if [ $rc -eq 0 ] && [ -s "$step" ]; then
        status=OK
    elif [ $rc -eq 124 ]; then
        status=TIMEOUT
    else
        reason=$(grep -m1 -a 'Error processing' "$log" \
                 | sed 's/.*\.json": //' | tr -d '\n' | cut -c1-400)
        [ -z "$reason" ] && reason=$(tail -2 "$log" | tr -d '\n' | cut -c1-400)
        # Export it anyway so the fault is visible, clearly marked as not FEM-grade.
        local diag="$dir/$name.DIAGNOSTIC.step"
        # MVB_MITRE_KEEP too: a refused mitre must still reach disk to be looked at, exactly like
        # a refused collision (Alf: "always produce step and 2D SVG together", even when it fails).
        ( ulimit -v $GEN_VMEM; MVB_LEAD_NO_VALIDATE=1 MVB_MITRE_KEEP=1 timeout 3600 \
              "$GEN" --real --fem --copper --segments $SEGMENTS -o "$diag" "$src" ) >> "$log" 2>&1
        if [ -s "$diag" ]; then status=REFUSED_DIAG; else status=FAIL; fi
    fi
    local dt=$((SECONDS - t0))
    # MESH GATE. "OK" used to mean only that a STEP reached disk, and the 2026-08-24 corpus
    # audit found 23 of 30 such "FEM-grade" files could not actually be meshed. A STEP that no
    # mesher will accept is not FEM-grade, so the sweep now has to PROVE it: import, fragment
    # into one conformal domain, tetrahedralise. Set SWEEP_NO_MESH_GATE=1 to build only.
    # OVERLAP GATE. Copper that interpenetrates copper (or the core) is a geometry defect
    # whatever the mesher does with it, and it is the one class a watertightness battery cannot
    # see: each solid is individually perfect. Runs before the mesh gate because an overlap
    # explains a mesh failure, not the other way round.
    local ovl=SKIPPED
    if [ "$status" = OK ] && [ -z "${SWEEP_NO_OVERLAP_GATE:-}" ] && [ -x "$INTERSECT" ]; then
        local olog="$dir/overlap.log"
        # OMFEM_INTERSECT_JOBS: the intersect tool now runs its fork-isolated pair booleans
        # through a pool -- the serial gate was ~50 min on a 178-solid toroid, which is what
        # produced every OVL_TIMEOUT row in the 2026-08-28 sweep.
        ( ulimit -v 50000000; OMFEM_INTERSECT_JOBS=${OMFEM_INTERSECT_JOBS:-4} \
              timeout "${OVERLAP_GATE_TIMEOUT:-3600}" "$INTERSECT" "$step" ) \
            > "$olog" 2>&1 < /dev/null
        case $? in
            0)   ovl=CLEAN ;;
            1)   ovl=OVERLAPS
                 reason=$(awk '/VERDICT: [0-9]/{sub(/^\[intersect\] VERDICT: /,""); print; exit}' "$olog" | cut -c1-200) ;;
            3)   ovl=INCONCLUSIVE ;;
            124) ovl=OVL_TIMEOUT ;;
            *)   ovl=OVL_ERROR ;;
        esac
    fi

    local mesh=SKIPPED
    if [ "$status" = OK ] && [ -z "${SWEEP_NO_MESH_GATE:-}" ] && [ -x "$STEPCHECK" ]; then
        local mlog="$dir/meshcheck.log"
        # --budget: stepcheck's own watchdog. It must sit ABOVE the honest runtime of the
        # slowest design that really does mesh, or the sweep reports a false NOT_MESHABLE:
        # 02_flyback (588 solids) meshes 588/588 in ~1473 s standalone and blew a 1800 s
        # budget here (2026-08-31) on ordinary run-to-run variance. 3000 s clears it with
        # margin and still fails fast on a genuinely hung boundary recovery, which runs for
        # hours. The outer `timeout` must stay above the budget so the watchdog is what fires.
        ( ulimit -v 50000000; timeout "${MESH_GATE_TIMEOUT:-3600}" \
              "$STEPCHECK" "$step" --mesh --budget "${MESH_GATE_BUDGET:-3000}" ) > "$mlog" 2>&1 < /dev/null
        case $? in
            0)   mesh=MESHABLE ;;
            124) mesh=MESH_TIMEOUT ;;
            *)   mesh=NOT_MESHABLE
                 reason=$(awk '/gmsh:/{sub(/^ *gmsh: /,""); print; exit}' "$mlog" | cut -c1-200)
                 [ -z "$reason" ] && reason=$(awk '/MESH /{print; exit}' "$mlog" | cut -c1-200) ;;
        esac
    fi
    printf '%s\t%s\t%s\t%s\t%ss\t%s\n' "$name" "$status" "$ovl" "$mesh" "$dt" "$reason" \
        >> "$out/RESULTS.tsv"
    printf '[%s/%s/%s] %s (%ss)\n' "$status" "$ovl" "$mesh" "$name" "$dt"
}
export -f run_one
export GEN SEGMENTS STEPCHECK INTERSECT GEN_VMEM SWEEP_RESUME

MAS=$ROOT/build/_deps/mas-src/examples
{
    for f in "$MAS"/[0-9][0-9]_*.json; do
        [ -e "$f" ] || continue
        b=$(basename "$f" .json)
        echo -e "$f\t$b"
    done
    for f in "$ROOT"/tests/mas_complete_fixtures/*.json; do
        [ -e "$f" ] || continue
        b="complete_$(basename "$f" .json)"
        echo -e "$f\t$b"
    done
} | xargs -P "$JOBS" -n 2 bash -c 'run_one "$0" "$1" "'"$OUT"'"'

echo "=== sweep done ==="
sort "$OUT/RESULTS.tsv"
