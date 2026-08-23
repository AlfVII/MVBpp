#!/usr/bin/env bash
# Real-winding FEM STEP sweep over the whole MAS corpus (numbered examples + the
# `complete` fixtures). One process per design, capped memory, bounded parallelism.
#
#   scripts/sweep_real_steps.sh <outdir> [jobs]
#
# Every design produces a STEP to look at:
#   <name>.step             — built clean, collision gate satisfied. FEM-grade.
#   <name>.DIAGNOSTIC.step  — the gate REFUSED this geometry; exported anyway with
#                             MVB_LEAD_NO_VALIDATE so the fault can be seen. The copper
#                             in it interpenetrates. Never feed one of these to a mesher.
# RESULTS.tsv carries the verdict and, for a refusal, the reason verbatim.
#
# NB every design gets its OWN directory: mvbpp_step_generator always writes
# <outdir>/magnetic.step before renaming, so parallel jobs sharing one directory
# clobber each other's intermediate file.
#
# Env passed through: MVB_ALLOW_SQUISH, MVB_PATH_DUMP.
set -u

ROOT=/home/alf/OpenMagnetics/MVB++
GEN=$ROOT/build/mvbpp_step_generator
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
: > "$OUT/RESULTS.tsv"

run_one() {
    local src="$1" name="$2" out="$3"
    local dir="$out/$name"
    # Start from an EMPTY directory. A design that used to build clean and now refuses (or the
    # reverse) would otherwise keep the previous run's artefact beside the current one, and the
    # pair gives no hint which is which except its timestamp -- 14_dab sat with a stale clean
    # .step next to a current .DIAGNOSTIC.step for a day.
    rm -rf "$dir"
    mkdir -p "$dir"
    local log="$dir/build.log"
    local step="$dir/$name.step"
    local t0=$SECONDS
    ( ulimit -v 16000000; timeout 3600 "$GEN" --real --fem --segments 0 -o "$step" "$src" ) \
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
        ( ulimit -v 16000000; MVB_LEAD_NO_VALIDATE=1 MVB_MITRE_KEEP=1 timeout 3600 \
              "$GEN" --real --fem --segments 0 -o "$diag" "$src" ) >> "$log" 2>&1
        if [ -s "$diag" ]; then status=REFUSED_DIAG; else status=FAIL; fi
    fi
    local dt=$((SECONDS - t0))
    printf '%s\t%s\t%ss\t%s\n' "$name" "$status" "$dt" "$reason" >> "$out/RESULTS.tsv"
    printf '[%s] %s (%ss)\n' "$status" "$name" "$dt"
}
export -f run_one
export GEN

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
