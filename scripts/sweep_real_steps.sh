#!/usr/bin/env bash
# Real-winding FEM STEP sweep over the whole MAS corpus (numbered examples + the
# `complete` fixtures). One process per design, capped memory, bounded parallelism —
# the point is a PASS/FAIL table with the real reason, not a pile of STEPs.
#
#   scripts/sweep_real_steps.sh <outdir> [jobs]
#
# NB every design gets its OWN directory: mvbpp_step_generator always writes
# <outdir>/magnetic.step before renaming, so parallel jobs sharing one directory
# clobber each other's intermediate file.
#
# Env passed through: MVB_ALLOW_SQUISH, MVB_LEAD_NO_VALIDATE, MVB_PATH_DUMP.
set -u

ROOT=/home/alf/OpenMagnetics/MVB++
GEN=$ROOT/build/mvbpp_step_generator
OUT=${1:?usage: sweep_real_steps.sh <outdir> [jobs]}
JOBS=${2:-3}

mkdir -p "$OUT"
: > "$OUT/RESULTS.tsv"

run_one() {
    local src="$1" name="$2" out="$3"
    local dir="$out/$name"
    mkdir -p "$dir"
    local log="$dir/build.log"
    local step="$dir/$name.step"
    local t0=$SECONDS
    ( ulimit -v 16000000; timeout 3600 "$GEN" --real --fem --segments 0 -o "$step" "$src" ) \
        > "$log" 2>&1
    local rc=$?
    local dt=$((SECONDS - t0))
    local status reason=""
    if [ $rc -eq 0 ] && [ -s "$step" ]; then
        status=OK
    elif [ $rc -eq 124 ]; then
        status=TIMEOUT
    else
        status=FAIL
        # The thrown reason, on one line, with the leading noise stripped.
        reason=$(grep -m1 -a 'Error processing' "$log" \
                 | sed 's/.*\.json": //' | tr -d '\n' | cut -c1-220)
        [ -z "$reason" ] && reason=$(tail -2 "$log" | tr -d '\n' | cut -c1-220)
    fi
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
