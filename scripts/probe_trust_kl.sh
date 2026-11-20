#!/usr/bin/env bash
# Determine whether --trust-kl is "bound" during the P+ stage of pplus_internal.
# Only concerned with 02_sad/experimental.mtz, so use rounds=1 cycles=1 to run the shortest path.
# D: trust_kl=0.001   E: trust_kl=0.0005 (control, same parameters as production job A)
# If both 02_sad/experimental.mtz are bit-for-bit identical => the KL cap was never bound (saturated); the parameter was consumed but had no effect.
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/.conda-unified/bin/python"
CCP4="/Applications/ccp4-9"
IDIR="$ROOT/validation_ten_cases/inputs/7E1D"
OUTROOT="$ROOT/work/optimize"
LOGDIR="$OUTROOT/logs"
REF="${OASIS_REFERENCE:-$HOME/.cache/oasis-sad/reference/reference-1tqw.mtz}"
mkdir -p "$OUTROOT" "$LOGDIR"

BASE=(
  --phase-cycles 100 --dm-cycles 40 --free-value 0
  --aoe 0.5 --triplet-weight 0.1 --kmin 0.03 --max-neighbors 600
  --build-model --no-auto-hand --input-columns --building-recipe pipeline
  --no-canonical-dm-seed --no-per-round-cleanup
  --dm-engine parrot --dm-solvent 0.5 --dm-mask-radius 3.0 --dm-anisotropy
  --trust-cleanup --no-trust-stability
  --feedback-cycles 2 --trace-resolution 2.0 --feedback-endpoints 2,3 --auto-fpp
  --copies 2 --max-per-target 48 --elements Se --wavelength 0.97852
  --fpp "Se=3.8357741149160507" --fcalc fc
  --dm-executable "$CCP4/bin/cparrot" --ccp4-root "$CCP4"
  --sad "$IDIR/selected_sad.mtz" --native "$IDIR/native.mtz"
  --sequence "$IDIR/building.fasta" --heavy "$IDIR/refined_heavy.pdb"
  --method pplus_internal --refinement-prior original_sad --reference "$REF"
  --no-trust-partial-cleanup --no-trust-work-weight --trust-wang-radius 14
  --rounds 1 --cycles 1 --seed 123
)

for kl in 0.001 0.0005; do
  name="7E1D_klprobe_${kl}"
  out="$OUTROOT/$name"
  [ -e "$out" ] && { echo "[skip] $name"; continue; }
  ( cd "$ROOT" && "$PY" -u run_sad_to_model.py "${BASE[@]}" --trust-kl "$kl" --output "$out" ) \
      >"$LOGDIR/$name.log" 2>&1 &
  echo "[launch] $name pid=$!"
done

wait
echo "[all-done] KL probe completed $(date '+%H:%M:%S')"
