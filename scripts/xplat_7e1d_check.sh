#!/usr/bin/env bash
# macOS-side cross-validation: reproduce the best configuration from the Linux-side 7E1D R-free optimization report (2026-11-19).
# A: exact reproduction of the Linux best configuration (pplus_internal, r14, c15, seed 123)
# B: pure seed variable (pplus_bundled, r5, c10, seed 123), compared with the macOS baseline s20261023 = 0.4167
# C: internal P+ + default seed (pplus_internal, r14, c15, seed 20261023), to separate method effects from seed effects
set -uo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PY="/Volumes/Fast_4T/GPT_WORK/OASIS4.3-Build/.conda-unified/bin/python"
CCP4="/Applications/ccp4-9"
IDIR="$ROOT/validation_ten_cases/inputs/7E1D"
OUTROOT="$ROOT/work/optimize"
LOGDIR="$OUTROOT/logs"
# The internal MLHL DM of pplus_internal needs a reference histogram (FC_BEST.F_phi.F/phi),
# This file is pinned by sha256 and is bit-for-bit identical across platforms, so it is not a source of cross-platform differences.
REF="${OASIS_REFERENCE:-$HOME/.cache/oasis-sad/reference/reference-1tqw.mtz}"
mkdir -p "$OUTROOT" "$LOGDIR"

COMMON=(
  --phase-cycles 100 --dm-cycles 40 --free-value 0
  --aoe 0.5 --triplet-weight 0.1 --kmin 0.03 --max-neighbors 600
  --build-model --no-auto-hand --input-columns
  --building-recipe pipeline
  --no-canonical-dm-seed --no-per-round-cleanup
  --dm-engine parrot --dm-solvent 0.5 --dm-mask-radius 3.0 --dm-anisotropy
  --trust-kl 0.0005 --trust-cleanup --no-trust-stability
  --feedback-cycles 2 --trace-resolution 2.0 --feedback-endpoints 2,3
  --auto-fpp
  --copies 2 --max-per-target 48 --elements Se --wavelength 0.97852
  --fpp "Se=3.8357741149160507" --fcalc fc
  --dm-executable "$CCP4/bin/cparrot" --ccp4-root "$CCP4"
  --sad "$IDIR/selected_sad.mtz" --native "$IDIR/native.mtz"
  --sequence "$IDIR/building.fasta" --heavy "$IDIR/refined_heavy.pdb"
)

launch() {
  local name="$1"; shift
  local out="$OUTROOT/$name"
  if [ -e "$out" ]; then echo "[skip] $name already exists"; return; fi
  ( cd "$ROOT" && "$PY" -u run_sad_to_model.py "${COMMON[@]}" "$@" --output "$out" ) \
      >"$LOGDIR/$name.log" 2>&1 &
  echo "[launch] $name pid=$! -> $LOGDIR/$name.log"
}

# A: Linux best configuration. Under pplus_internal, refinement_prior must be original_sad (dm_phifom is bundled-only),
#    trust-partial/work-weight also does not apply to the internal path.
#    The Wang radius only takes effect for trust/bundled; internal must use 14 (otherwise validate rejects it directly).
launch 7E1D_internal_r14_c15_s123 \
  --method pplus_internal --refinement-prior original_sad --reference "$REF" \
  --no-trust-partial-cleanup --no-trust-work-weight --trust-wang-radius 14 \
  --rounds 14 --cycles 15 --seed 123

# B: only the seed is changed (everything else is identical to the macOS baseline 0.4167)
launch 7E1D_bundled_r5_c10_s123 \
  --method pplus_bundled --refinement-prior dm_phifom \
  --trust-partial-cleanup --trust-work-weight --trust-wang-radius 8 \
  --rounds 5 --cycles 10 --seed 123

# C: internal P+ + default seed, to separate method effects from seed effects
launch 7E1D_internal_r14_c15_s20261023 \
  --method pplus_internal --refinement-prior original_sad --reference "$REF" \
  --no-trust-partial-cleanup --no-trust-work-weight --trust-wang-radius 14 \
  --rounds 14 --cycles 15 --seed 20261023

echo "[done] Three jobs launched in parallel, waiting for all to complete"
wait
echo "[all-done] All three jobs finished $(date '+%Y-%m-%d %H:%M:%S')"
