#!/usr/bin/env bash
# End-to-end check of the in-house engines on the 6ZPP calibration case.
# handdm replaces the external CCP4 Parrot density modification;
# minirefmac replaces the external CCP4 REFMAC refinement.
set -euo pipefail
cd "$(dirname "$0")/../.."   # repository root

PY="${PHASEBUILD_PYTHON:-/opt/miniconda3/bin/python3.12}"
INS=validation_ten_cases/inputs/6ZPP
OUT=work/inhouse/6ZPP_parrot_legacy
ROUNDS="${ROUNDS:-2}"
CYCLES="${CYCLES:-10}"

# Refuse to overwrite an existing run directory instead of deleting it:
# the previous results are the comparison baseline and must not be destroyed.
if [ -e "$OUT" ]; then
  echo "output directory already exists: $OUT" >&2
  echo "set OUT to a fresh directory to run again" >&2
  exit 2
fi

export PHASEBUILD_PYTHON="$PY"
"$PY" -u run_sad_to_model.py \
  --phase-cycles 100 --dm-cycles 40 --rounds "$ROUNDS" --cycles "$CYCLES" --free-value 0 \
  --seed 20261023 --aoe 0.5 --triplet-weight 0.1 --kmin 0.03 --max-neighbors 600 \
  --build-model --no-auto-hand --input-columns \
  --building-recipe legacy \
  --no-canonical-dm-seed --no-per-round-cleanup \
  --dm-engine parrot --dm-solvent 0.5 --dm-mask-radius 3.0 --dm-anisotropy \
  --trust-kl 0.0005 --trust-cleanup --no-trust-stability --trust-wang-radius 8.0 \
  --feedback-cycles 0 \
  --refine-engine refmac \
  --method pplus_bundled --copies 1 --max-per-target 48 \
  --elements Se --wavelength 0.97918 --fpp "Se=3.841352735389343" --fcalc fc \
  --auto-fpp \
  --ccp4-root /Applications/ccp4-9 \
  --handdm-script /Volumes/Fast_4T/Workbuddy_Oasis/oasis-dm/src/handdm_v2.py \
  --handdm-hist-lib /Volumes/Fast_4T/Workbuddy_Oasis/oasis-dm/hist.lib \
  --minirefmac-root /Volumes/Fast_4T/Workbuddy_Oasis/oasis-refine \
  --sad "$INS/selected_sad.mtz" \
  --native "$INS/native.mtz" \
  --sequence "$INS/building.fasta" \
  --heavy "$INS/refined_heavy.pdb" \
  --output "$OUT"
