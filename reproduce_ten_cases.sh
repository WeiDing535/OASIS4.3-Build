#!/usr/bin/env bash
# =============================================================================
# OASIS4.3 ten-case SAD phasing reproduction script (relative-path version)
# -----------------------------------------------------------------------------
# Usage:
#   bash reproduce_ten_cases.sh <case ID>      e.g.: bash reproduce_ten_cases.sh 3O2E
#   bash reproduce_ten_cases.sh all           serially reproduce all 10 cases
#
# Available case IDs: 3O2E 5D9D 5HHK 6HF7 6ZPP 7E1D 7ENM 7TBS 7TFQ 7TJ1
#
# Dependencies (OASIS does not bundle any CCP4 components; all are user-supplied with specified paths):
#   1) Python interpreter: default python3.12, overridable via the environment variable PHASEBUILD_PYTHON;
#      Requires Python 3.12+ with gemmi 0.7.5 installed.
#   2) User-supplied CCP4: set CCP4_ROOT (default /Applications/ccp4-9); the directory must contain
#      bin/phaser, bin/cbuccaneer, bin/refmac5, bin/cparrot and lib/data.
#   3) Input data: under the bundled validation_ten_cases/inputs/<PDB>/
#      selected_sad.mtz / native.mtz / building.fasta / refined_heavy.pdb
#      (7TFQ also has observations.mtz). These are 18 MB public example data, bundled with the package.
#
# Output: bundled work/reproduce/<case ID>/ (work/ is ignored by .gitignore and not committed).
#
# Parameter source: validation_ten_cases/verify_3fa5b86_20261116/cases/<PDB>/command.txt.
# Old command.txt contained absolute paths from the dev machine (/Volumes/Fast_4T/...) and some paths were
# truncated by spaces and became invalid; this script always uses $(dirname "$0") to locate the package root and uses paths relative to the package root throughout.
# =============================================================================
set -euo pipefail

# ---------- Package root location and global paths (all relative to package root, no absolute dev-machine paths) ----------
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY="${PHASEBUILD_PYTHON:-python3.12}"
RUNNER="$ROOT/run_sad_to_model.py"
CCP4_ROOT="${CCP4_ROOT:-/Applications/ccp4-9}"
DM_EXEC="$CCP4_ROOT/bin/cparrot"
INS="$ROOT/validation_ten_cases/inputs"
OUT="$ROOT/work/reproduce"

ALL_CASES=(3O2E 5D9D 5HHK 6HF7 6ZPP 7E1D 7ENM 7TBS 7TFQ 7TJ1)

echo "[Environment] Python: $PY (requires Python 3.12+ and gemmi 0.7.5)" >&2

# ---------- Pre-checks ----------
if ! command -v "$PY" >/dev/null 2>&1; then
    echo "[Error] Python interpreter not found: ${PY}. Please install Python 3.12+ and set PHASEBUILD_PYTHON." >&2
    exit 1
fi
if [ ! -f "$RUNNER" ]; then
    echo "[Error] Main program not found: ${RUNNER} (should be in the package root)." >&2
    exit 1
fi
if [ ! -x "$DM_EXEC" ] || [ ! -x "$CCP4_ROOT/bin/phaser" ]; then
    echo "[Error] CCP4 incomplete or wrong path: ${CCP4_ROOT} (missing bin/cparrot or bin/phaser)." >&2
    echo "       Please set CCP4_ROOT to point to the user-supplied CCP4 installation directory." >&2
    exit 1
fi
# ---------- Common parameters (identical for all 10 cases, from command.txt) ----------
COMMON=(
    --phase-cycles 100 --dm-cycles 40 --rounds 5 --cycles 10 --free-value 0
    --seed 20261023
    --aoe 0.5 --triplet-weight 0.1 --kmin 0.03 --max-neighbors 600
    --build-model --no-auto-hand --input-columns
    --building-recipe pipeline
    --no-canonical-dm-seed --no-per-round-cleanup
    --dm-engine parrot --dm-solvent 0.5 --dm-mask-radius 3.0 --dm-anisotropy
    --trust-kl 0.0005 --trust-cleanup --no-trust-stability --trust-wang-radius 8.0
    --feedback-cycles 2 --trace-resolution 2.0 --feedback-endpoints 2,3
    --auto-fpp
)

# ---------- Single-case execution ----------
run_case() {
    local pdb="$1"
    local idir="$INS/$pdb"
    local odir="$OUT/$pdb"
    # Default: refinement prior and trust cleanup switches for the 9 pplus_bundled cases in command.txt
    local REFPRIOR=dm_phifom
    local -a TRUST=(--trust-partial-cleanup --trust-work-weight)

    case "$pdb" in
        3O2E)
            # Special case: K edge, --fcalc fc, --dm-solvents auto,0.5
            COPIES=1 METHOD=pplus_bundled MAXPT=48 ELEMENTS=K WAVELENGTH=1.5418
            FPP="K=1.067139540598837" FCALC=fc EXTRA=(--dm-solvents auto,0.5)
            ;;
        5D9D)
            COPIES=1 METHOD=pplus_bundled MAXPT=48 ELEMENTS=Hg WAVELENGTH=0.9839
            FPP="Hg=9.743289557438965" FCALC=fc EXTRA=()
            ;;
        5HHK)
            # Special case: --max-per-target 13, copies 2
            COPIES=2 METHOD=pplus_bundled MAXPT=13 ELEMENTS=Se WAVELENGTH=0.97918
            FPP="Se=3.841352735389343" FCALC=fc EXTRA=()
            ;;
        6HF7)
            # Special case: --max-per-target 24, dual elements Tb,Mg, copies 3
            COPIES=3 METHOD=pplus_bundled MAXPT=24 ELEMENTS=Tb,Mg WAVELENGTH=1.6486
            FPP="Tb=10.595324959486323,Mg=0.20230138019255747" FCALC=fc EXTRA=()
            ;;
        6ZPP)
            COPIES=1 METHOD=pplus_bundled MAXPT=48 ELEMENTS=Se WAVELENGTH=0.97918
            FPP="Se=3.841352735389343" FCALC=fc EXTRA=()
            ;;
        7E1D)
            COPIES=2 METHOD=pplus_bundled MAXPT=48 ELEMENTS=Se WAVELENGTH=0.97852
            FPP="Se=3.8357741149160507" FCALC=fc EXTRA=()
            ;;
        7ENM)
            COPIES=1 METHOD=pplus_bundled MAXPT=48 ELEMENTS=Se WAVELENGTH=0.979
            FPP="Se=3.839830253513282" FCALC=fc EXTRA=()
            ;;
        7TBS)
            COPIES=1 METHOD=pplus_bundled MAXPT=48 ELEMENTS=Se WAVELENGTH=0.97937
            FPP="Se=3.8429606462715618" FCALC=fc EXTRA=()
            ;;
        7TFQ)
            # Special case: --method pplus_trust + observation_prior,
            #       refinement-prior=original_sad, disable trust-partial/work-weight
            COPIES=1 METHOD=pplus_trust MAXPT=48 ELEMENTS=Se WAVELENGTH=0.9796
            FPP="Se=0.4992114320528286" FCALC=fc
            REFPRIOR=original_sad
            TRUST=(--no-trust-partial-cleanup --no-trust-work-weight)
            EXTRA=(--observation-prior "$idir/observations.mtz")
            ;;
        7TJ1)
            COPIES=3 METHOD=pplus_bundled MAXPT=48 ELEMENTS=Se WAVELENGTH=0.9796
            FPP="Se=0.4992114320528286" FCALC=fc EXTRA=()
            ;;
        *)
            echo "[Error] Unknown case ID: ${pdb}. Available: ${ALL_CASES[*]}" >&2
            return 2
            ;;
    esac

    for f in selected_sad.mtz native.mtz building.fasta refined_heavy.pdb; do
        if [ ! -f "$idir/$f" ]; then
            echo "[Error] $pdb missing input file: $idir/$f" >&2
            return 1
        fi
    done

    mkdir -p "$odir"
    echo "===== Reproducing ${pdb}: method=$METHOD copies=$COPIES elements=$ELEMENTS wl=$WAVELENGTH max-per-target=$MAXPT =====" >&2

    local cmd=(
        "$PY" -u "$RUNNER"
        "${COMMON[@]}"
        --method "$METHOD"
        --copies "$COPIES"
        --max-per-target "$MAXPT"
        --elements "$ELEMENTS" --wavelength "$WAVELENGTH" --fpp "$FPP"
        --fcalc "$FCALC"
        --refinement-prior "$REFPRIOR"
        "${TRUST[@]+"${TRUST[@]}"}"
        "${EXTRA[@]+"${EXTRA[@]}"}"
        --dm-executable "$DM_EXEC"
        --ccp4-root "$CCP4_ROOT"
        --sad "$idir/selected_sad.mtz"
        --native "$idir/native.mtz"
        --sequence "$idir/building.fasta"
        --heavy "$idir/refined_heavy.pdb"
        --output "$odir"
    )
    "${cmd[@]}"
    echo "===== $pdb completed, output directory: $odir =====" >&2
}

# ---------- Entry point ----------
if [ "$#" -ne 1 ]; then
    echo "Usage: bash $(basename "$0") <case ID|all>" >&2
    echo "Case IDs: ${ALL_CASES[*]}" >&2
    exit 1
fi

if [ "$1" = "all" ]; then
    for pdb in "${ALL_CASES[@]}"; do
        run_case "$pdb"
    done
else
    run_case "$1"
fi
