#!/bin/bash
# ==========================================================================
# run_pipeline.sh  —  NK-IOSOE end-to-end driver (estimation + shock plots)
#
# This is the single source of truth for the run sequence. The SLURM scripts
# in cluster/ just set #SBATCH headers and call this with a stage argument.
# It can also be run directly on an interactive node / your laptop.
#
# USAGE:
#   ./run_pipeline.sh all          # moments -> compile -> estimate -> shocks
#   ./run_pipeline.sh estimate     # moments -> compile -> estimate
#   ./run_pipeline.sh shocks       # shocks only (needs Data/smm_estimates.csv)
#   ./run_pipeline.sh moments|compile     # individual stages
#
# ENV OVERRIDES:
#   JULIA=/path/to/julia     (default: julia in PATH)
#   NTHREADS=48              (default: nproc)
#
# WHY THIS ORDER:
#   1. moments   compute_data_moments.jl   -> Data/{sectoral,aggregate}_moments.csv
#   2. compile   main_SOE_gap.jl           -> mod/nk_iosoe_context.jls + Jacobians
#                                             (the heavy artifacts, built HERE not shipped)
#   3. estimate  run_smm_estimation.jl     -> Data/smm_estimates.csv  (the estimated theta)
#   4. shocks    run_all_shocks.jl         -> figures/tables, AUTO-loading smm_estimates.csv
#                                             so every IRF uses the ESTIMATED parametrization
# ==========================================================================
set -euo pipefail

STAGE="${1:-all}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JD="$ROOT/Model/julia_dynare"
JULIA="${JULIA:-julia}"
NTHREADS="${NTHREADS:-$(nproc 2>/dev/null || echo 1)}"

# Threading: Julia threads for the estimation; BLAS pinned to 1 to avoid contention.
export JULIA_NUM_THREADS="$NTHREADS"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

# Collect all shock figures/tables into results/ (instead of the local Dropbox/
# Overleaf path the scripts default to, which does not exist on the cluster).
export OVERLEAF_ROOT="${OVERLEAF_ROOT:-$ROOT/results}"
mkdir -p "$OVERLEAF_ROOT"

echo "======================================================"
echo "  NK-IOSOE pipeline   stage=$STAGE"
echo "  root      : $ROOT"
echo "  julia     : $($JULIA --version 2>/dev/null || echo "$JULIA")"
echo "  threads   : $NTHREADS"
echo "  results   : $OVERLEAF_ROOT"
echo "======================================================"

_run() { echo; echo ">>> $*"; "$@"; }

stage_moments() {
  _run "$JULIA" --project="$ROOT" "$JD/compute_data_moments.jl"
}
stage_compile() {
  _run "$JULIA" --project="$ROOT" "$JD/main_SOE_gap.jl"
}
stage_estimate() {
  _run "$JULIA" --project="$ROOT" --threads="$NTHREADS" --heap-size-hint=100G \
       "$JD/run_smm_estimation.jl"
}
stage_shocks() {
  # run_all_shocks.jl auto-loads Data/smm_estimates.csv inside each shock script.
  if [ ! -f "$ROOT/Data/smm_estimates.csv" ]; then
    echo "WARNING: Data/smm_estimates.csv not found — shocks will fall back to"
    echo "         calibrated defaults. Run the 'estimate' stage first for the"
    echo "         estimated parametrization."
  fi
  [ -f "$ROOT/Data/sectoral_moments.csv" ] || stage_moments
  _run "$JULIA" --project="$ROOT" "$JD/run_all_shocks.jl"
}

case "$STAGE" in
  moments)  stage_moments ;;
  compile)  stage_compile ;;
  estimate) stage_moments; stage_compile; stage_estimate ;;
  shocks)   stage_shocks ;;
  all)      stage_moments; stage_compile; stage_estimate; stage_shocks ;;
  *) echo "Unknown stage '$STAGE'. Use: all | estimate | shocks | moments | compile"; exit 2 ;;
esac

echo
echo "======================================================"
echo "  Stage '$STAGE' finished OK."
echo "  Estimated parameters : Data/smm_estimates.csv"
echo "  Moment fit           : Data/smm_results.csv"
echo "  Shock figures/tables : $OVERLEAF_ROOT  (+ Model/julia_dynare/figures, .../tables)"
echo "======================================================"
