#!/bin/bash
#SBATCH --job-name=nk_kappa
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=20:00:00
#SBATCH --output=kappa_%j.out
#SBATCH --error=kappa_%j.err

# ==========================================================================
# SMM estimation of the SECTORAL PRICE RIGIDITIES kappa_i (theta[39:51]).
#
#   sbatch cluster/run_kappa_smm.sh affine          # headline: level + dispersion
#   sbatch cluster/run_kappa_smm.sh free 2.0        # eleven free kappa_i, ridge 2.0
#   sbatch cluster/run_kappa_smm.sh scale           # level only
#
# WHY --cpus-per-task=4 AND --threads=1, NOT 48.
#
# This job runs Julia SINGLE-THREADED on purpose. Two independent reasons, both
# already paid for in this project:
#
#   1. Under multithreading roughly 36% of evaluations return a STALE decision
#      rule from the per-thread Klein cache (networks-dpmp-klein-cache-stale-R).
#      A stale rule is a silently wrong objective, and for a parameter whose
#      whole question is "is it identified", a run polluted by stale evaluations
#      cannot answer it. --threads=1 is the standing recommendation for any run
#      whose result will be believed.
#   2. kappa_i is written through set_param!, which writes into the PER-THREAD
#      _SMM_PARAMS vector. That is thread-safe for the per-evaluation path, but
#      it is the same machinery that made SMM_KAPPA_SCALE require threads=1, and
#      there is no reason to relitigate it here.
#
# The 4 CPUs are for BLAS headroom and the OS, not for Julia threads. This job
# is therefore SLOW per evaluation but correct — budget wall time, not cores.
# If it needs to be faster, run several rungs as separate array tasks rather
# than threading one.
#
# FIRST TIME ONLY, on the login node:
#   julia --project=cluster cluster/setup_cluster.jl
# ==========================================================================

set -euo pipefail

MODE="${1:-affine}"
SHRINK="${2:-0.0}"

case "$MODE" in
    off|scale|affine|free) ;;
    *) echo "ERROR: mode must be one of off|scale|affine|free (got '$MODE')"; exit 1 ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JULIA_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$JULIA_DIR/../../Data"
cd "$JULIA_DIR"

echo "======================================================"
echo "  NK-IOSOE — sectoral price rigidity kappa_i"
echo "======================================================"
echo "Job ID:        ${SLURM_JOB_ID:-local}"
echo "Node:          ${SLURMD_NODENAME:-$(hostname)}"
echo "Mode:          SMM_KAPPA_MODE=$MODE"
[ "$MODE" = "free" ] && echo "Ridge:         SMM_KAPPA_SHRINK=$SHRINK"
echo "Started:       $(date)"
echo "======================================================"

module load Julia/1.11.6 2>/dev/null || \
module load julia/1.11.6 2>/dev/null || \
module load julia 2>/dev/null || \
echo "WARNING: no Julia module found — assuming julia is in PATH"
echo "Julia version: $(julia --version)"

export JULIA_NUM_THREADS=1
export JULIA_CPU_TARGET="generic"
export JULIA_DEPOT_PATH="$SCRIPT_DIR/.julia_depot:$HOME/.julia"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1
mkdir -p "$SCRIPT_DIR/.julia_depot"

# --------------------------------------------------------------------------
# Prerequisites. The moment CSVs MUST carry the price-rigidity block (85-109);
# smm_estimation.jl hard-errors otherwise, but failing here costs seconds
# instead of a queue slot.
# --------------------------------------------------------------------------
for f in "$JULIA_DIR/mod/nk_iosoe_context.jls" \
         "$JULIA_DIR/mod/dynare_endo_names.csv" \
         "$DATA_DIR/sectoral_moments.csv" \
         "$DATA_DIR/aggregate_moments.csv"; do
    [ -f "$f" ] || { echo "ERROR: required file not found: $f"; exit 1; }
done
head -1 "$DATA_DIR/sectoral_moments.csv" | grep -q "ratio_dPH_dY" || {
    echo "ERROR: sectoral_moments.csv predates the price-rigidity block (2026-09-02)."
    echo "       Regenerate locally:  bash run_kappa_ladder.sh 0"
    exit 1
}
grep -q "^rbar_dPH," "$DATA_DIR/aggregate_moments.csv" || {
    echo "ERROR: aggregate_moments.csv is missing rbar_dPH. Run: bash run_kappa_ladder.sh 0"
    exit 1
}
echo "Required files: OK (price-rigidity moment block present)"

julia --project="$SCRIPT_DIR" -e 'import Pkg; Pkg.instantiate(); println("  Packages: OK")'

# --------------------------------------------------------------------------
# The run. SMM_WARM_START=1 picks up the free dimensions from the checkpoint —
# CMA-ES lands in different basins from the same seed, so a cold start can
# return a WORSE point than one already on disk.
# --------------------------------------------------------------------------
echo ""
echo "=== SMM: kappa mode '$MODE', single-threaded ==="
echo ""

export SMM_KAPPA_MODE="$MODE"
[ "$MODE" = "free" ] && export SMM_KAPPA_SHRINK="$SHRINK"
export SMM_WARM_START="${SMM_WARM_START:-1}"

julia --project="$SCRIPT_DIR" --threads=1 --heap-size-hint=24G \
    "$JULIA_DIR/run_smm_estimation.jl"
JULIA_EXIT=$?

# --------------------------------------------------------------------------
# Stamp the results with the mode AND the job id. Without the mode in the name
# an affine run and a free run overwrite each other and the comparison table
# silently reads whichever finished last.
# --------------------------------------------------------------------------
TAG="kappa_${MODE}"
[ "$MODE" = "free" ] && TAG="${TAG}_sh${SHRINK}"
for f in smm_results smm_estimates smm_checkpoint; do
    src="$JULIA_DIR/estimation_results/${f}.csv"
    [ -f "$src" ] && cp "$src" "$JULIA_DIR/estimation_results/${f}_${TAG}.csv" && \
        echo "Saved: ${f}_${TAG}.csv"
done

# The comparison exhibit, built from whatever is now on disk.
julia --project="$SCRIPT_DIR" --threads=1 "$JULIA_DIR/compare_kappa_micro.jl" || \
    echo "(comparison script failed — run it by hand)"

echo ""
echo "======================================================"
echo "  Job Completed"
echo "======================================================"
echo "Exit code:   $JULIA_EXIT"
echo "Finished:    $(date)"
echo "Runtime:     $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m $((SECONDS % 60))s"
sacct -j "${SLURM_JOB_ID:-0}" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true
