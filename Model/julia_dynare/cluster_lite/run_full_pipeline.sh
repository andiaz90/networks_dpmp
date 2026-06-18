#!/bin/bash
#SBATCH --job-name=nkiosoe_full
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=2-00:00:00
#SBATCH --output=nkiosoe_full_%j.out
#SBATCH --error=nkiosoe_full_%j.err
# ==========================================================================
# FULL pipeline in one job: data moments -> compile model -> SMM estimation
# -> shock plots (with the estimated parametrization).
#
# FIRST TIME ONLY (login node):  julia --project=. cluster/setup_cluster.jl
# SUBMIT:                        sbatch cluster/run_full_pipeline.sh
# ==========================================================================
set -euo pipefail

echo "================================================================"
echo "  NK-IOSOE full pipeline"
echo "  Job ID:    ${SLURM_JOB_ID:-local}   Node: $(hostname)"
echo "  CPUs:      ${SLURM_CPUS_PER_TASK:-$(nproc)}   Mem: ${SLURM_MEM_PER_NODE:-?} MB"
echo "  Started:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================"

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)}"
[ -f "$ROOT/run_pipeline.sh" ] || ROOT="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)"
cd "$ROOT"
echo "Working directory: $(pwd)"

# Julia module: cascade -> EasyBuild fallback -> hard fail.
if ! command -v module &>/dev/null && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
fi
if command -v module &>/dev/null; then
    module load Julia/1.11.6 2>/dev/null || module load Julia/1.11 2>/dev/null || \
    module load julia/1.11.6 2>/dev/null || module load julia 2>/dev/null || true
fi
if ! command -v julia &>/dev/null; then
    JULIA_FALLBACK_BIN="/repositorio/modules/.local/easybuild/software/Julia/1.11.6-linux-x86_64/bin"
    [ -x "${JULIA_FALLBACK_BIN}/julia" ] && export PATH="${JULIA_FALLBACK_BIN}:$PATH"
fi
command -v julia &>/dev/null || { echo "ERROR: Julia not found."; exit 1; }
echo "Julia: $(julia --version)"

export NTHREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export JULIA_NUM_THREADS="$NTHREADS"
export JULIA_CPU_TARGET="generic"
export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export GKSwstype=100
export JULIA=julia
mkdir -p "$ROOT/.julia_depot"

# Wall-clock self-limit for the estimation stage (30 min under the allocation).
if [ -z "${SMM_MAX_HOURS:-}" ]; then
    LIM=$(squeue -h -j "${SLURM_JOB_ID:-0}" -o "%l" 2>/dev/null || echo "")
    case "$LIM" in
        *-*)   D=${LIM%%-*}; HMS=${LIM#*-}; H=${HMS%%:*}; export SMM_MAX_HOURS=$(awk "BEGIN{printf \"%.2f\", $D*24 + $H - 0.5}") ;;
        *:*:*) H=${LIM%%:*}; M=${LIM#*:}; M=${M%%:*};     export SMM_MAX_HOURS=$(awk "BEGIN{printf \"%.2f\", $H + $M/60 - 0.5}") ;;
        *)     export SMM_MAX_HOURS=47.5 ;;
    esac
fi
echo "Threads: $NTHREADS | SMM_MAX_HOURS=$SMM_MAX_HOURS"

julia --startup-file=no --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'

# Pre-flight smoke test before the long run.
echo ""; echo "-- smoke test --"
if ! julia --startup-file=no --project="$ROOT" "$ROOT/Model/julia_dynare/smoke_test.jl"; then
    echo "ABORTING: smoke test failed — pipeline NOT launched."; exit 1
fi

echo ""; echo "-- pipeline --"
set +e
stdbuf -oL -eL bash "$ROOT/run_pipeline.sh" all
JULIA_EXIT=$?
set -e

for f in smm_results smm_estimates smm_checkpoint; do
  [ -f "$ROOT/Data/${f}.csv" ] && cp "$ROOT/Data/${f}.csv" "$ROOT/Data/${f}_${SLURM_JOB_ID}.csv"
done
echo ""
echo "Exit code: $JULIA_EXIT   Finished: $(date)   Elapsed: ${SECONDS}s"
sacct -j "${SLURM_JOB_ID:-0}" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true
exit $JULIA_EXIT
