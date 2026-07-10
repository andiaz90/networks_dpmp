#!/bin/bash
#SBATCH --job-name=nkiosoe_smm
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=180G
#SBATCH --time=2-00:00:00
#SBATCH --output=nkiosoe_smm_%j.out
#SBATCH --error=nkiosoe_smm_%j.err
# ==========================================================================
# SMM estimation: data moments -> compile model -> CMA-ES estimation.
# Produces Model/julia_dynare/estimation_results/ (estimates, checkpoint,
# best_sol.txt, min_loss.txt, progress log). Copy that folder back to your
# machine (same relative path) and main_SOE_gap.jl uses it automatically:
#   scp -r <cluster>:.../nk_iosoe_cluster/Model/julia_dynare/estimation_results Model/julia_dynare/
# Run run_shocks.sh
# afterwards to generate the IRF plots with these estimates.
#
# FIRST TIME ONLY (login node):  julia --project=. cluster/setup_cluster.jl
# SUBMIT:                        sbatch cluster/run_estimation.sh
# OVERRIDE walltime/cpus:        sbatch --time=12:00:00 --cpus-per-task=16 cluster/run_estimation.sh
#
# MONITOR:
#   tail -f nkiosoe_smm_<jobid>.out        # live: eval count, best obj, fit decomposition
#                                          #   [Y= PH= L= Agg= Rk= CY= NL=], fail count,
#                                          #   ms/eval, Klein cache hit rate
#   tail -f Model/julia_dynare/estimation_results/smm_progress_log.csv   # trajectory, 1 row / 50 evals
#   cat Model/julia_dynare/estimation_results/min_loss.txt               # current best loss
#   cat Model/julia_dynare/estimation_results/best_sol.txt               # current best theta
#   squeue -u $USER ; sacct -j <jobid> --format=JobID,MaxRSS,Elapsed,State
#
# The estimation writes Model/julia_dynare/estimation_results/smm_checkpoint.csv (+
# best_sol.txt/min_loss.txt) LIVE on EVERY improvement, so a killed job loses
# almost nothing — resubmit and it warm-starts.
# NOTE 2026-07-08: θ is now 36 params (kappaw added) — an old 35-length
# checkpoint is auto-rejected and the run cold-starts from default θ₀.
# ==========================================================================
set -euo pipefail

# --- banner ---------------------------------------------------------------- #
echo "================================================================"
echo "  NK-IOSOE SMM estimation"
echo "  Job ID:    ${SLURM_JOB_ID:-local}"
echo "  Node:      $(hostname)"
echo "  CPUs:      ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "  Memory:    ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "  Started:   $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================================"

# --- root resolution (SLURM may copy the script to /var/spool/slurmd) ------ #
ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)}"
# If SLURM_SUBMIT_DIR isn't the bundle root, fall back to parent-of-cluster/.
if [ ! -f "$ROOT/run_pipeline.sh" ]; then
    ROOT="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)"
fi
cd "$ROOT"
echo "Working directory: $(pwd)"

# --- Julia module: cascade, then EasyBuild fallback, then hard-fail -------- #
if ! command -v module &>/dev/null && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
fi
if command -v module &>/dev/null; then
    module load Julia/1.11.6 2>/dev/null || module load Julia/1.11 2>/dev/null || \
    module load julia/1.11.6 2>/dev/null || module load julia/1.11 2>/dev/null || \
    module load julia 2>/dev/null || true
fi
if ! command -v julia &>/dev/null; then
    JULIA_FALLBACK_BIN="/repositorio/modules/.local/easybuild/software/Julia/1.11.6-linux-x86_64/bin"
    [ -x "${JULIA_FALLBACK_BIN}/julia" ] && export PATH="${JULIA_FALLBACK_BIN}:$PATH"
fi
if ! command -v julia &>/dev/null; then
    echo "ERROR: Julia not found (tried module names + ${JULIA_FALLBACK_BIN:-<none>})."
    exit 1
fi
echo "Julia: $(julia --version)"

# --- environment ----------------------------------------------------------- #
export NTHREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export JULIA_NUM_THREADS="$NTHREADS"
export JULIA_CPU_TARGET="generic"                       # heterogeneous nodes
export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export GKSwstype=100                                     # headless GR (no X11)
export JULIA=julia
mkdir -p "$ROOT/.julia_depot"

# Wall-clock self-limit: 30 min under the SLURM allocation so the optimiser
# stops itself and flushes the final checkpoint before SLURM SIGKILLs the job.
# Derive hours from the SLURM time limit when available; else default to 47.5h
# (matches the 2-day header above). Override by exporting SMM_MAX_HOURS yourself.
if [ -z "${SMM_MAX_HOURS:-}" ]; then
    LIM=$(squeue -h -j "${SLURM_JOB_ID:-0}" -o "%l" 2>/dev/null || echo "")
    case "$LIM" in
        *-*)   D=${LIM%%-*}; HMS=${LIM#*-}; H=${HMS%%:*}; export SMM_MAX_HOURS=$(awk "BEGIN{printf \"%.2f\", $D*24 + $H - 0.5}") ;;
        *:*:*) H=${LIM%%:*}; M=${LIM#*:}; M=${M%%:*};     export SMM_MAX_HOURS=$(awk "BEGIN{printf \"%.2f\", $H + $M/60 - 0.5}") ;;
        *)     export SMM_MAX_HOURS=47.5 ;;
    esac
fi
echo "Threads: $NTHREADS | depot: $JULIA_DEPOT_PATH | SMM_MAX_HOURS=$SMM_MAX_HOURS"

# --- instantiate (fast if setup_cluster.jl already ran) -------------------- #
julia --startup-file=no --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'

# --- pre-flight smoke test (abort before burning the allocation) ----------- #
echo ""
echo "-- smoke test --------------------------------------------------"
if ! julia --startup-file=no --project="$ROOT" "$ROOT/Model/julia_dynare/smoke_test.jl"; then
    echo "ABORTING: smoke test failed — estimation NOT launched. Fix and resubmit."
    exit 1
fi

# --- launch (line-buffered so tail -f streams; capture exit under set +e ---- #
# --- so the post-run archive always runs even on a nonzero Julia exit) ------ #
echo ""
echo "-- estimation --------------------------------------------------"
set +e
stdbuf -oL -eL bash "$ROOT/run_pipeline.sh" estimate
JULIA_EXIT=$?
set -e

# --- post-run: job-stamped copies + accounting ----------------------------- #
ESTDIR="$ROOT/Model/julia_dynare/estimation_results"
for f in smm_results.csv smm_estimates.csv smm_checkpoint.csv best_sol.txt min_loss.txt; do
    base="${f%.*}"; ext="${f##*.}"
    [ -f "$ESTDIR/$f" ] && cp "$ESTDIR/$f" "$ESTDIR/${base}_${SLURM_JOB_ID}.${ext}"
done
if [ -f "$ESTDIR/min_loss.txt" ]; then
    echo ""
    echo "  Best loss: $(head -1 "$ESTDIR/min_loss.txt")   (see $ESTDIR/best_sol.txt)"
    echo "  Copy results back:  scp -r <cluster>:$ESTDIR  Model/julia_dynare/"
fi
echo ""
echo "================================================================"
echo "  Exit code: $JULIA_EXIT   Finished: $(date)   Elapsed: ${SECONDS}s"
echo "================================================================"
sacct -j "${SLURM_JOB_ID:-0}" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true

exit $JULIA_EXIT
