#!/bin/bash
#SBATCH --job-name=nkiosoe_shocks
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=nkiosoe_shocks_%j.out
#SBATCH --error=nkiosoe_shocks_%j.err
# ==========================================================================
# Shock plots only. Recompiles the model per shock and AUTO-loads
# Model/julia_dynare/estimation_results/smm_estimates.csv, so every IRF uses the
# ESTIMATED parametrization. Run AFTER run_estimation.sh.
#
# SUBMIT:  sbatch cluster/run_shocks.sh
# Subset:  edit the run_pipeline.sh call below, or run on a login node:
#          ./run_pipeline.sh shocks            (all four)
#          julia --project=. Model/julia_dynare/run_all_shocks.jl oil mfg
# ==========================================================================
set -euo pipefail

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)}"
[ -f "$ROOT/run_pipeline.sh" ] || ROOT="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)"
cd "$ROOT"
echo "Job ${SLURM_JOB_ID:-local} on $(hostname)  started $(date)"

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

export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
mkdir -p "$ROOT/.julia_depot"
export JULIA=julia
export NTHREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export JULIA_CPU_TARGET="generic"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 BLIS_NUM_THREADS=1
export GKSwstype=100   # headless GR (no X11) — required for plotting on compute nodes

julia --startup-file=no --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'

bash "$ROOT/run_pipeline.sh" shocks

echo "Finished $(date)  (elapsed ${SECONDS}s)"
sacct -j "${SLURM_JOB_ID:-0}" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true
