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
# Data/smm_estimates.csv, so every IRF uses the ESTIMATED parametrization.
# Run AFTER run_estimation.sh (needs Data/smm_estimates.csv).
#
# SUBMIT:  sbatch cluster/run_shocks.sh
# Subset:  edit the run_pipeline.sh call below, or run on a login node:
#          ./run_pipeline.sh shocks            (all four)
#          julia --project=. Model/julia_dynare/run_all_shocks.jl oil mfg
# ==========================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)"
cd "$ROOT"
echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME  started $(date)"

module load Julia/1.11.6 2>/dev/null || module load julia/1.11.6 2>/dev/null || \
module load julia 2>/dev/null || echo "WARNING: no Julia module — assuming julia in PATH"

export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
mkdir -p "$ROOT/.julia_depot"
export JULIA=julia
export NTHREADS=$SLURM_CPUS_PER_TASK

julia --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'

bash "$ROOT/run_pipeline.sh" shocks

echo "Finished $(date)  (elapsed ${SECONDS}s)"
sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true
