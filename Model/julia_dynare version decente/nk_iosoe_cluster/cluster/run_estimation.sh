#!/bin/bash
#SBATCH --job-name=nkiosoe_smm
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=1-12:00:00
#SBATCH --output=nkiosoe_smm_%j.out
#SBATCH --error=nkiosoe_smm_%j.err
# ==========================================================================
# SMM estimation only: data moments -> compile model -> CMA-ES estimation.
# Produces Data/smm_estimates.csv (the estimated theta). Run the shocks job
# afterwards to generate the IRF plots with these estimates.
#
# FIRST TIME ONLY (login node): julia --project=. cluster/setup_cluster.jl
# SUBMIT:                       sbatch cluster/run_estimation.sh
# ==========================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)"
cd "$ROOT"
echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME  ($SLURM_CPUS_PER_TASK cpus)  started $(date)"

module load Julia/1.11.6 2>/dev/null || module load julia/1.11.6 2>/dev/null || \
module load julia 2>/dev/null || echo "WARNING: no Julia module — assuming julia in PATH"

export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
mkdir -p "$ROOT/.julia_depot"
export JULIA=julia
export NTHREADS=$SLURM_CPUS_PER_TASK

julia --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'

bash "$ROOT/run_pipeline.sh" estimate

for f in smm_results smm_estimates smm_checkpoint; do
  [ -f "$ROOT/Data/${f}.csv" ] && cp "$ROOT/Data/${f}.csv" "$ROOT/Data/${f}_${SLURM_JOB_ID}.csv"
done
echo "Finished $(date)  (elapsed ${SECONDS}s)"
sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true
