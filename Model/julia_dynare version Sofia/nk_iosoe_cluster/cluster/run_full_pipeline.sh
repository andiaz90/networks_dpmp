#!/bin/bash
#SBATCH --job-name=nkiosoe_full
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --output=nkiosoe_full_%j.out
#SBATCH --error=nkiosoe_full_%j.err
# ==========================================================================
# FULL pipeline in one job: data moments -> compile model -> SMM estimation
# -> shock plots (with the estimated parametrization).
#
# FIRST TIME ONLY (login node, before sbatch):
#   julia --project=. cluster/setup_cluster.jl
#
# SUBMIT:
#   sbatch cluster/run_full_pipeline.sh
# ==========================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$(dirname "${BASH_SOURCE[0]}")")" && pwd)"
cd "$ROOT"

echo "Job $SLURM_JOB_ID on $SLURMD_NODENAME  ($SLURM_CPUS_PER_TASK cpus)  started $(date)"

# Match to your cluster's Julia module name/version.
module load Julia/1.11.6 2>/dev/null || module load julia/1.11.6 2>/dev/null || \
module load julia 2>/dev/null || echo "WARNING: no Julia module — assuming julia in PATH"

# Local depot inside the package keeps the install self-contained + reusable.
export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
mkdir -p "$ROOT/.julia_depot"

export JULIA=julia
export NTHREADS=$SLURM_CPUS_PER_TASK

# Make sure packages are present (instant if setup_cluster.jl already ran).
julia --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'

bash "$ROOT/run_pipeline.sh" all

# Job-stamped copies of the key estimation outputs.
for f in smm_results smm_estimates smm_checkpoint; do
  [ -f "$ROOT/Data/${f}.csv" ] && cp "$ROOT/Data/${f}.csv" "$ROOT/Data/${f}_${SLURM_JOB_ID}.csv"
done

echo "Finished $(date)  (elapsed ${SECONDS}s)"
sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true
