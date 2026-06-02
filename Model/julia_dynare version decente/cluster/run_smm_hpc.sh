#!/bin/bash
#SBATCH --job-name=nk_iosoe_smm
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=14:00:00
#SBATCH --output=smm_%j.out
#SBATCH --error=smm_%j.err

# ==========================================================================
# SLURM script: SMM estimation for NK-IOSOE Chile model (Julia/Dynare.jl)
#
# USAGE:
#   sbatch run_smm_hpc.sh
#   sbatch --cpus-per-task=48 run_smm_hpc.sh   # override CPU count
#
# FIRST TIME ONLY — run setup on the login node before submitting:
#   julia --project="$SCRIPT_DIR" cluster/setup_cluster.jl
# ==========================================================================

set -euo pipefail

echo "======================================================"
echo "  NK-IOSOE SMM Estimation — HPC Job"
echo "======================================================"
echo "Job ID:       $SLURM_JOB_ID"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Memory:       $SLURM_MEM_PER_NODE MB"
echo "Node:         $SLURMD_NODENAME"
echo "Started:      $(date)"
echo "======================================================"

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# SCRIPT_DIR is the cluster/ subdirectory; the Julia code is one level up
JULIA_DIR="$(dirname "$SCRIPT_DIR")"
cd "$JULIA_DIR"
echo "Working directory: $(pwd)"

# --------------------------------------------------------------------------
# Julia module — adjust version to what your cluster provides
# --------------------------------------------------------------------------
module load Julia/1.11.6 2>/dev/null || \
module load julia/1.11.6 2>/dev/null || \
module load julia 2>/dev/null || \
echo "WARNING: no Julia module found — assuming julia is in PATH"

echo "Julia version: $(julia --version)"

# --------------------------------------------------------------------------
# Thread and BLAS settings
# NOTE: We use Julia Threads (not Distributed), so --threads=N, NOT -p N.
#       BLAS threads must be 1 to avoid contention with Julia threads.
# --------------------------------------------------------------------------
export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK
export JULIA_CPU_TARGET="generic"
export JULIA_DEPOT_PATH="$SCRIPT_DIR/.julia_depot:$HOME/.julia"

# Prevent BLAS/LAPACK from spawning competing threads
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export BLIS_NUM_THREADS=1

echo "Julia threads: $JULIA_NUM_THREADS"
echo "Julia depot:   $JULIA_DEPOT_PATH"
mkdir -p "$SCRIPT_DIR/.julia_depot"

# --------------------------------------------------------------------------
# Verify required files exist before launching
# --------------------------------------------------------------------------
CONTEXT_FILE="$JULIA_DIR/mod/nk_iosoe_context.jls"
ENDO_NAMES="$JULIA_DIR/mod/dynare_endo_names.csv"
SEC_MOMENTS="$JULIA_DIR/../../../Data/sectoral_moments.csv"
AGG_MOMENTS="$JULIA_DIR/../../../Data/aggregate_moments.csv"

for f in "$CONTEXT_FILE" "$ENDO_NAMES" "$SEC_MOMENTS" "$AGG_MOMENTS"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: required file not found: $f"
        echo "Run prepare_cluster_package.sh locally first."
        exit 1
    fi
done
echo "Required files: OK"

# --------------------------------------------------------------------------
# Instantiate packages (fast if already installed, slow only first time)
# --------------------------------------------------------------------------
echo "Checking Julia packages..."
julia --project="$SCRIPT_DIR" -e "
import Pkg
Pkg.instantiate()
println(\"  Packages: OK\")
"

# --------------------------------------------------------------------------
# Launch SMM estimation
# --------------------------------------------------------------------------
echo ""
echo "=== Starting SMM estimation with $JULIA_NUM_THREADS threads ==="
echo ""

julia \
    --project="$SCRIPT_DIR" \
    --threads="$JULIA_NUM_THREADS" \
    --heap-size-hint=100G \
    "$JULIA_DIR/run_smm_estimation.jl"

JULIA_EXIT=$?

# --------------------------------------------------------------------------
# Copy results with job ID stamp
# --------------------------------------------------------------------------
DATA_DIR="$JULIA_DIR/../../../Data"
for f in smm_results smm_estimates smm_checkpoint; do
    src="$DATA_DIR/${f}.csv"
    if [ -f "$src" ]; then
        cp "$src" "$DATA_DIR/${f}_${SLURM_JOB_ID}.csv"
        echo "Saved: ${f}_${SLURM_JOB_ID}.csv"
    fi
done

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo ""
echo "======================================================"
echo "  Job Completed"
echo "======================================================"
echo "Exit code:   $JULIA_EXIT"
echo "Finished:    $(date)"
echo "Runtime:     $((SECONDS / 3600))h $((SECONDS % 3600 / 60))m $((SECONDS % 60))s"
echo ""

sacct -j "$SLURM_JOB_ID" --format=JobID,JobName,MaxRSS,Elapsed,State 2>/dev/null || true

exit $JULIA_EXIT
