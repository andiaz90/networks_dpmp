#!/bin/bash
# ==========================================================================
# prepare_cluster_package.sh
#
# Run this LOCALLY (on your Mac/laptop) before uploading to the cluster.
# It bundles all files the cluster needs into a single tarball.
#
# USAGE:
#   cd /path/to/Networks-DPMP/Model/julia_dynare
#   bash cluster/prepare_cluster_package.sh
#
# OUTPUT:
#   nk_iosoe_cluster_YYYYMMDD.tar.gz  (in the julia_dynare directory)
#
# WHAT IT CHECKS:
#   1. nk_iosoe_context.jls exists (requires main_SOE_gap.jl to have run)
#   2. sectoral_moments.csv and aggregate_moments.csv exist
#      (requires compute_data_moments.jl to have run)
#   3. Compiled Jacobian files exist in mod/NK_SOE_lev_gap2/model/julia/
# ==========================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JULIA_DIR="$(dirname "$SCRIPT_DIR")"
REPO_ROOT="$(dirname "$(dirname "$JULIA_DIR")")"
DATA_DIR="$REPO_ROOT/Data"
DATE=$(date +%Y%m%d)
BUNDLE="nk_iosoe_cluster_${DATE}.tar.gz"

echo "======================================================"
echo "  NK-IOSOE Cluster Package Preparation"
echo "======================================================"
echo "Julia dir:  $JULIA_DIR"
echo "Data dir:   $DATA_DIR"
echo "Output:     $JULIA_DIR/$BUNDLE"
echo ""

# --------------------------------------------------------------------------
# 1. Check prerequisites
# --------------------------------------------------------------------------
echo "--- Checking prerequisites ---"

CONTEXT="$JULIA_DIR/mod/nk_iosoe_context.jls"
ENDO_NAMES="$JULIA_DIR/mod/dynare_endo_names.csv"
JACOBIAN="$JULIA_DIR/mod/NK_SOE_lev_gap2/model/julia/SparseDynamicG1!.jl"
SEC_MOM="$DATA_DIR/sectoral_moments.csv"
AGG_MOM="$DATA_DIR/aggregate_moments.csv"

MISSING=0
for f in "$CONTEXT" "$ENDO_NAMES" "$JACOBIAN" "$SEC_MOM" "$AGG_MOM"; do
    if [ -f "$f" ]; then
        SIZE=$(du -sh "$f" | cut -f1)
        echo "  OK  $SIZE  $(basename $f)"
    else
        echo "  MISSING: $f"
        MISSING=1
    fi
done

if [ "$MISSING" -eq 1 ]; then
    echo ""
    echo "ERROR: Missing required files. Before running this script:"
    echo "  1. Run main_SOE_gap.jl to generate nk_iosoe_context.jls"
    echo "  2. Run compute_data_moments.jl to generate sectoral/aggregate_moments.csv"
    exit 1
fi
echo ""

# --------------------------------------------------------------------------
# 2. Create staging directory
# --------------------------------------------------------------------------
STAGE="$JULIA_DIR/_cluster_stage"
rm -rf "$STAGE"
mkdir -p "$STAGE/julia_code/mod"
mkdir -p "$STAGE/Data"
mkdir -p "$STAGE/cluster"

echo "--- Copying Julia source files ---"
# Core estimation files
for f in \
    run_smm_estimation.jl \
    smm_estimation.jl \
    smm_model_moments.jl \
    steady_ntwsoe.jl \
    steady_ntwsoe_system.jl \
    utils.jl
do
    cp "$JULIA_DIR/$f" "$STAGE/julia_code/"
    echo "  $f"
done

echo ""
echo "--- Copying compiled Dynare model ---"
# NK_SOE_lev_gap2 compiled files (needed for Klein solver + Jacobians)
mkdir -p "$STAGE/julia_code/mod/NK_SOE_lev_gap2/model/julia"
mkdir -p "$STAGE/julia_code/mod/NK_SOE_lev_gap2/model/json"
cp "$JULIA_DIR/mod/NK_SOE_lev_gap2/model/julia/"*.jl \
   "$STAGE/julia_code/mod/NK_SOE_lev_gap2/model/julia/"
cp "$JULIA_DIR/mod/NK_SOE_lev_gap2/model/json/"*.json \
   "$STAGE/julia_code/mod/NK_SOE_lev_gap2/model/json/"
echo "  NK_SOE_lev_gap2 model/julia/ ($(ls "$STAGE/julia_code/mod/NK_SOE_lev_gap2/model/julia/" | wc -l | tr -d ' ') files)"

# .mod source files and generated parameter file
for f in \
    NK_SOE_lev_gap2.mod \
    NK_SOE_lev_gap2_smm.mod \
    definition_block_io.mod \
    definition_block_lab.mod \
    definition_block_nsec.mod \
    solution_block.mod \
    params_jl.mod \
    dynare_endo_names.csv
do
    [ -f "$JULIA_DIR/mod/$f" ] && cp "$JULIA_DIR/mod/$f" "$STAGE/julia_code/mod/" && echo "  $f"
done

echo ""
echo "--- Copying Dynare context ---"
# Serialized Dynare context (contains compiled model + steady state)
cp "$JULIA_DIR/mod/nk_iosoe_context.jls" "$STAGE/julia_code/mod/"
SIZE=$(du -sh "$JULIA_DIR/mod/nk_iosoe_context.jls" | cut -f1)
echo "  nk_iosoe_context.jls  ($SIZE)"

# Copy smm context too if it exists
[ -f "$JULIA_DIR/mod/nk_iosoe_smm_context.jls" ] && \
    cp "$JULIA_DIR/mod/nk_iosoe_smm_context.jls" "$STAGE/julia_code/mod/" && \
    echo "  nk_iosoe_smm_context.jls"

echo ""
echo "--- Copying data moments ---"
cp "$DATA_DIR/sectoral_moments.csv"  "$STAGE/Data/"
cp "$DATA_DIR/aggregate_moments.csv" "$STAGE/Data/"
echo "  sectoral_moments.csv"
echo "  aggregate_moments.csv"

# Copy any existing checkpoint (warm start for the cluster)
[ -f "$DATA_DIR/smm_checkpoint.csv" ] && \
    cp "$DATA_DIR/smm_checkpoint.csv" "$STAGE/Data/" && \
    echo "  smm_checkpoint.csv (warm start)"

echo ""
echo "--- Copying cluster scripts ---"
cp "$SCRIPT_DIR/run_smm_hpc.sh"      "$STAGE/cluster/"
cp "$SCRIPT_DIR/setup_cluster.jl"    "$STAGE/cluster/"
cp "$SCRIPT_DIR/Project.toml"        "$STAGE/cluster/"
[ -f "$SCRIPT_DIR/Manifest.toml" ] && cp "$SCRIPT_DIR/Manifest.toml" "$STAGE/cluster/"
echo "  run_smm_hpc.sh"
echo "  setup_cluster.jl"
echo "  Project.toml"

# --------------------------------------------------------------------------
# 3. Write a cluster-side layout README
# --------------------------------------------------------------------------
cat > "$STAGE/README_LAYOUT.txt" << 'EOF'
NK-IOSOE Cluster Package Layout
================================

Extract with:
    tar -xzf nk_iosoe_cluster_YYYYMMDD.tar.gz

Expected directory tree after extraction:

    nk_iosoe_cluster/
    ├── julia_code/          ← Julia source files
    │   ├── run_smm_estimation.jl
    │   ├── smm_estimation.jl
    │   ├── smm_model_moments.jl
    │   ├── utils.jl
    │   ├── steady_ntwsoe.jl
    │   ├── steady_ntwsoe_system.jl
    │   └── mod/
    │       ├── nk_iosoe_context.jls        ← compiled model (15 MB)
    │       ├── dynare_endo_names.csv
    │       ├── params_jl.mod
    │       ├── NK_SOE_lev_gap2.mod
    │       └── NK_SOE_lev_gap2/model/julia/  ← compiled Jacobians
    ├── Data/
    │   ├── sectoral_moments.csv
    │   └── aggregate_moments.csv
    └── cluster/
        ├── run_smm_hpc.sh       ← SLURM submission script
        ├── setup_cluster.jl     ← one-time package install
        └── Project.toml         ← Julia package spec

STEPS ON CLUSTER:
  1. cd nk_iosoe_cluster
  2. julia --project=cluster cluster/setup_cluster.jl    # first time only
  3. sbatch cluster/run_smm_hpc.sh

OUTPUT FILES (written to Data/):
  smm_results.csv      — full moment-fit table
  smm_estimates.csv    — estimated parameters
  smm_checkpoint.csv   — warm-start checkpoint
  smm_results_JOBID.csv, smm_estimates_JOBID.csv  — job-stamped copies
EOF

# --------------------------------------------------------------------------
# 4. Build tarball
# --------------------------------------------------------------------------
echo ""
echo "--- Building tarball ---"
cd "$JULIA_DIR"
STAGE_NAME="nk_iosoe_cluster"
mv "$STAGE" "$STAGE_NAME"
tar -czf "$BUNDLE" "$STAGE_NAME"
rm -rf "$STAGE_NAME"

BUNDLE_SIZE=$(du -sh "$BUNDLE" | cut -f1)
echo ""
echo "======================================================"
echo "  Bundle ready: $BUNDLE  ($BUNDLE_SIZE)"
echo "======================================================"
echo ""
echo "Upload to cluster:"
echo "  scp $JULIA_DIR/$BUNDLE user@cluster.host:~/"
echo ""
echo "On cluster:"
echo "  tar -xzf $BUNDLE"
echo "  cd nk_iosoe_cluster"
echo "  julia --project=cluster cluster/setup_cluster.jl   # first time only"
echo "  sbatch cluster/run_smm_hpc.sh"
echo ""
