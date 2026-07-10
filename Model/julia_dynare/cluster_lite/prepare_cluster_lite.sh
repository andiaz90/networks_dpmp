#!/bin/bash
# ==========================================================================
# prepare_cluster_lite.sh
#
# Build the LIGHTWEIGHT cluster bundle (estimation + shock plots). Unlike the
# old prepare_cluster_package.sh, this ships NO heavy compiled artifacts
# (no nk_iosoe_context.jls, no Jacobians, no precomputed moments): those are
# generated ON THE CLUSTER. It ships only source + the raw data needed to run.
#
# RUN LOCALLY from the repo:
#   cd /path/to/Networks-DPMP/Model/julia_dynare
#   bash cluster/prepare_cluster_lite.sh
#
# OUTPUT (next to this script's repo root): nk_iosoe_cluster_lite_YYYYMMDD.tar.gz
# ==========================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # cluster/
JD="$(dirname "$SCRIPT_DIR")"                                # Model/julia_dynare
REPO_ROOT="$(dirname "$(dirname "$JD")")"                    # Networks-DPMP
DATA="$REPO_ROOT/Data"
DATE=$(date +%Y%m%d)
STAGE="$JD/nk_iosoe_cluster"
BUNDLE="$JD/nk_iosoe_cluster_lite_${DATE}.tar.gz"

rm -rf "$STAGE"
mkdir -p "$STAGE/Model/julia_dynare/mod" "$STAGE/Data/computed" "$STAGE/cluster"

# --- Julia source (no compiled artifacts) ---
JL_FILES=(compute_data_moments.jl main_SOE_gap.jl run_smm_estimation.jl
  smm_estimation.jl smm_model_moments.jl smm_inference.jl steady_ntwsoe.jl steady_ntwsoe_system.jl
  utils.jl smoke_test.jl figs_SOE_gap.jl plot_scripts.jl run_dynare_subprocess.jl
  run_all_shocks.jl shock_plots_common.jl
  oil_shock_analysis.jl agr_shock_analysis.jl agrmin_shock_analysis.jl
  min_shock_analysis.jl mfg_shock_analysis.jl plot_chiib_comparison.jl)
for f in "${JL_FILES[@]}"; do cp "$JD/$f" "$STAGE/Model/julia_dynare/"; done

# --- mod/ source only ---
MOD_FILES=(NK_SOE_lev_gap2.mod definition_block_nsec.mod
  definition_block_io.mod definition_block_lab.mod solution_block.mod run_dynare_model.jl)
for f in "${MOD_FILES[@]}"; do cp "$JD/mod/$f" "$STAGE/Model/julia_dynare/mod/"; done

# --- raw data needed to run (incl. fpa); moments are recomputed on cluster ---
DATA_FILES=(sector_calibration.csv IO_2021_chile.csv fpa_vector_few_industries_chile.csv
  count_workers_by_sector.csv deflactor_pib.csv pib_sectorial_bc.csv
  reer_chile_bis.xlsx datos_CCNN_mayo2025.xlsx)
for f in "${DATA_FILES[@]}"; do cp "$DATA/$f" "$STAGE/Data/"; done
cp "$DATA/computed/export_shares_chile.csv" "$STAGE/Data/computed/"

# --- runtime/config files (must sit beside this script in cluster/) ---
cp "$SCRIPT_DIR/Project.toml"          "$STAGE/Project.toml"
cp "$SCRIPT_DIR/setup_cluster.jl"      "$STAGE/cluster/"
cp "$SCRIPT_DIR/setup_cluster.sh"      "$STAGE/cluster/"
cp "$SCRIPT_DIR/run_full_pipeline.sh"  "$STAGE/cluster/"
cp "$SCRIPT_DIR/run_estimation.sh"     "$STAGE/cluster/"
cp "$SCRIPT_DIR/run_shocks.sh"         "$STAGE/cluster/"
cp "$SCRIPT_DIR/run_pipeline.sh"       "$STAGE/run_pipeline.sh"
cp "$SCRIPT_DIR/README_CLUSTER.md"     "$STAGE/README_CLUSTER.md"
chmod +x "$STAGE"/run_pipeline.sh "$STAGE"/cluster/*.sh

cd "$JD"
tar -czf "$BUNDLE" -C "$JD" "nk_iosoe_cluster"
rm -rf "$STAGE"
echo "Bundle: $BUNDLE  ($(du -sh "$BUNDLE" | cut -f1))"
