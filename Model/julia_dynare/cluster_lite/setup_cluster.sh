#!/bin/bash
# ==========================================================================
# setup_cluster.sh — one-shot login-node setup: instantiate + smoke test.
#
#   ./cluster/setup_cluster.sh             # setup + verify only
#   ./cluster/setup_cluster.sh --submit    # setup + verify + sbatch estimation
#
# Run this on the LOGIN node (package downloads need internet; compute nodes
# typically have none). Safe to re-run: instantiate is a no-op when the
# environment is already consistent. If a stale Manifest.toml (e.g. from
# extracting a new bundle over an old directory) breaks package loading, it
# is removed and the environment re-resolved automatically — job 7196 failure
# mode. The estimation is only submitted if the smoke test PASSES.
# ==========================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
echo "Bundle root: $ROOT"

# --- Julia: module cascade -> EasyBuild fallback -> hard fail --------------- #
if ! command -v module &>/dev/null && [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck source=/dev/null
    source /etc/profile.d/modules.sh
fi
if command -v module &>/dev/null; then
    module load Julia/1.11.6 2>/dev/null || module load Julia/1.11 2>/dev/null || \
    module load julia/1.11.6 2>/dev/null || module load julia 2>/dev/null || true
fi
if ! command -v julia &>/dev/null; then
    JB="/repositorio/modules/.local/easybuild/software/Julia/1.11.6-linux-x86_64/bin"
    [ -x "$JB/julia" ] && export PATH="$JB:$PATH"
fi
command -v julia &>/dev/null || { echo "ERROR: Julia not found."; exit 1; }
echo "Julia: $(julia --version)"

# --- same depot as the SLURM jobs (run_estimation.sh) ----------------------- #
export JULIA_DEPOT_PATH="$ROOT/.julia_depot:${HOME}/.julia"
export JULIA_CPU_TARGET="generic"
mkdir -p "$ROOT/.julia_depot"

_instantiate() {
    julia --startup-file=no --project="$ROOT" -e 'import Pkg; Pkg.instantiate()'
}
_smoke() {
    julia --startup-file=no --project="$ROOT" "$ROOT/Model/julia_dynare/smoke_test.jl"
}

echo ""
echo "-- instantiate -------------------------------------------------"
_instantiate

echo ""
echo "-- smoke test --------------------------------------------------"
if ! _smoke; then
    if [ -f "$ROOT/Manifest.toml" ]; then
        echo ""
        echo "Smoke test failed — retrying once with a fresh Manifest.toml"
        echo "(a stale manifest from a previous bundle is the usual cause)."
        rm -f "$ROOT/Manifest.toml"
        _instantiate
        echo ""
        echo "-- smoke test (retry) -------------------------------------"
        _smoke || { echo "ERROR: smoke test still failing — do NOT submit. See output above."; exit 1; }
    else
        echo "ERROR: smoke test failed with a fresh environment — do NOT submit."
        exit 1
    fi
fi

echo ""
echo "================================================================"
echo "  Setup complete — environment verified."
if [ "${1:-}" = "--submit" ]; then
    JOB=$(sbatch --parsable "$ROOT/cluster/run_estimation.sh")
    echo "  Submitted estimation: job $JOB"
    echo "  Monitor:  tail -f nkiosoe_smm_${JOB}.out"
    echo "            tail -f Data/smm_progress_log.csv"
else
    echo "  Submit with:  sbatch cluster/run_estimation.sh"
    echo "  (or rerun this script with --submit)"
fi
echo "================================================================"
