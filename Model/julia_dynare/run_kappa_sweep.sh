#!/usr/bin/env bash
# run_kappa_sweep.sh
# =================
# DIAGNOSTIC sweep over SMM_KAPPA_SCALE, a global multiplier on the twelve
# calibrated sectoral Rotemberg price-adjustment costs.
#
# THE QUESTION. std(PH_i) is 9.1% of the SMM objective and corr(Y_i,PH_i) is
# 38.8% — 48% between them. They may share ONE mechanism: model sectoral prices
# barely move, so the output-price correlation collapses toward zero. At the
# 2026-08-25 estimate the model is 3.6x too smooth in Agriculture, 2.9x in
# Manufacturing, 2.4x in Finance. Mining — the only sector with NO Phillips
# curve — is the only one that matches.
#
# WHAT WOULD SETTLE IT. If both blocks improve TOGETHER as kappa falls, they
# share a mechanism and a free kappa scale is worth building. If only std(PH_i)
# moves while corr(Y_i,PH_i) sits still, they do not, and the price-stickiness
# route is dead — which is just as useful to know.
#
# Read it sceptically. The errors go BOTH ways: Construction and Business
# Services are about 2x too VOLATILE, so a single global scale is blunt and the
# blocks may trade off against each other.
#
# This runs SMM_REPORT_ONLY at the current theta-hat — NO estimation, one model
# solve per point, about a minute each. Nothing is re-optimised, so the blocks
# move only because kappa moved.
#
# USAGE
#     ./run_kappa_sweep.sh                 # default grid
#     ./run_kappa_sweep.sh 0.1 0.25 1 4    # custom grid
#
# OUTPUT
#     logs/kappa_<s>.log      per point
#     logs/kappa_profile.csv  obj + the two blocks + std(GDP) + rbar

set -uo pipefail
cd "$(dirname "$0")"
mkdir -p logs

GRID=${*:-"0.10 0.25 0.50 1.00 2.00 4.00"}
WARM=${SMM_WARM_START:-smm_estimates_nonmkt.csv}
SUMMARY=logs/kappa_profile.csv

echo "=== kappa scale sweep: $GRID ==="
echo "    theta held FIXED at $WARM (report-only; nothing is re-estimated)"
echo "kappa_scale,objective,block_PH,block_corrYPH,std_GDP,rbar_model,status" > "$SUMMARY"

for s in $GRID; do
    LOG="logs/kappa_${s}.log"
    printf '  kappa x %-6s -> %s ... ' "$s" "$LOG"
    START=$(date +%s)

    SMM_REPORT_ONLY=1 SMM_FREE_EPSY=0 \
    SMM_KAPPA_SCALE="$s" SMM_WARM_START="$WARM" \
    julia --threads=1 --project=. run_smm_estimation.jl > "$LOG" 2>&1
    RC=$?
    ELAPSED=$(( $(date +%s) - START ))

    if [ $RC -ne 0 ]; then
        printf 'FAILED (rc=%d, %ds)\n' "$RC" "$ELAPSED"
        echo "$s,NA,NA,NA,NA,NA,failed" >> "$SUMMARY"
        continue
    fi

    # obj(theta0) is the report-only objective at the fixed theta.
    OBJ=$(grep -m1 'obj(θ₀)' "$LOG" | sed 's/.*= *//')

    # Block totals come from the "-- <name> ... block  X.XXX  YY.Y% --" headers.
    # Match on the block NAME, never on a column index: moments have been
    # inserted twice in this codebase and every positional read broke silently.
    blk() { grep -m1 -- "-- $1" "$LOG" | sed 's/.*block *//; s/ .*//'; }
    B_PH=$(blk 'std(PH_i)')
    B_CY=$(blk 'corr(Y_i,PH_i)')

    # std(GDP) and rbar: last-but-one field is the MODEL value. The rbar label
    # contains a comma inside quotes, so count from the END, not the start.
    SGDP=$(grep -m1 'std(GDP) ' "$LOG" | awk '{print $3}')
    RBAR=$(grep -m1 'rbar:' "$LOG" | awk '{print $(NF-2)}')

    STATUS=ok
    grep -q 'STALE decision rule' "$LOG" && STATUS=stale_solver
    grep -q 'model failed'        "$LOG" && STATUS=solve_failed

    echo "$s,$OBJ,$B_PH,$B_CY,$SGDP,$RBAR,$STATUS" >> "$SUMMARY"
    printf 'obj=%-10s PH=%-7s corrYPH=%-7s (%ds)\n' "$OBJ" "$B_PH" "$B_CY" "$ELAPSED"
done

echo
echo "=== profile ==="
column -s, -t "$SUMMARY"
echo
echo "Baseline is kappa_scale = 1.00 (obj 6.5303, PH 0.597, corrYPH 2.532)."
echo "VERDICT RULE: if PH and corrYPH fall TOGETHER, build the free scale."
echo "If only PH moves, the two blocks do not share a mechanism — stop here."
