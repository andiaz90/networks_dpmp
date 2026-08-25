#!/usr/bin/env bash
# =============================================================================
#  run_epsY_sweep.sh — profile the SMM objective in epsY
# =============================================================================
#
# WHY A SWEEP RATHER THAN ESTIMATING epsY
# ---------------------------------------
# epsY shifts the STEADY STATE, so estimating it inside CMA-ES forces a
# nonlinear SS re-solve on nearly every evaluation: the SS cache is keyed on
# (epsY, epsM, etastar) and candidates differ in epsY by tiny continuous
# amounts, so it never hits. Pinning epsY and sweeping it solves the SS ONCE per
# grid point, the cache then hits 100%, and each run costs what a normal run
# costs.
#
# It also produces a better object than a point estimate: the objective PROFILE
# in epsY shows whether epsY is identified at all, and is a figure for the paper.
#
# USAGE
#   ./run_epsY_sweep.sh                 # default grid
#   ./run_epsY_sweep.sh 0.8 0.9 1.0     # custom grid
#
# OUTPUT
#   logs/epsY_<v>.log                       full console output per grid point
#   estimation_results/smm_{estimates,results}_epsY<v>.csv
#   logs/epsY_profile.csv                   the profile: epsY, obj, key params
#
# NOTES
#   * SMM_STRICT_SENS=0 so one fragile grid point cannot abort the sweep. Any
#     point that trips it is flagged in the summary — treat its estimate as
#     untrustworthy, not as a data point on the profile.
#   * Grid points are independent. To parallelise, run this script with
#     disjoint grids in separate terminals.

set -uo pipefail
cd "$(dirname "$0")"

GRID=("$@")
if [ ${#GRID[@]} -eq 0 ]; then
    GRID=(0.6 0.7 0.8 0.9 1.0 1.2 1.4)
fi

mkdir -p logs
SUMMARY="logs/epsY_profile.csv"
echo "epsY,objective,cl,rho_om,rho_A,rho_zeta,sigma_zeta,kappaw,lambda_A,lambda_om,rbar_model,status" > "$SUMMARY"

printf '\n=== epsY sweep: %s ===\n\n' "${GRID[*]}"

for e in "${GRID[@]}"; do
    LOG="logs/epsY_${e}.log"
    printf '  epsY = %-5s -> %s ... ' "$e" "$LOG"
    START=$(date +%s)

    SMM_FREE_EPSY=0 SMM_EPSY="$e" SMM_TAG="epsY${e}" SMM_STRICT_SENS=0 \
        julia --threads=1 --project=. run_smm_estimation.jl > "$LOG" 2>&1
    RC=$?

    ELAPSED=$(( $(date +%s) - START ))

    EST="estimation_results/smm_estimates_epsY${e}.csv"
    if [ $RC -ne 0 ] || [ ! -f "$EST" ]; then
        printf 'FAILED (rc=%d, %ds)\n' "$RC" "$ELAPSED"
        echo "$e,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,failed" >> "$SUMMARY"
        continue
    fi

    # obj_hat sits in column 3 of the first data row of smm_estimates.csv.
    OBJ=$(awk -F, 'NR==2{print $3}' "$EST")
    getp() { awk -F, -v k="$1" '$1==k{printf "%s", $2}' "$EST"; }
    # rbar from smm_results.csv, by NAME so a future moment insertion cannot
    # shift it. Careful: the moment label "rbar: avg pairwise corr(Y_i,Y_j)"
    # contains a COMMA inside its quotes, so -F, splits it and the field indices
    # move. Count from the END instead: the last three fields are always
    # data, model, diff — so the model value is $(NF-1).
    RBAR=$(grep -i 'rbar' "estimation_results/smm_results_epsY${e}.csv" \
           | head -1 | awk -F, '{print $(NF-1)}')

    STATUS=ok
    grep -q "FLAT DIRECTIONS" "$LOG" && STATUS=flat_untrustworthy
    grep -q "STALE decision rule" "$LOG" && STATUS=stale_solver

    echo "$e,$OBJ,$(getp ilabcosts),$(getp rho_om),$(getp rho_A),$(getp rho_zeta),$(getp sigma_zeta),$(getp kappaw),$(getp lambda_A),$(getp lambda_om),$RBAR,$STATUS" >> "$SUMMARY"
    printf 'obj=%-10s %-20s (%ds)\n' "$OBJ" "$STATUS" "$ELAPSED"
done

printf '\n=== profile ===\n\n'
column -t -s, "$SUMMARY"
printf '\nFull logs in logs/ ; profile in %s\n' "$SUMMARY"
printf 'Points marked stale_solver or flat_untrustworthy are NOT valid profile points.\n\n'
