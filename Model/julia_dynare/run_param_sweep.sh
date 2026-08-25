#!/usr/bin/env bash
# =============================================================================
#  run_param_sweep.sh — profile the SMM objective in a pinned elasticity
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
#   ./run_param_sweep.sh epsM                    # default grid for epsM
#   ./run_param_sweep.sh epsY 0.8 0.9 1.0        # custom grid
#
# epsM IS THE PROPAGATION PARAMETER. The variance decomposition of 2026-08-22
# showed the excess comovement (rbar 0.417 model vs 0.221 data) comes from
# network propagation of the twelve INDEPENDENT sectoral TFP shocks: on their
# own they give rbar = 0.585, above the all-shocks total. epsY could not fix
# that because it is the M-vs-V-vs-L elasticity — substitution among a sector's
# own input TYPES. Substitution ACROSS the twelve supplying sectors is epsM,
# pinned at 0.20. That is the one to sweep.
#
# Expect a TRADE-OFF: Atalay's low eps_m is exactly what makes sectoral shocks
# matter (83% -> 21% in his paper), so raising epsM to cut rbar will also cut
# the sectoral-shock share of aggregate volatility. Watch std(Y_i), std(GDP)
# and the 54/41 TFP-vs-external GDP split, not just rbar.
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

PARAM="${1:-epsM}"; shift || true
case "$PARAM" in
    epsY) ENVVAR=SMM_EPSY; DEFAULT_GRID=(0.6 0.7 0.8 0.9 1.0 1.2 1.4) ;;
    epsM) ENVVAR=SMM_EPSM; DEFAULT_GRID=(0.05 0.10 0.20 0.30 0.40 0.50) ;;
    *) echo "unknown parameter '$PARAM' (expected epsY or epsM)"; exit 2 ;;
esac

GRID=("$@")
[ ${#GRID[@]} -eq 0 ] && GRID=("${DEFAULT_GRID[@]}")

mkdir -p logs
SUMMARY="logs/${PARAM}_profile.csv"
echo "${PARAM},objective,cl,rho_om,rho_A,rho_zeta,sigma_zeta,kappaw,lambda_A,lambda_om,rbar_model,status" > "$SUMMARY"

printf '\n=== %s sweep: %s ===\n\n' "$PARAM" "${GRID[*]}"

for e in "${GRID[@]}"; do
    LOG="logs/${PARAM}_${e}.log"
    printf '  %s = %-5s -> %s ... ' "$PARAM" "$e" "$LOG"
    START=$(date +%s)

    # epsY is pinned in BOTH sweeps: the 2026-08-22 profile showed it is not
    # identified (flat 0.8-1.4), so it stays at 0.80 unless it is the swept
    # parameter itself.
    env SMM_FREE_EPSY=0 "$ENVVAR=$e" SMM_TAG="${PARAM}${e}" SMM_STRICT_SENS=0 \
        julia --threads=1 --project=. run_smm_estimation.jl > "$LOG" 2>&1
    RC=$?

    ELAPSED=$(( $(date +%s) - START ))

    EST="estimation_results/smm_estimates_${PARAM}${e}.csv"
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
    RBAR=$(grep -i 'rbar' "estimation_results/smm_results_${PARAM}${e}.csv" \
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
