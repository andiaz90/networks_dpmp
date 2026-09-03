#!/bin/bash
# ==========================================================================
# run_kappa_ladder.sh — sectoral price rigidity kappa_i, estimated
#
# The ladder from PLAN_estimate_sectoral_kappa_2026-09-02.md, in order. Each
# rung is a NESTED restriction of the next, so every outcome is reportable —
# including "kappa is not identified by these moments", which is the outcome
# the lambda_rho experiment produced for sectoral persistence and which is
# worth knowing cheaply rather than after a 14-hour cluster job.
#
# STOP AT THE FIRST RUNG THAT FAILS. Rung 1 exists precisely so that rungs
# 3-5 are not run against a broken wiring.
#
# USAGE
#   bash run_kappa_ladder.sh 1        # one rung
#   bash run_kappa_ladder.sh 1 2 3    # several, in order
#   bash run_kappa_ladder.sh all      # rungs 0-6
#   bash run_kappa_ladder.sh 7        # LOCAL exploratory run (see rung 7)
#
#   KAPPA_LOCAL_HOURS=0.5 KAPPA_LOCAL_MODE=scale bash run_kappa_ladder.sh 7
#
# THREADS. Every rung runs --threads=1. Two independent reasons: ~36% of
# evaluations return a stale decision rule under threading (see
# networks-dpmp-klein-cache-stale-R), and kappa is written through the
# per-thread _SMM_PARAMS vector. A kappa result from a threaded run is not
# trustworthy; the full CMA-ES belongs on the cluster (cluster/run_kappa_smm.sh).
# ==========================================================================
set -euo pipefail

JD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$JD"
J="julia --project=. --threads=1"
LOGS="$JD/logs"; mkdir -p "$LOGS"
STAMP=$(date +%Y%m%d_%H%M%S)

rung() { echo; echo "=============================================================="; \
         echo "  RUNG $1 — $2"; echo "=============================================================="; }

run_rung() {
case "$1" in

# --------------------------------------------------------------------------
0)  rung 0 "Regenerate the moment CSVs with the price-rigidity block (85-109)"
    # Moments 85-109 do not exist in any CSV written before 2026-09-02, and
    # smm_estimation.jl HARD-ERRORS rather than falling back — estimating a
    # price-rigidity parameter with no price-rigidity moment is the failure this
    # whole exercise exists to avoid.
    rm -f ../../Data/sectoral_moments.csv ../../Data/aggregate_moments.csv
    $J compute_data_moments.jl 2>&1 | tee "$LOGS/kappa_r0_moments_$STAMP.log"
    python3 ../../Data/build_reallocation_calibration.py
    echo
    echo ">> Read section 9d in the log above BEFORE going on. It prints the sign"
    echo "   check of autocorr(PH_i) and std(dPH_i)/std(dY_i) against kappa^cal."
    ;;

# --------------------------------------------------------------------------
1)  rung 1 "Wiring verification — nesting, read-back, non-degeneracy, SS invariance"
    # SMM_REPORT_ONLY stops before CMA-ES, so this costs one model solve. The
    # four checks are in smm_estimation.jl under "KAPPA BLOCK VERIFICATION";
    # each corresponds to a bug this codebase has already had in some other
    # parameter. If any of them fails, NOTHING below this line is meaningful.
    SMM_KAPPA_MODE=affine SMM_REPORT_ONLY=1 SMM_WARM_START=0 \
        $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r1_verify_$STAMP.log"
    echo
    echo ">> Required in the log: [1] nesting  [2] read-back  [3] non-degeneracy  [4] SS invariance, all OK."
    echo ">> Also check obj(theta0) here EQUALS the obj(theta0) of an SMM_KAPPA_MODE=off run:"
    echo "   at the seed the two configurations are the same model."
    ;;

# --------------------------------------------------------------------------
2)  rung 2 "Baseline at the CALIBRATED kappa, with the new moments scored"
    # The question this answers, and it is worth answering before freeing
    # anything: does the calibrated (US-pattern) kappa vector already fit the
    # price-rigidity block? If it does, there is nothing to estimate. If it
    # misses, the per-sector table says WHICH sectors are driving it — and that
    # is the substantive result, free of any estimation.
    SMM_KAPPA_MODE=off SMM_REPORT_ONLY=1 SMM_WARM_START=1 \
        $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r2_baseline_$STAMP.log"
    echo
    echo ">> Look at the block lines for 97-108 and 109 in the fit table."
    ;;

# --------------------------------------------------------------------------
3)  rung 3 "Pre-flight sensitivity scan — is kappa identified AT ALL?"
    # THE CHEAPEST WAY TO BE TOLD NO. log(kappaV) was pinned on exactly this
    # evidence: |d obj| = 1.8e-4 against 1.0e+1 for the lambdas, flat to four
    # orders of magnitude. If theta[39]/theta[40] come back at that level, stop
    # here and report flatness — it is a result, and it costs minutes.
    #
    # Accept only if |d obj| is within about one order of magnitude of cl (4.5e-1).
    SMM_KAPPA_MODE=affine SMM_MAX_EVALS=1 SMM_WARM_START=1 \
        $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r3_sens_$STAMP.log"
    grep -E "sensitivity|KAPPA BLOCK|\[1\]|\[2\]|\[3\]|\[4\]" "$LOGS/kappa_r3_sens_$STAMP.log" || true
    ;;

# --------------------------------------------------------------------------
4)  rung 4 "PROFILE lambda_kappa on a grid (level free, dispersion fixed at each point)"
    # PROFILING BEATS ESTIMATING, the epsY lesson. A grid gives the objective
    # PROFILE — which shows whether the parameter is identified and is a figure
    # worth putting in the paper — rather than a point estimate from a search
    # that spent its budget somewhere else.
    #
    # lambda_kappa = 0 is HOMOGENEOUS stickiness, 1 is the calibrated
    # cross-section. Both are on the grid on purpose: they are the two nested
    # tests, and the profile reads off the likelihood-ratio-style comparison.
    for LK in 0.00 0.25 0.50 0.75 1.00 1.25 1.50 2.00; do
        echo "--- lambda_kappa = $LK ---"
        SMM_KAPPA_MODE=scale SMM_LAMBDA_KAPPA="$LK" SMM_WARM_START=1 \
            $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r4_prof_lk${LK}_$STAMP.log"
        cp estimation_results/smm_estimates.csv "estimation_results/smm_estimates_lk${LK}.csv" 2>/dev/null || true
    done
    echo ">> Plot obj against lambda_kappa. A flat profile means the dispersion is not identified;"
    echo "   a minimum at 0 rejects heterogeneous stickiness; a minimum near 1 confirms the calibration."
    ;;

# --------------------------------------------------------------------------
5)  rung 5 "Joint CMA-ES, level + dispersion (the headline estimate)"
    SMM_KAPPA_MODE=affine SMM_WARM_START=1 \
        $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r5_affine_$STAMP.log"
    cp estimation_results/smm_estimates.csv estimation_results/smm_estimates_kappa_affine.csv
    ;;

# --------------------------------------------------------------------------
6)  rung 6 "Eleven free log kappa_i — the sector-by-sector micro comparison"
    # ROBUSTNESS COLUMN, NOT THE HEADLINE. Eleven parameters against 38 usable
    # quarters (44 in the common window, less the pandemic) will overfit. The
    # ridge makes the shrinkage an explicit, reported choice; run the ladder of
    # values and show the estimate is not an artefact of one of them.
    for SH in 0.0 0.5 2.0 10.0; do
        echo "--- ridge = $SH ---"
        SMM_KAPPA_MODE=free SMM_KAPPA_SHRINK="$SH" SMM_WARM_START=1 \
            $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r6_free_sh${SH}_$STAMP.log"
        cp estimation_results/smm_estimates.csv "estimation_results/smm_estimates_kappa_free_sh${SH}.csv" 2>/dev/null || true
    done
    $J compare_kappa_micro.jl 2>&1 | tee "$LOGS/kappa_r6_compare_$STAMP.log"
    ;;

# --------------------------------------------------------------------------
7)  rung 7 "LOCAL exploratory run — bounded, scratch output, safe to interrupt"
    # For watching the model behave, not for producing a number to report.
    #
    # THREE THINGS MAKE THIS LAPTOP-SAFE, and all three matter:
    #
    #  (a) SMM_ESTIMATION_DIR=estimation_results_local. Every write — checkpoint,
    #      best_sol.txt, min_loss.txt, and smm_estimates.csv at the end — goes to
    #      a scratch folder. Without this, an under-converged laptop run
    #      overwrites the stored theta AND becomes "the newest estimate" that
    #      main_SOE_gap.jl picks up for every figure in the paper. The warm start
    #      still reads the REAL checkpoint (fallback in smm_estimation.jl).
    #
    #  (b) SMM_FREE_EPSY=0. epsY shifts the STEADY STATE, so estimating it forces
    #      a nonlinear SS re-solve on essentially every evaluation and the
    #      single-slot SS cache never hits. Pinning it is the codebase's own
    #      standing advice ("the fast path AND the better one") and it is what
    #      makes a laptop run feasible at all. kappa itself does NOT touch the
    #      SS — verification check [4] — so the kappa block costs nothing here.
    #      Free dims: 8 baseline + 2 kappa = 10.
    #
    #  (c) SMM_MAX_HOURS + the stagnation stop. Both unwind gracefully through
    #      the same path and SAVE the best theta found, so Ctrl-C is not needed
    #      and nothing is lost if you walk away.
    #
    # Watch the ms/eval column in the first few lines. If it is above ~2000 the
    # run will not get anywhere useful in an evening — take it to the cluster.
    HOURS="${KAPPA_LOCAL_HOURS:-1.5}"
    MODE="${KAPPA_LOCAL_MODE:-affine}"
    echo "  mode=$MODE  wall-clock cap=${HOURS}h  epsY PINNED  output -> estimation_results_local/"
    SMM_KAPPA_MODE="$MODE" \
    SMM_ESTIMATION_DIR=estimation_results_local \
    SMM_FREE_EPSY=0 \
    SMM_MAX_HOURS="$HOURS" \
    SMM_STAG_EVALS=400 \
    SMM_WARM_START=1 \
        $J run_smm_estimation.jl 2>&1 | tee "$LOGS/kappa_r7_local_${MODE}_$STAMP.log"
    echo
    echo ">> Results in estimation_results_local/ — NOT in estimation_results/."
    echo ">> Compare against the calibrated kappa:"
    echo "   SMM_ESTIMATION_DIR=estimation_results_local $J compare_kappa_micro.jl"
    echo ">> Nothing here is a reportable estimate. It tells you which direction the"
    echo "   level and the dispersion move, and how fast an evaluation is."
    ;;

*)  echo "unknown rung: $1"; exit 1 ;;
esac
}

if [ "${1:-}" = "all" ]; then set -- 0 1 2 3 4 5 6; fi
[ $# -eq 0 ] && { sed -n '1,25p' "$0"; exit 1; }
for r in "$@"; do run_rung "$r"; done
echo; echo "Done. Logs in $LOGS"
