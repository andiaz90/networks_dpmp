"""
run_smm_estimation_v2.jl
========================
Entry point for optimized SMM estimation (v2).

HOW TO RUN:
  julia --threads=auto --project=. run_smm_estimation_v2.jl

KEY IMPROVEMENTS OVER v1:
  1. SS caching with warm-start tolerance (3-5x fewer nlsolve calls)
  2. Tighter, economically motivated bounds (fewer BK failures → fewer wasted evals)
  3. Proportional weighting (percentage errors matter equally)
  4. 5x weight on rank correlations (drives the sectoral narrative)
  5. Objective decomposition at pre-flight and final results
"""

using CSV, DataFrames, Printf, Serialization, Dynare

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

# Include all modules at TOP LEVEL (world-age safety)
include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))   # needed for resolve_first_order!, set_param!, etc.
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))       # needed for build_baseline, default_theta0, etc.
include(joinpath(SCRIPT_DIR, "smm_estimation_v2.jl"))    # v2 additions


function _main_v2()
    @printf "\n%s\n" repeat("=", 60)
    @printf "  SMM ESTIMATION v2 — NK-IOSOE Chile Model\n"
    @printf "%s\n\n" repeat("=", 60)

    # Step 1: Check data
    sec_ok = isfile(joinpath(DATA_DIR, "sectoral_moments.csv"))
    agg_ok = isfile(joinpath(DATA_DIR, "aggregate_moments.csv"))
    if !sec_ok || !agg_ok
        error("Data moment files not found in $(DATA_DIR). Run bootstrap_csv.jl or compute_data_moments.jl first.")
    end
    @printf "--- Step 1: Data moments OK ---\n\n"

    # Step 2: Load context
    main_ctx = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls")
    smm_ctx  = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_smm_context.jls")
    context_file = isfile(main_ctx) ? main_ctx : smm_ctx
    !isfile(context_file) && error("No Dynare context. Run main_SOE_gap.jl first.")

    @printf "--- Step 2: Loading context from %s ---\n\n" basename(context_file)
    context = deserialize(context_file)

    ss_check = context.results.model_results[1].trends.endogenous_steady_state
    @printf "  %d endogenous variables\n" length(ss_check)
    length(ss_check) < 400 && error("Context has $(length(ss_check)) vars, expected ≥400.")

    # Load endo_names from CSV
    endo_names_file = joinpath(SCRIPT_DIR, "mod", "dynare_endo_names.csv")
    endo_names_override = if isfile(endo_names_file)
        en = CSV.read(endo_names_file, DataFrame)
        @printf "  endo_names: %d variables (from CSV)\n\n" nrow(en)
        String.(en.variable)
    else
        @printf "  endo_names: using symboltable\n\n"
        nothing
    end

    # Step 3: Run
    @printf "--- Step 3: Running SMM v2 ---\n\n"
    θ_hat, obj_hat, moments_hat = Base.invokelatest(
        smm_run_v2, context; endo_names_override=endo_names_override)

    if isnan(obj_hat)
        @printf "\n  Estimation did not complete.\n\n"
        return
    end

    @printf "\n%s\n  ESTIMATION COMPLETE\n%s\n\n" repeat("=",60) repeat("=",60)
    @printf "  Objective: %.6f\n" obj_hat
    @printf "  Estimates: %s\n" joinpath(DATA_DIR, "smm_estimates.csv")
    @printf "  Results:   %s\n\n" joinpath(DATA_DIR, "smm_results_v2.csv")
end

Base.invokelatest(_main_v2)
