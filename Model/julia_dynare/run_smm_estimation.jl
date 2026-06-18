"""
run_smm_estimation.jl
=====================
Entry point for SMM estimation of the NK-IOSOE Chile model.

HOW TO RUN:
  julia --threads=auto --project=. run_smm_estimation.jl

OUTPUT:
  Data/smm_estimates.csv    — estimated θ (auto-loaded by main_SOE_gap.jl)
  Data/smm_results.csv      — full moment-fit table
  Data/smm_checkpoint.csv   — warm-start checkpoint for subsequent runs
"""

using CSV, DataFrames, Printf, Serialization, Dynare

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

# Include all modules at TOP LEVEL (Julia 1.12 world-age safety)
include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))
include(joinpath(SCRIPT_DIR, "smm_inference.jl"))


function _main()
    @printf "\n%s\n  SMM ESTIMATION — NK-IOSOE Chile Model\n%s\n\n" repeat("=",60) repeat("=",60)

    # Step 1: Check data moments
    sec_ok = isfile(joinpath(DATA_DIR, "sectoral_moments.csv"))
    agg_ok = isfile(joinpath(DATA_DIR, "aggregate_moments.csv"))
    if !sec_ok || !agg_ok
        error("""
        Data moment files not found in $(DATA_DIR).
        Run bootstrap_csv.jl first:
          julia --project=. bootstrap_csv.jl
        Or compute from scratch:
          julia --project=. compute_data_moments.jl
        """)
    end
    @printf "--- Step 1: Data moments OK ---\n\n"

    # Step 2: Load Dynare context.
    # Single unified model now (NK_SOE_lev_gap2.mod); the separate _smm context
    # was retired. main_SOE_gap.jl compiles the unified model to this context.
    ctx_path = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls")
    !isfile(ctx_path) && error("""
        No Dynare context found at $(ctx_path). Run main_SOE_gap.jl first:
          julia --project=. main_SOE_gap.jl
        """)

    @printf "--- Step 2: Loading context from %s ---\n\n" basename(ctx_path)
    context = deserialize(ctx_path)

    n_endo = length(context.results.model_results[1].trends.endogenous_steady_state)
    @printf "  %d endogenous variables\n" n_endo
    n_endo < 400 && error("Context has $(n_endo) variables, expected ≥400.")

    # Load endo_names from CSV (bypasses potentially broken symboltable)
    endo_names_file = joinpath(SCRIPT_DIR, "mod", "dynare_endo_names.csv")
    endo_names_override = if isfile(endo_names_file)
        en = CSV.read(endo_names_file, DataFrame)
        @printf "  endo_names: %d variables (from CSV)\n\n" nrow(en)
        String.(en.variable)
    else
        @printf "  endo_names: using symboltable\n\n"
        nothing
    end

    # Step 3: Run estimation
    @printf "--- Step 3: Running SMM estimation ---\n\n"
    theta_hat, obj_hat, moments_hat = Base.invokelatest(
        smm_run, context; endo_names_override=endo_names_override)

    isnan(obj_hat) && (@printf "\n  Estimation did not complete.\n\n"; return)

    @printf "\n%s\n  COMPLETE — obj=%.6f\n%s\n\n" repeat("=",60) obj_hat repeat("=",60)
end

Base.invokelatest(_main)
