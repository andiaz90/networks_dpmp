"""
run_smm_estimation.jl
=====================
Entry point for SMM estimation of the NK-IOSOE Chile model.

HOW TO RUN:
  julia --threads=auto --project=. run_smm_estimation.jl

FIXES APPLIED (relative to original):
  1. HP filter applied to model moments (matches data HP-filtered log deviations)
  2. Monetary policy shock (eps_i) activated in Σe (was missing)
  3. Lag-1 covariance formula corrected (spurious B·Σ·R' term removed)
  4. Tighter, economically motivated parameter bounds
  5. Proportional weighting + 5x weight on rank correlations
  6. SS caching with warm-start tolerance
  7. Pre-allocated scratch arrays, in-place Lyapunov, @views

ESTIMATED PARAMETERS (23-element vector θ):
  1   ilabcosts        aggregate labour adjustment cost
  2   epsY             elast. of subst. in production
  3   epsM             elast. of subst. between materials
  4   log(kappaV)      log of import price adj. cost
  5   rho_om           AR persistence, goods-services shock
  6   sigma_om         std dev, goods-services shock
  7   rho_A            AR persistence, TFP shocks (common)
  8–19 isigma_tfp_i   std dev of sector i TFP shock
  20  rho_pvstar       AR persistence, import price shock
  21  sigma_pvstar     std dev, import price shock
  22  rho_xi           AR persistence, preference shock
  23  sigma_xi         std dev, preference shock

OUTPUT:
  Data/smm_estimates.csv    — estimated θ (auto-loaded by main_SOE_gap.jl)
  Data/smm_results.csv      — full moment-fit table
  Data/smm_checkpoint.csv   — CMA-ES warm-start checkpoint
"""

# =========================================================================== #
#  Top-level: packages + includes (no world-age issues)                       #
# =========================================================================== #

using CSV, DataFrames, Printf, Serialization, Dynare

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

# Include all modules at TOP LEVEL so they share the same world
include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))        # build_baseline, default_theta0, etc.
include(joinpath(SCRIPT_DIR, "smm_estimation_v2.jl"))      # smm_run_v2 + optimized moments


# =========================================================================== #
#  _main() — wrapped for Julia 1.12 world-age compatibility                  #
# =========================================================================== #

function _main()

    @printf "\n%s\n" repeat("=", 60)
    @printf "  SMM ESTIMATION — NK-IOSOE Chile Model\n"
    @printf "  (HP-filtered moments, monetary shock active, corrected Γ₁)\n"
    @printf "%s\n\n" repeat("=", 60)

    # ---- Step 1: Check data moments ----------------------------------------
    sec_ok = isfile(joinpath(DATA_DIR, "sectoral_moments.csv"))
    agg_ok = isfile(joinpath(DATA_DIR, "aggregate_moments.csv"))

    if !sec_ok || !agg_ok
        error("""
        Data moment files not found in $(DATA_DIR).
        Run bootstrap_csv.jl first to generate them from existing .mat files:
          julia --project=. bootstrap_csv.jl
        Or run compute_data_moments.jl to recompute from raw Excel/CSV sources.
        """)
    end
    @printf "--- Step 1: Data moments OK ---\n"
    @printf "  %s\n" joinpath(DATA_DIR, "sectoral_moments.csv")
    @printf "  %s\n\n" joinpath(DATA_DIR, "aggregate_moments.csv")

    # ---- Step 2: Load the Dynare context -----------------------------------
    main_ctx = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls")
    smm_ctx  = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_smm_context.jls")
    context_file = isfile(main_ctx) ? main_ctx : smm_ctx

    if !isfile(context_file)
        error("""
        No Dynare context found.  Run main_SOE_gap.jl first:
          julia --project=. main_SOE_gap.jl
        This creates $(main_ctx) which contains the first-order decision rule.
        """)
    end

    @printf "--- Step 2: Loading context ---\n"
    @printf "  %s\n\n" context_file
    context = deserialize(context_file)

    ss_check = context.results.model_results[1].trends.endogenous_steady_state
    n_endo   = length(ss_check)
    @printf "  Context loaded: %d endogenous variables\n" n_endo

    if n_endo < 400
        @printf "\n  ERROR: context has only %d variables (expected ~491).\n" n_endo
        error("Context has $(n_endo) variables, expected ≥400.")
    end

    # Read endo_names from the CSV written by the Dynare subprocess.
    MOD_DIR = joinpath(SCRIPT_DIR, "mod")
    endo_names_file = joinpath(MOD_DIR, "dynare_endo_names.csv")
    if isfile(endo_names_file)
        endo_names_csv = CSV.read(endo_names_file, DataFrame)
        endo_names_override = String.(endo_names_csv.variable)
        @printf "  endo_names loaded from CSV: %d variables\n\n" length(endo_names_override)
    else
        endo_names_override = nothing
        @printf "  endo_names CSV not found — will use symboltable\n\n"
    end

    # ---- Step 3: Run SMM estimation ----------------------------------------
    @printf "--- Step 3: Running SMM estimation ---\n\n"
    theta_hat, obj_hat, moments_hat = Base.invokelatest(
        smm_run_v2, context; endo_names_override=endo_names_override)

    # ---- Step 4: Summary ---------------------------------------------------
    if isnan(obj_hat)
        @printf "\n  Estimation did not complete — see messages above.\n\n"
        return
    end

    @printf "\n%s\n  ESTIMATION COMPLETE\n%s\n\n" repeat("=",60) repeat("=",60)
    @printf "  Objective at θ̂:  %.6f\n" obj_hat
    @printf "  Estimates:  %s\n" joinpath(DATA_DIR, "smm_estimates.csv")
    @printf "  Results:    %s\n" joinpath(DATA_DIR, "smm_results.csv")
    @printf "\n  Apply estimates: re-run main_SOE_gap.jl (EXERCISE=0)\n"
    @printf "  smm_estimates.csv is loaded automatically.\n\n"

end  # function _main()

Base.invokelatest(_main)
