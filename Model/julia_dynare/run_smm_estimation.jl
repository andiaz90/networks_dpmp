"""
run_smm_estimation.jl
=====================
Entry point for SMM estimation of the NK-IOSOE Chile model.

HOW TO RUN:
  julia --project=. run_smm_estimation.jl

WHAT IT DOES:
  Step 1 — Verifies data moments exist (run bootstrap_csv.jl first if missing)
  Step 2 — Loads the compiled Dynare context from mod/nk_iosoe_smm_context.jls
            (created automatically when main_SOE_gap.jl runs)
  Step 3 — Runs SMM estimation (CMA-ES optimizer, ~100k model evaluations)
  Step 4 — Saves Data/smm_estimates.csv for main_SOE_gap.jl to load


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
#  Top-level: packages + includes (no world-age issues here)                  #
# =========================================================================== #

using CSV, DataFrames, Printf, Serialization, Dynare

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

# Include all modules at TOP LEVEL so they share the same world
include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))


# =========================================================================== #
#  _main() — wrapped for Julia 1.12 world-age compatibility                  #
# =========================================================================== #

function _main()

    @printf "\n%s\n" repeat("=", 60)
    @printf "  SMM ESTIMATION — NK-IOSOE Chile Model\n"
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
    # Strategy: prefer nk_iosoe_context.jls (main context, written by
    # main_SOE_gap.jl with stoch_simul).  Even though stoch_simul's gees call
    # fails on ARM Mac and corrupts the symboltable (to 6 names), the
    # model_results[1] still has 491 variables AND lre.g1_1/g1_2 are populated
    # by the first-order perturbation step BEFORE gees is called.  Those
    # non-zero decision rule matrices are what the cached-fallback path in
    # resolve_first_order! needs.
    #
    # The SMM context (nk_iosoe_smm_context.jls, no stoch_simul) has a clean
    # symboltable but zero lre.g1_1/g1_2 (first-order solution never computed).
    #
    # endo_names bypass (from dynare_endo_names.csv) avoids both symboltable
    # issues regardless of which context is loaded.
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

    ctx_path = context_file

    @printf "--- Step 2: Loading context ---\n"
    @printf "  %s\n\n" ctx_path
    context = deserialize(ctx_path)

    # Note: context.symboltable may have fewer names than model results (Dynare.jl issue:
    # symboltable update fails when gees errors during check; even though
    # context.results.model_results[1] correctly has 491 variables).
    # Validate using the model results, not the symbol table.
    ss_check = context.results.model_results[1].trends.endogenous_steady_state
    n_endo   = length(ss_check)
    @printf "  Context loaded: %d endogenous variables (from model results)\n" n_endo

    if n_endo < 400
        @printf "\n  ERROR: context has only %d variables (expected ~491).\n" n_endo
        @printf "  Fix: delete context file and re-run:\n"
        @printf "    rm %s\n" ctx_path
        @printf "    julia --project=. main_SOE_gap.jl\n"
        @printf "    julia --project=. run_smm_estimation.jl\n\n"
        error("Context has $(n_endo) variables, expected ≥400.")
    end

    # Read endo_names from the CSV written by the Dynare subprocess.
    # We do NOT use Dynare.get_endogenous(context.symboltable) because it may be inconsistent.
    # the symboltable may be inconsistent (6 names vs 491 in the results).
    MOD_DIR_smm = joinpath(SCRIPT_DIR, "mod")
    endo_names_file = joinpath(MOD_DIR_smm, "dynare_endo_names.csv")
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
    # Pass endo_names_override so smm_run uses CSV names instead of broken symboltable
    theta_hat, obj_hat, moments_hat = smm_run(context; endo_names_override=endo_names_override)

    # ---- Step 4: Summary ---------------------------------------------------
    if isnan(obj_hat)
        # Estimation did not complete.
        # smm_run already printed the reason above — nothing more to do here.
        @printf "\n  Estimation did not complete — see messages above.\n"
        @printf "  main_SOE_gap.jl will use hard-coded defaults until\n"
        @printf "  Data/smm_estimates.csv is provided from an Intel/x86 run.\n\n"
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
