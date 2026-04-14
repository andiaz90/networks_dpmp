"""
run_smm_estimation.jl
=====================
Entry point for SMM estimation of the NK-IOSOE Chile model.

HOW TO RUN:
  julia --project=. run_smm_estimation.jl

WHAT IT DOES:
  Step 1 — Compute data moments (skipped if CSVs already exist)
  Step 2 — Run main_SOE_gap.jl with EXERCISE=0 (Baseline)
            This computes the steady state, writes params_jl.mod,
            runs Dynare, and saves the compiled context to
            mod/nk_iosoe_context.jls
  Step 3 — Load the context and run SMM estimation (CMA-ES)
            Estimated parameters saved to Data/smm_estimates.csv
  Step 4 — Re-run main_SOE_gap.jl (EXERCISE=0) with the estimates
            to verify the moment fit

ESTIMATED PARAMETERS (23-element vector θ):
  1   ilabcosts        aggregate labour adjustment cost
  2   epsY             elast. of subst. in production
  3   epsM             elast. of subst. between materials
  4   log(kappaV)      log of import price adj. cost
  5   rho_om           AR persistence, goods-services shock
  6   sigma_om         std dev, goods-services shock
  7   rho_A            AR persistence, TFP shocks
  8–19 isigma_tfp_i   std dev, sector i TFP shock
  20  rho_pvstar       AR persistence, import price shock
  21  sigma_pvstar     std dev, import price shock
  22  rho_xi           AR persistence, preference shock
  23  sigma_xi         std dev, preference shock

OUTPUT:
  Data/smm_estimates.csv    — estimated θ (read by main_SOE_gap.jl)
  Data/smm_results.csv      — full moment-fit table
  Data/smm_checkpoint.csv   — rolling best-so-far (warm start)

WARM START:
  If Data/smm_estimates.csv or Data/smm_checkpoint.csv already exists,
  the optimizer starts from the previous best solution.
  Delete those files to restart from scratch.
"""

SCRIPT_DIR = @__DIR__

# =========================================================================== #
#  STEP 1: Check data moments                                                 #
# =========================================================================== #

using CSV, DataFrames, Printf

DATA_DIR = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

sec_ok = isfile(joinpath(DATA_DIR, "sectoral_moments.csv"))
agg_ok = isfile(joinpath(DATA_DIR, "aggregate_moments.csv"))

if !sec_ok || !agg_ok
    @printf "--- Step 1: Computing data moments ---\n"
    include(joinpath(SCRIPT_DIR, "compute_data_moments.jl"))
else
    @printf "--- Step 1: Data moments already available ---\n"
    @printf "  %s\n" joinpath(DATA_DIR, "sectoral_moments.csv")
    @printf "  %s\n\n" joinpath(DATA_DIR, "aggregate_moments.csv")
end

# =========================================================================== #
#  STEP 2: Compile model with EXERCISE=0 (Baseline)                          #
# =========================================================================== #

context_file = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls")

if isfile(context_file)
    @printf "--- Step 2: Compiled Dynare context already exists ---\n"
    @printf "  %s\n\n" context_file
    @printf "  (Delete this file and re-run to recompile with new parameters)\n\n"
else
    @printf "--- Step 2: Compiling model (EXERCISE=0, Baseline) ---\n"
    @printf "  This runs main_SOE_gap.jl and takes ~3 minutes on first run.\n\n"

    # Run main_SOE_gap.jl with EXERCISE=0 via subprocess to compile the model
    julia_exe = joinpath(Sys.BINDIR, "julia")
    project   = dirname(Base.active_project())
    main_script = joinpath(SCRIPT_DIR, "main_SOE_gap.jl")

    # Override EXERCISE to 0 via environment variable
    ENV["SMM_EXERCISE"] = "0"
    run(`$julia_exe --project=$project $main_script`)
    delete!(ENV, "SMM_EXERCISE")

    if !isfile(context_file)
        error("""
        Context file not found after running main_SOE_gap.jl:
          $context_file
        Something went wrong — check the output above for errors.
        """)
    end
    @printf "\n--- Model compiled and context saved ---\n\n"
end

# =========================================================================== #
#  STEP 3: Run SMM estimation                                                 #
# =========================================================================== #

@printf "--- Step 3: SMM Estimation ---\n\n"

using Dynare, Serialization

context = deserialize(context_file)
@printf "  Context loaded: %d endogenous variables\n\n" length(Dynare.get_endogenous(context.symboltable))

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))

theta_hat, obj_hat, moments_hat = smm_run(context)

# =========================================================================== #
#  STEP 4: Summary                                                            #
# =========================================================================== #

@printf "\n%s\n" repeat("=", 60)
@printf "  ESTIMATION COMPLETE\n"
@printf "%s\n\n" repeat("=", 60)
@printf "  Objective at θ̂:  %.6f\n" obj_hat
@printf "  Estimates saved: %s\n" joinpath(DATA_DIR, "smm_estimates.csv")
@printf "\n  Next step: re-run the model with the estimates applied:\n"
@printf "    julia --project=. main_SOE_gap.jl   (EXERCISE=0)\n\n"
@printf "  smm_estimates.csv is automatically loaded by main_SOE_gap.jl.\n\n"
