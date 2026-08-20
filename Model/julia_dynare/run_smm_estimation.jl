"""
run_smm_estimation.jl
=====================
Entry point for SMM estimation of the NK-IOSOE Chile model.

HOW TO RUN:
  julia --threads=auto --project=. run_smm_estimation.jl

OUTPUT (all in Model/julia_dynare/estimation_results/ — identical layout on
the cluster; copy that folder back to your machine and main_SOE_gap.jl picks
up the newest θ automatically):
  smm_estimates.csv    — estimated θ (auto-loaded by main_SOE_gap.jl)
  smm_results.csv      — full 60-moment fit table
  smm_checkpoint.csv   — warm-start checkpoint (updated live, every improvement)
  best_sol.txt         — current best θ, human-readable (name  value)
  min_loss.txt         — objective at the current best θ
  smm_progress_log.csv — objective trajectory, one row per 50 evaluations
  smm_inference.csv    — standard errors (after a completed run)
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
    # CONTEXT comes from the top-level @dynare below, NOT from deserialisation
    # (2026-08-19). See the long note there: a deserialised context cannot be
    # solved, because Dynare's solver calls runtime-generated functions held in
    # the module-level Dynare.DFunctions, which only @dynare populates.
    context = _DYNARE_CONTEXT[]
    context === nothing && error("""
        Dynare context not built. The top-level @dynare in this file must run
        before _main(); see the DFunctions note above.
        """)

    @printf "--- Step 2: Model loaded via @dynare (DFunctions populated) ---\n\n"

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

# =========================================================================== #
#  BUILD THE DYNARE CONTEXT AT TOP LEVEL  (2026-08-19)                        #
# =========================================================================== #
# Why not deserialize mod/nk_iosoe_context.jls, as this script did until today:
#
# Dynare's compute_first_order_solution! ultimately calls
# Dynare.DFunctions.dynamic_resid! / dynamic_g1!, which are RuntimeGenerated
# functions living in a MODULE-LEVEL namespace. They are populated only when
# @dynare parses a .mod file IN THIS PROCESS. Deserialising a context restores
# the data but not those functions, so the call hit whatever model DFunctions
# happened to hold — a 6-equation leftover from Dynare's own precompilation:
#
#     AssertionError: length(residual) == 6      (our model has 591)
#
# With Dynare's solver unusable, resolve_first_order! fell through to the
# hand-rolled Klein path (broken: see smm_model_moments.jl) and then to its last
# resort, returning the decision rule already in the context — at θ_baseline.
# That returns success while ignoring θ, which is why the objective was exactly
# constant in all 7 free parameters and why every estimation run since July
# "converged" instantly having estimated nothing.
#
# Running @dynare here costs ~30-60 s once, against multi-hour estimation runs.
# It must be at TOP LEVEL: inside a function it hits world-age errors (the same
# reason run_dynare_subprocess.jl does it this way).
#
# NOTE this recompiles the model, so params_jl.mod must already be current —
# i.e. run main_SOE_gap.jl first, exactly as before.
const _DYNARE_CONTEXT = Ref{Any}(nothing)

# EVERY statement below must stay at TOP LEVEL, each its own statement.
# `@dynare` runs the Dynare preprocessor at MACRO-EXPANSION time, not at
# runtime. Julia expands an entire block (let/function/begin) before executing
# any of it, so wrapping this in `let` meant the preprocessor ran while the
# working directory was still the caller's — "Could not open file:
# NK_SOE_lev_gap2.mod" — and @dynare then returned a Vector instead of a
# Context. run_dynare_subprocess.jl keeps these as separate top-level
# statements for exactly this reason; do not "tidy" them into a block.
_SMM_MOD_DIR = joinpath(SCRIPT_DIR, "mod")

isfile(joinpath(_SMM_MOD_DIR, "params_jl.mod")) || error("""
    mod/params_jl.mod not found. Run main_SOE_gap.jl first — it writes the
    parameter file that NK_SOE_lev_gap2.mod includes.
    """)

_SMM_VERBOSE = get(ENV, "DYNARE_VERBOSE", "0") in ("1", "true")
_SMM_LOG     = joinpath(_SMM_MOD_DIR, "dynare_smm_load.log")
_SMM_OLD_PWD = pwd()

@printf "--- Loading model via @dynare (populates Dynare.DFunctions) ---\n"
@printf "    output → %s%s\n" basename(_SMM_LOG) (_SMM_VERBOSE ? "  [DYNARE_VERBOSE=1: inline]" : "")

cd(_SMM_MOD_DIR)

_SMM_STDOUT = stdout
_SMM_STDERR = stderr
_SMM_LOG_IO = _SMM_VERBOSE ? nothing : open(_SMM_LOG, "w")
if !_SMM_VERBOSE
    redirect_stdout(_SMM_LOG_IO)
    redirect_stderr(_SMM_LOG_IO)
end

_SMM_CTX = try
    @dynare "NK_SOE_lev_gap2"
finally
    if !_SMM_VERBOSE
        redirect_stdout(_SMM_STDOUT)
        redirect_stderr(_SMM_STDERR)
        close(_SMM_LOG_IO)
    end
    cd(_SMM_OLD_PWD)
end

_SMM_CTX isa Dynare.Context || error("""
    @dynare returned $(typeof(_SMM_CTX)), not a Dynare.Context.
    The preprocessor almost certainly failed — see $(basename(_SMM_LOG)).
    """)
_DYNARE_CONTEXT[] = _SMM_CTX

@printf "    done — %d endogenous, %d exogenous\n\n" length(
    _SMM_CTX.results.model_results[1].trends.endogenous_steady_state) _SMM_CTX.models[1].exogenous_nbr

Base.invokelatest(_main)
