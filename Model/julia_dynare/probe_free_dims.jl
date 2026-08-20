"""
probe_free_dims.jl
==================
Diagnose why CMA-ES terminates immediately and why every Jacobian column is zero.

SYMPTOMS (2026-08-19 run):
  - 126 evaluations in 18.4 s against a 4 h / 300,000-eval budget
  - best_obj frozen at the θ₀ value; all 7 free params returned unchanged
  - smm_inference: std err exactly 0.0 and t = Inf for all 7 FREE params
  - Klein cache = 0%  (so the solver IS re-solving each call)
  - identical signature in the July run (100 evals, 20 s)

Both the optimizer stall and the zero Jacobian are explained if the objective
does not respond to the free parameters. This walks each free dimension across
its whole admissible range and prints the objective, so we can tell:

  (a) objective moves        -> optimizer/termination problem
  (b) objective flat in ALL  -> θ never reaches the model (wiring bug)
  (c) objective flat in SOME -> those dims are genuinely unidentified

Run AFTER main_SOE_gap.jl (needs mod/nk_iosoe_context.jls):
    julia --project=. probe_free_dims.jl
"""

using CSV, DataFrames, Printf, Serialization, Dynare, LinearAlgebra

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))

function _main()
    ctx_path = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls")
    isfile(ctx_path) || error("No context. Run main_SOE_gap.jl first.")
    context = deserialize(ctx_path)

    en_file = joinpath(SCRIPT_DIR, "mod", "dynare_endo_names.csv")
    endo_names = isfile(en_file) ? String.(CSV.read(en_file, DataFrame).variable) :
                                   Dynare.get_endogenous(context.symboltable)

    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)
    W = build_weighting_matrix(data_moments)

    # θ₀ exactly as the run left it. Nothing moved, so smm_estimates.csv IS θ₀,
    # which avoids re-deriving it (and any drift from smm_run's seeding block).
    est = CSV.read(joinpath(ESTIMATION_DIR, "smm_estimates.csv"), DataFrame)
    nrow(est) == N_THETA || error("smm_estimates.csv has $(nrow(est)) rows, expected $(N_THETA)")
    θ0 = Float64.(est.value)
    SMM_PIN[] = copy(θ0)   # pinned entries read from here, exactly as in smm_run

    obj(θ) = begin
        m, ok = smm_model_moments(θ, context, baseline, endo_names)
        (!ok || any(isnan, m)) && return NaN
        ψ = data_moments .- m
        dot(ψ, W * ψ)
    end

    o0 = obj(θ0)
    @printf "\n%s\n  FREE-DIMENSION PROBE\n%s\n\n" repeat("=", 66) repeat("=", 66)
    @printf "  obj(θ₀) = %.6f\n\n" o0

    for p in FREE_THETA
        lo, hi = LB[p], UB[p]
        grid = [lo + f*(hi-lo) for f in (0.0, 0.02, 0.10, 0.25, 0.50, 0.75, 1.0)]
        @printf "  θ[%d]  %-20s  (θ₀ = %.4f,  range [%.4g, %.4g])\n" p PARAM_LABELS[p] θ0[p] lo hi
        vals = Float64[]
        for g in grid
            θ = copy(θ0); θ[p] = g
            v = obj(θ); push!(vals, v)
            @printf "      %-12.5g -> %s\n" g (isnan(v) ? "   FAILED" : @sprintf("%12.6f  (Δ = %+.6f)", v, v - o0))
        end
        fin = filter(isfinite, vals)
        spread = isempty(fin) ? 0.0 : maximum(fin) - minimum(fin)
        @printf "      spread over the full range: %.6f  %s\n\n" spread (
            spread < 1e-8 ? "<<<< FLAT — this parameter does not reach the model" :
            spread < 1e-3 ? "<<   nearly flat — weakly identified" : "")
    end

    @printf "%s\n" repeat("=", 66)
    @printf "  If EVERY dimension is flat, θ is not reaching the solver and the\n"
    @printf "  problem is wiring, not optimisation. If some move, CMA-ES is\n"
    @printf "  stalling on the flat ones and the free set should be trimmed.\n"
    @printf "%s\n\n" repeat("=", 66)
end

Base.invokelatest(_main)
