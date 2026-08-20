"""
probe_chain.jl
==============
The objective is flat in all 7 free parameters even after syncing
_SMM_PARAMS -> context.work.params. So the break is somewhere else in the chain

    set_param!  ->  _SMM_PARAMS  ->  Jacobian g1  ->  T,R  ->  Γ  ->  moments

This tests each link separately for ONE parameter (kappaw, θ[36]) at two very
different values (0 = flexible wages, 400 = very sticky). kappaw is chosen
because main_SOE_gap.jl demonstrably responds to it (objective 13.40 / 14.86 /
15.82 across 0 / 100 / 400), so any flatness here is definitely a bug and not
economics.

    julia --project=. probe_chain.jl
"""

using CSV, DataFrames, Printf, Serialization, Dynare, LinearAlgebra

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))

_sum(A) = isempty(A) ? NaN : sum(abs, A)

function _main()
    context = deserialize(joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls"))
    en = joinpath(SCRIPT_DIR, "mod", "dynare_endo_names.csv")
    endo_names = String.(CSV.read(en, DataFrame).variable)
    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)
    W = build_weighting_matrix(data_moments)

    est = CSV.read(joinpath(ESTIMATION_DIR, "smm_estimates.csv"), DataFrame)
    θ0 = Float64.(est.value)
    SMM_PIN[] = copy(θ0)

    @printf "\n%s\n  CHAIN PROBE — kappaw (θ[36])\n%s\n\n" repeat("=",70) repeat("=",70)

    idx = param_idx(context, "kappaw")
    @printf "  param_idx(context, \"kappaw\") = %s\n" string(idx)
    idx === nothing && (@printf "  >>> set_param! is a SILENT NO-OP: the name is not in the table.\n"; return)

    for kw in (0.0, 400.0)
        @printf "\n  ---- kappaw = %.1f ----\n" kw
        # LINK 1: set_param! -> _SMM_PARAMS
        set_param!(context, "kappaw", kw)
        got = get_param_val(context, "kappaw")
        @printf "   1. set_param! -> get_param_val   : %.4f  %s\n" got (
            isapprox(got, kw) ? "ok" : ">>> NOT WRITTEN")
        @printf "      _SMM_PARAMS[tid][idx]         : %.4f\n" _SMM_PARAMS[_tid()][][idx]
        wp = try context.work.params catch; Float64[] end
        @printf "      context.work.params[idx]      : %s\n" (
            isempty(wp) ? "unreadable" : @sprintf("%.4f", wp[idx]))

        # LINK 2: params -> Jacobian
        Gj = _eval_dynamic_jacobian(context)
        @printf "   2. sum|dynamic Jacobian|        : %.10e\n" _sum(Gj)

        # LINK 3: Jacobian -> decision rule
        ok, T, R, _ = resolve_first_order!(context)
        @printf "   3. resolve_first_order! ok=%-5s  sum|T|=%.10e  sum|R|=%.10e\n" ok _sum(T) _sum(R)

        # LINK 4: full objective through the estimator's own entry point
        θ = copy(θ0); θ[36] = kw
        m, okm = smm_model_moments(θ, context, baseline, endo_names)
        if okm && !any(isnan, m)
            ψ = data_moments .- m
            @printf "   4. smm_model_moments -> obj     : %.8f\n" dot(ψ, W*ψ)
        else
            @printf "   4. smm_model_moments FAILED (ok=%s)\n" okm
        end
    end

    @printf "\n%s\n" repeat("=",70)
    @printf "  Read DOWN the chain: the first link whose number does NOT differ\n"
    @printf "  between the two kappaw values is where θ stops propagating.\n"
    @printf "  stale-decision-rule fallbacks so far: %d\n" _STALE_DR_USES[]
    @printf "%s\n\n" repeat("=",70)
end

Base.invokelatest(_main)
