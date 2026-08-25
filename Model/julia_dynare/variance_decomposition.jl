"""
variance_decomposition.jl
=========================
Which shocks generate the model's sectoral output volatility — and, above all,
which generate its CROSS-SECTORAL COMOVEMENT?

WHY THIS EXISTS
---------------
The 84-moment fit has one large, stubborn miss that is not a fitting problem but
a structural one:

    rbar = mean_{i<j} corr(Y_i, Y_j)     model ~0.42     data 0.221

The model produces about twice the cross-sectoral comovement in the Chilean
data. The first hypothesis was that input substitution was too low — a
complementary network propagates sector-specific shocks, so raising epsY should
cut rbar. The epsY sweep of 2026-08-22 REFUTED that: across epsY = 0.7 … 1.4,
rbar moves only 0.457 → 0.417 → 0.416 → 0.423 → 0.427. It bottoms out around
0.416 and turns back up, nowhere near 0.221.

So the comovement comes from shocks that hit every sector at once, not from
propagation of sector-specific ones. The candidates are all PINNED, which is why
the estimation could never fix this: sigma_pvstar = 0.054 (fit as an AR(1) to
the terms of trade), the world copper price shock (eps_pc alone moves GDP 2.1%
on impact), and the monetary and goods-reallocation shocks.

WHAT IT DOES
------------
Shocks are independent and the model is linear, so the covariance matrix is
additive across shock groups:  Gamma = sum_g Gamma_g. Solving the model with only
group g active therefore gives that group's exact contribution, and the variance
shares sum to 100% — which the script checks.

For each group it reports:
  * share of var(Y_i) for every sector, and of var(GDP)
  * rbar CONDITIONAL on that group alone — "if only these shocks existed, how
    much would sectors comove?"

rbar is a correlation, so it does NOT decompose additively. Reading it per group
is still exactly the diagnostic wanted: a group with a small variance share but
rbar near 1 is a common factor, and if its share is non-trivial it will drag the
total up.

The mechanism is the same one-line trick Stage 2 used: baseline.active_exo_idx
is read in exactly ONE place, to build Sigma_e. Hand the moment function a
baseline whose active_exo_idx lists only one group and everything else is
silenced, with T and R untouched — one model solve, one Lyapunov pass per group.

HOW TO RUN
    julia --threads=1 --project=. variance_decomposition.jl
    VD_TAG=epsY0.8 julia --threads=1 --project=. variance_decomposition.jl

VD_TAG picks which estimate to decompose: VD_TAG=epsY0.8 reads
estimation_results/smm_estimates_epsY0.8.csv. Default is the untagged file.
Use the best VALID sweep point — epsY = 0.8, obj 8.406.

OUTPUT
    tables/variance_decomposition[_TAG].csv
"""

using CSV, DataFrames, Printf, Dynare, Statistics, LinearAlgebra

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

# Include order must match run_smm_estimation.jl — the canonical
# smm_model_moments lives in smm_estimation.jl and must come last.
include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))

const VD_TAG = get(ENV, "VD_TAG", "")


"""
    load_theta_hat(context) -> Vector{Float64}

Read the full 38-element θ from a saved estimates CSV, by NAME. Falls back to
`default_theta0` for anything the file does not carry.
"""
function load_theta_hat(context)
    θ = default_theta0(context)
    fname = isempty(VD_TAG) ? "smm_estimates.csv" : "smm_estimates_$(VD_TAG).csv"
    path  = joinpath(SCRIPT_DIR, "estimation_results", fname)
    isfile(path) || error("""
        $(fname) not found in estimation_results/.
        Run the estimation first, or set VD_TAG to a sweep point that exists,
        e.g. VD_TAG=epsY0.8.
        """)
    est = CSV.read(path, DataFrame)
    d   = Dict(String.(est.param) .=> Float64.(est.value))
    n   = 0
    for (k, name) in enumerate(CSV_PARAM_NAMES)
        if haskey(d, name) && isfinite(d[name])
            θ[k] = d[name]; n += 1
        end
    end
    @printf "  θ from %s (%d of %d entries matched by name)\n" fname n length(CSV_PARAM_NAMES)
    return θ
end


function _main()
    @printf "\n%s\n  VARIANCE DECOMPOSITION BY SHOCK GROUP\n%s\n\n" repeat("=",64) repeat("=",64)

    context = _DYNARE_CONTEXT[]
    context === nothing && error("Dynare context not built — see the @dynare block at the bottom.")

    endo_file  = joinpath(SCRIPT_DIR, "mod", "dynare_endo_names.csv")
    endo_names = isfile(endo_file) ? String.(CSV.read(endo_file, DataFrame).variable) :
                                     Dynare.get_endogenous(context.symboltable)
    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)
    nsec     = baseline.nsec

    θ = load_theta_hat(context)
    SMM_PIN[] = copy(θ)     # the objective reads pinned entries from here
    @printf "  epsY = %.3f  epsM = %.3f  lambda_A = %.3f  lambda_om = %.3f\n\n" θ[2] θ[3] θ[37] θ[38]

    # ---- Shock groups ------------------------------------------------------ #
    exo_names = smm_exo_names(context)
    exo_pos   = Dict(nm => k for (k, nm) in enumerate(exo_names))
    want(nms) = [exo_pos[nm] for nm in nms if haskey(exo_pos, nm)]

    groups = [
        ("sectoral TFP",      want(["epsA_$(i)"   for i in 1:nsec])),
        ("sectoral demand",   want(["eps_om_$(i)" for i in 1:nsec])),
        ("external + copper", want(["eps_pvstar", "eps_pc"])),
        ("monetary",          want(["eps_i"])),
        ("preference (zeta)", want(["eps_zeta"])),
        ("goods realloc.",    want(["eps_omg"])),
    ]

    # Every active shock must land in exactly one group, or the shares will not
    # sum to the total and the decomposition is silently wrong.
    allg = vcat([g[2] for g in groups]...)
    let act = Set(baseline.active_exo_idx), cov = Set(allg)
        length(allg) == length(cov) || error("a shock appears in more than one group")
        if act != cov
            miss  = sort(collect(setdiff(act, cov)))
            extra = sort(collect(setdiff(cov, act)))
            error("""
                Shock groups do not tile the active set.
                In active_exo_idx but no group: $(isempty(miss) ? "none" : [exo_names[i] for i in miss])
                In a group but not active:      $(isempty(extra) ? "none" : [exo_names[i] for i in extra])
                Fix the group list above (or active_shock_indices in smm_model_moments.jl).
                """)
        end
    end
    @printf "  %d active shocks tile %d groups.\n\n" length(allg) length(groups)

    # ---- Total, then each group ------------------------------------------- #
    m_tot, ok = smm_model_moments(θ, context, baseline, endo_names)
    (ok && !any(isnan, m_tot[1:nsec])) || error("model failed at θ̂ with all shocks on")

    # Look the moment indices up by NAME. Hardcoding 78 and 37 would survive
    # today and break silently the next time a moment is inserted — which has
    # already happened twice in this codebase.
    _midx(nm) = let k = findfirst(==(nm), MOMENT_NAMES)
        k === nothing && error("moment \"$(nm)\" not found in MOMENT_NAMES — update this script")
        k
    end
    I_RBAR = _midx("rbar: avg pairwise corr(Y_i,Y_j)")
    I_GDP  = _midx("std(GDP)")
    @printf "  moment indices: rbar=%d  std(GDP)=%d  (of %d)\n\n" I_RBAR I_GDP N_MOMENTS
    var_tot_Y   = m_tot[1:nsec] .^ 2
    var_tot_GDP = m_tot[I_GDP]^2

    ng    = length(groups)
    varY  = zeros(nsec, ng)
    varG  = zeros(ng)
    rbarg = fill(NaN, ng)

    for (gi, (gname, idx)) in enumerate(groups)
        b_g = merge(baseline, (active_exo_idx = idx,))
        m_g, ok_g = smm_model_moments(θ, context, b_g, endo_names)
        if !ok_g
            @printf "  WARNING: model failed with only '%s' active — column left at zero.\n" gname
            continue
        end
        varY[:, gi] .= m_g[1:nsec] .^ 2
        varG[gi]     = m_g[I_GDP]^2
        rbarg[gi]    = m_g[I_RBAR]
    end

    # Additivity check: independent shocks + linear model ⇒ variances sum.
    err_Y = maximum(abs.(sum(varY, dims=2)[:] .- var_tot_Y) ./ max.(var_tot_Y, 1e-30))
    err_G = abs(sum(varG) - var_tot_GDP) / max(var_tot_GDP, 1e-30)
    @printf "  Additivity check: max relative error %.2e (sectors), %.2e (GDP)\n" err_Y err_G
    if max(err_Y, err_G) > 1e-6
        @printf "  %s\n" repeat("!", 70)
        @printf "  Variances do NOT sum to the total. Either a shock is missing from the\n"
        @printf "  groups, or Sigma_e is not being rebuilt per group. Shares below are\n"
        @printf "  NOT trustworthy.\n"
        @printf "  %s\n" repeat("!", 70)
    end
    println()

    # ---- Report ------------------------------------------------------------ #
    names12 = ["Agriculture","Mining","Manufacturing","Utilities","Construction",
               "Trade/Hotels","Transport/Comm","Finance","Real Estate",
               "Business Serv.","Personal Serv.","Public Admin."][1:nsec]

    @printf "  SHARE OF var(Y_i) BY SHOCK GROUP (%%)\n\n"
    @printf "  %-16s" "sector"
    for (gname, _) in groups; @printf " %18s" first(gname, 17); end
    println()
    @printf "  %s\n" repeat("-", 16 + 18*ng)
    for i in 1:nsec
        @printf "  %-16s" names12[i]
        tot = max(sum(varY[i, :]), 1e-30)
        for gi in 1:ng; @printf " %17.1f " (100*varY[i,gi]/tot); end
        println()
    end
    @printf "  %s\n" repeat("-", 16 + 18*ng)
    @printf "  %-16s" "GDP"
    for gi in 1:ng; @printf " %17.1f " (100*varG[gi]/max(sum(varG),1e-30)); end
    println()

    @printf "\n\n  rbar CONDITIONAL ON EACH GROUP ALONE   (data 0.2210, model all-shocks %.4f)\n\n" m_tot[I_RBAR]
    @printf "  %-20s %10s %12s\n" "group" "rbar" "% var(GDP)"
    @printf "  %s\n" repeat("-", 46)
    for (gi, (gname, _)) in enumerate(groups)
        @printf "  %-20s %10.4f %11.1f%%\n" gname rbarg[gi] (100*varG[gi]/max(sum(varG),1e-30))
    end
    @printf "  %s\n" repeat("-", 46)
    @printf "  A group with HIGH rbar and a non-trivial variance share is a common\n"
    @printf "  factor and is what drags the total comovement above the data.\n"

    # ---- COMMON-FACTOR TEST ------------------------------------------------ #
    #
    # WHY. Two sweeps have now failed to move rbar: epsY over 0.7-1.4 (0.416 to
    # 0.457) and epsM over 0.05-0.50 (0.387 to 0.448), against data 0.221. Both
    # are input-substitution elasticities, so if comovement came from cost
    # propagation through Gamma, one of them would have bitten. Neither did.
    #
    # THE HYPOTHESIS. Final demand in NK_SOE_lev_gap2.mod is Cobb-Douglas at both
    # tiers — line 275 for the goods/services composite, lines 334-335 for the
    # sector split:
    #
    #     Cg_j = gammag_j*exp(om_j)/norm_g * (p_g/P_j) * C_g
    #
    # Unit elasticity, so each sector's NOMINAL expenditure is a fixed share of
    # total nominal spending and, in logs,
    #
    #     log C_j = log share_j + log(P*C) - log P_j
    #
    # The term log(P*C) loads with coefficient ONE on all twelve sectors and does
    # not pass through Gamma at all. If that is what sets rbar, no substitution
    # elasticity could ever have touched it, and the om_i shocks — the only
    # objects in the model that move the shares — should be the only force
    # pushing rbar down. The decomposition above already shows exactly that
    # (sectoral demand: rbar = -0.028).
    #
    # THE TEST. Take the model's own HP-filtered covariance matrix and ask how
    # one-dimensional the sectoral output block is:
    #   * PC1 share of the 12x12 sectoral output CORRELATION matrix, against the
    #     33% the same calculation gives on Chilean data;
    #   * each sector's loading on that first component;
    #   * corr(Y_i, C) — the direct correlation with aggregate consumption.
    # A PC1 share near 1 with uniformly high loadings confirms a single common
    # factor with near-unit loadings, i.e. the budget constraint, not the network.
    #
    # Costs one extra model solve. MOMENT_DIAG_SINK must be disarmed afterwards.
    @printf "\n\n  %s\n  COMMON-FACTOR TEST\n  %s\n\n" repeat("=", 62) repeat("=", 62)

    MOMENT_DIAG_SINK[] = true
    _, ok_diag = smm_model_moments(θ, context, baseline, endo_names)
    diag_out = MOMENT_DIAG_SINK[]
    MOMENT_DIAG_SINK[] = nothing        # ALWAYS disarm

    if !ok_diag || !(diag_out isa NamedTuple)
        @printf "  Model failed on the diagnostic solve — test skipped.\n"
    else
        Γd = diag_out.Gamma
        nm = diag_out.names
        pos = Dict(n => k for (k, n) in enumerate(nm))

        yidx = [pos["Y_$(i)"] for i in 1:nsec]
        sd   = [sqrt(max(Γd[k, k], 0.0)) for k in yidx]
        Rm   = [ (sd[a] * sd[b] < 1e-15) ? 0.0 :
                 clamp(Γd[yidx[a], yidx[b]] / (sd[a] * sd[b]), -1.0, 1.0)
                 for a in 1:nsec, b in 1:nsec ]

        ev     = eigen(Symmetric(Rm))
        ord    = sortperm(ev.values, rev = true)
        lam    = ev.values[ord]
        pc1shr = lam[1] / sum(lam)
        # Loading of sector i on PC1, scaled so it reads as a correlation.
        load1  = ev.vectors[:, ord[1]] .* sqrt(max(lam[1], 0.0))
        load1 .*= sign(sum(load1))      # sign of an eigenvector is arbitrary

        kC = get(pos, "C", 0)
        corrC = [ kC == 0 ? NaN :
                  let d = sd[i] * sqrt(max(Γd[kC, kC], 0.0))
                      d < 1e-15 ? 0.0 : clamp(Γd[yidx[i], kC] / d, -1.0, 1.0)
                  end for i in 1:nsec ]

        # BENCHMARK. For an equicorrelated n x n matrix with common correlation
        # r, the first eigenvalue is 1 + (n-1)r, so PC1's share is exactly
        # (1 + (n-1)r)/n. Evaluating that at the DATA's rbar and at the MODEL's
        # rbar gives the two reference points without importing a number from
        # the data script — and the gap between the model's ACTUAL PC1 share and
        # its own equicorrelation benchmark measures how far from
        # one-dimensional the model's comovement really is. A single common
        # factor with near-equal loadings sits right on the benchmark.
        eqpc1(r) = (1 + (nsec - 1) * r) / nsec
        d_rbar   = 0.221
        @printf "  PC1 share of the 12x12 sectoral output correlation matrix\n\n"
        @printf "    model, actual                       %.3f\n" pc1shr
        @printf "    model, equicorrelation at rbar=%.3f  %.3f\n" m_tot[I_RBAR] eqpc1(m_tot[I_RBAR])
        @printf "    DATA,  equicorrelation at rbar=%.3f  %.3f\n\n" d_rbar eqpc1(d_rbar)
        @printf "    (the data's own PC1 share is printed by compute_data_moments.jl;\n"
        @printf "     do not reuse the 33%% figure — that one is the 24-column Y-and-P\n"
        @printf "     panel, not this 12-column output block.)\n\n"
        @printf "  %-16s %12s %14s\n" "sector" "PC1 loading" "corr(Y_i, C)"
        @printf "  %s\n" repeat("-", 44)
        for i in 1:nsec
            @printf "  %-16s %12.3f %14.3f\n" names12[i] load1[i] corrC[i]
        end
        @printf "  %s\n" repeat("-", 44)
        @printf "  %-16s %12.3f %14.3f\n" "mean |loading|" mean(abs.(load1)) mean(filter(isfinite, abs.(corrC)))
        @printf "\n  eigenvalue shares: %s\n" join([@sprintf("%.2f", l/sum(lam)) for l in lam[1:min(4,end)]], "  ")
        @printf "\n  READING IT.\n"
        @printf "  CONFIRMS the Cobb-Douglas-demand diagnosis: PC1 share close to its own\n"
        @printf "    equicorrelation benchmark, loadings uniformly high and similar across\n"
        @printf "    sectors, corr(Y_i,C) high nearly everywhere. One factor with near-unit\n"
        @printf "    loadings is the budget constraint, and no substitution elasticity can\n"
        @printf "    reach it — which is exactly why the epsY and epsM sweeps both failed.\n"
        @printf "  REFUTES it: loadings spread out or changing sign, low corr(Y_i,C), or a\n"
        @printf "    PC1 share well BELOW the benchmark. Then comovement is many-dimensional,\n"
        @printf "    the demand block is not the culprit, and this diagnosis is wrong.\n"
        @printf "  Mining is the sector to watch either way — the decomposition puts it at\n"
        @printf "    100%% own TFP, so it is the natural dissenter from any common factor.\n"

        mkpath(joinpath(SCRIPT_DIR, "tables"))
        CSV.write(joinpath(SCRIPT_DIR, "tables",
                           "common_factor_test$(isempty(VD_TAG) ? "" : "_" * VD_TAG).csv"),
                  DataFrame(sector = names12, pc1_loading = load1, corr_Y_C = corrC))
    end

    # ---- Save -------------------------------------------------------------- #
    out = DataFrame(sector = [names12; "GDP"])
    for (gi, (gname, _)) in enumerate(groups)
        col = [100 .* varY[:, gi] ./ max.(sum(varY, dims=2)[:], 1e-30);
               100 * varG[gi] / max(sum(varG), 1e-30)]
        out[!, Symbol(replace(gname, " " => "_", "." => "", "+" => "and"))] = col
    end
    rb = DataFrame(group = [g[1] for g in groups], rbar_alone = rbarg,
                   share_var_GDP_pct = 100 .* varG ./ max(sum(varG), 1e-30))

    mkpath(joinpath(SCRIPT_DIR, "tables"))
    sfx = isempty(VD_TAG) ? "" : "_" * VD_TAG
    f1 = joinpath(SCRIPT_DIR, "tables", "variance_decomposition$(sfx).csv")
    f2 = joinpath(SCRIPT_DIR, "tables", "variance_decomposition_rbar$(sfx).csv")
    CSV.write(f1, out); CSV.write(f2, rb)
    @printf "\n  -> %s\n  -> %s\n\n" f1 f2
end


# =========================================================================== #
#  BUILD THE DYNARE CONTEXT AT TOP LEVEL                                      #
# =========================================================================== #
# Same statement-by-statement layout as run_smm_estimation.jl: @dynare expands
# at MACRO-EXPANSION time so it cannot sit inside a let/begin block, and a
# deserialised context cannot be solved because Dynare.DFunctions is only
# populated by an in-process @dynare.
const _DYNARE_CONTEXT = Ref{Any}(nothing)

_VD_MOD_DIR = joinpath(SCRIPT_DIR, "mod")

isfile(joinpath(_VD_MOD_DIR, "params_jl.mod")) ||
    error("mod/params_jl.mod not found. Run main_SOE_gap.jl first.")

_VD_LOG     = joinpath(_VD_MOD_DIR, "dynare_vd_load.log")
_VD_OLD_PWD = pwd()

@printf "\n--- Loading model via @dynare ---\n    output → %s\n" basename(_VD_LOG)

cd(_VD_MOD_DIR)

_VD_STDOUT = stdout
_VD_STDERR = stderr
_VD_IO     = open(_VD_LOG, "w")
redirect_stdout(_VD_IO)
redirect_stderr(_VD_IO)

_VD_CTX = try
    @dynare "NK_SOE_lev_gap2"
finally
    redirect_stdout(_VD_STDOUT)
    redirect_stderr(_VD_STDERR)
    close(_VD_IO)
    cd(_VD_OLD_PWD)
end

_VD_CTX isa Dynare.Context ||
    error("@dynare returned $(typeof(_VD_CTX)), not a Context — see $(basename(_VD_LOG)).")
_DYNARE_CONTEXT[] = _VD_CTX

@printf "    done — %d endogenous, %d exogenous\n" length(
    _VD_CTX.results.model_results[1].trends.endogenous_steady_state) _VD_CTX.models[1].exogenous_nbr

Base.invokelatest(_main)
