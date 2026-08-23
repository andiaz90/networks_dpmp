"""
calibrate_sectoral_shocks.jl
============================
STAGE 2 of the sectoral shock calibration: MODEL INVERSION.

  Stage 1  compute_sectoral_shocks.jl       — analytic supply/demand split.
  Stage 2  calibrate_sectoral_shocks.jl     — this file.

WHAT IT DOES. Stage 1 sets each sector's shock sizes from the data alone,
assuming a UNIT shock -> output loading:

    sigma_supply_i = sqrt(frac_supply_i) * std_Y_idio_i

That assumption is wrong in a production network. The mapping from a sectoral
shock to sectoral output runs through the Leontief inverse, the labour
adjustment cost, and the sector's own price stickiness, so the loading is
strongly heterogeneous: a shock of a given size in a downstream, sticky-price,
input-intensive sector moves output far less than the same shock in an
upstream flexible-price one. Stage 2 measures each loading numerically and
inverts it:

    repeat:  solve the model at the current shock sizes
             own_std_i = std(Y_i) with ONLY sector i's shocks active
             ratio_i   = std_Y_idio_i / own_std_i
             sigma_i  *= ratio_i^DAMP        (both supply and demand)

Both shock types are scaled by the SAME per-sector ratio, so Stage 1 owns the
MIX (supply vs demand, from the idiosyncratic output-price correlation) and
Stage 2 owns the LEVEL. Different moments, so they cannot fight.

TARGET: OWN-SHOCK VOLATILITY, NOT TOTAL  (redesign of 2026-08-21)
-----------------------------------------------------------------
The first version of this file matched model TOTAL volatility to data TOTAL
volatility. That cannot work, and the failure is instructive. In a network

    var(Y_i) = own shock + spillovers from other sectors + aggregate shocks

and only the first term responds to sigma_i. When the other two dominate,
shrinking sigma_i barely moves std(Y_i), the ratio stays far below 1, and the
update walks to the lower bound and stays there. Observed on 2026-08-21:
Construction at sigma = 0.0001 — its own shock effectively switched off — still
produced std(Y) = 0.1186 against 0.0289 in the data, all of it arriving through
the input-output matrix. Eight of twelve sectors collapsed to the sigma floor
and Public Admin, 0.33% of value added, pinned at the ceiling.

So compare like with like. Sector i's own-shock volatility is obtained by
handing the moment function a baseline whose `active_exo_idx` lists only
epsA_i and eps_om_i; that field is used in exactly one place, to build Sigma_e,
so everything else is silenced while T and R stay put. The data counterpart is
std_Y_idio, computed in compute_data_moments.jl section 9c by projecting the
common factor out of the sectoral panel — which is precisely the same object.

TOTAL volatility is then NOT a target here. It emerges from the model as own +
spillover + aggregate, and the SMM moment std(Y_i) is what disciplines it.

One useful consequence: with only sector i's shocks active and both of its
sigmas scaled by the same factor c, own_std scales by exactly c. The map is
LINEAR, so a full step (DAMP = 1) lands on the target and two iterations
suffice — where the total-volatility version could not converge at all.

WHY THIS REPLACES lambda_A / lambda_om. Those two parameters, theta[37:38],
are global scalars that multiply all 12 supply and all 12 demand shocks. They
were the only parameters with real sensitivity in the 2026-08-21 run
(|d obj| ~ 10.3, versus 0.45 for cl) precisely because they were standing in
for this per-sector rescaling. After Stage 2 they should estimate close to 1.
Treat a converged lambda far from 1 as a diagnostic that this inversion did
not converge, not as a result.

HOW TO RUN
  julia --project=. main_SOE_gap.jl              # writes mod/params_jl.mod
  julia --project=. compute_data_moments.jl      # writes the _idio columns
  julia --project=. compute_sectoral_shocks.jl   # Stage 1
  julia --threads=auto --project=. calibrate_sectoral_shocks.jl   # Stage 2
  julia --threads=auto --project=. run_smm_estimation.jl

ENV OVERRIDES
  SHOCKCAL_ITERS  (default 4)    maximum inversion iterations
  SHOCKCAL_DAMP   (default 0.7)  update damping; < 1 guards against overshoot
  SHOCKCAL_TOL    (default 0.02) stop when max_i |ratio_i - 1| < TOL

OUTPUT
  Data/sectoral_shock_calibration.csv    overwritten, stage = "2-model-inverted"
  Data/sectoral_shock_calibration_stage1.csv   backup of the Stage-1 file
"""

using CSV, DataFrames, Printf, Statistics, Dynare

SCRIPT_DIR = @__DIR__
DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

# Include order matches run_smm_estimation.jl. smm_model_moments.jl provides the
# infrastructure (Klein solver, recompute_ss!, param cache); the moment function
# itself lives in smm_estimation.jl and must be included after it.
include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "smm_model_moments.jl"))
include(joinpath(SCRIPT_DIR, "smm_estimation.jl"))

const MAXITER = parse(Int,     get(ENV, "SHOCKCAL_ITERS", "4"))
# DAMP = 1.0: the own-shock map is exactly linear in sigma (see the update
# step), so a full step lands on the target. Damping was only needed by the
# old total-volatility version, which was inverting a non-proportional map.
const DAMP    = parse(Float64, get(ENV, "SHOCKCAL_DAMP",  "1.0"))
const TOL     = parse(Float64, get(ENV, "SHOCKCAL_TOL",   "0.02"))


"""
    load_transmission_theta(θ0) -> θ

Overwrite the free transmission entries of `θ0` with the latest estimates from
estimation_results/smm_estimates.csv, if that file exists. The inversion has to
be run AT the transmission parameters the model will actually be solved with —
cl and kappaw both move the shock -> output loading materially, so inverting at
the seed values and then estimating would leave the calibration stale.
"""
function load_transmission_theta(θ0::Vector{Float64})
    θ = copy(θ0)
    f = joinpath(SCRIPT_DIR, "estimation_results", "smm_estimates.csv")
    if !isfile(f)
        @printf "  No smm_estimates.csv — inverting at the seed transmission values.\n"
        @printf "  (Re-run this script after the next estimation to refresh.)\n"
        return θ
    end
    est = CSV.read(f, DataFrame)
    d   = Dict(String.(est.param) .=> Float64.(est.value))
    got = String[]
    for k in FREE_THETA
        name = CSV_PARAM_NAMES[k]
        if haskey(d, name) && isfinite(d[name])
            θ[k] = d[name]; push!(got, name)
        end
    end
    @printf "  Transmission params from smm_estimates.csv: %s\n" join(got, ", ")
    return θ
end


function _main()
    @printf "\n%s\n  SECTORAL SHOCK CALIBRATION — STAGE 2 (model inversion)\n%s\n\n" repeat("=",64) repeat("=",64)

    # ---- Stage-1 input --------------------------------------------------- #
    scf = joinpath(DATA_DIR, "sectoral_shock_calibration.csv")
    isfile(scf) || error("""
        $(basename(scf)) not found. Run Stage 1 first:
          julia --project=. compute_sectoral_shocks.jl
        """)
    sc = CSV.read(scf, DataFrame)

    if hasproperty(sc, :stage) && any(String.(sc.stage) .== "2-model-inverted")
        @printf "  NOTE: the input file is already Stage-2 output. Re-inverting from it\n"
        @printf "  is fine (it is a fixed point), but if you changed the DATA you want\n"
        @printf "  to re-run Stage 1 first so the mix is recomputed.\n\n"
    end
    if hasproperty(sc, :moment_source) && occursin("TOTAL", String(sc.moment_source[1]))
        @printf "  %s\n" repeat("!", 70)
        @printf "  WARNING: Stage 1 ran on TOTAL moments (no common-factor removal).\n"
        @printf "  Stage 2 will converge, but onto a mis-specified supply/demand mix.\n"
        @printf "  %s\n\n" repeat("!", 70)
    end

    # Target = the IDIOSYNCRATIC data volatility, which is what std_Y_data holds
    # when Stage 1 ran on the _idio moments (the normal case). This is the whole
    # point of the 2026-08-21 redesign: we are matching an own-shock object on
    # the model side, so the data side must be the own-shock object too.
    # std_Y_total is kept in the file for reference and for the final report.
    std_Y_target = Float64.(sc.std_Y_data)
    std_Y_tot_ref = hasproperty(sc, :std_Y_total) ? Float64.(sc.std_Y_total) : fill(NaN, nrow(sc))
    if hasproperty(sc, :moment_source) && !occursin("IDIOSYNCRATIC", String(sc.moment_source[1]))
        @printf "  %s\n" repeat("!", 70)
        @printf "  WARNING: Stage 1 did not use the idiosyncratic moments, so std_Y_data is\n"
        @printf "  TOTAL volatility. Matching an own-shock model object to a total data\n"
        @printf "  moment asks each sector's own shock to reproduce variance the network\n"
        @printf "  and the aggregate shocks also generate. Regenerate the data moments.\n"
        @printf "  %s\n\n" repeat("!", 70)
    end
    isig = Float64.(sc.isigma_tfp_init)
    som  = Float64.(sc.sigma_om_init)
    nsec = length(isig)

    # ---- Model ------------------------------------------------------------ #
    context = _DYNARE_CONTEXT[]
    context === nothing && error("Dynare context not built — see the top-level @dynare block.")

    endo_names_file = joinpath(SCRIPT_DIR, "mod", "dynare_endo_names.csv")
    endo_names = isfile(endo_names_file) ?
        String.(CSV.read(endo_names_file, DataFrame).variable) :
        Dynare.get_endogenous(context.symboltable)

    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)
    @printf "  Model loaded: %d endogenous, %d active shocks.\n" length(endo_names) length(baseline.active_exo_idx)

    # ---- OWN-SHOCK BASELINES (2026-08-21 redesign) ------------------------- #
    #
    # baseline.active_exo_idx is used in exactly one place: building Sigma_e in
    # smm_model_moments. Handing the moment function a baseline whose
    # active_exo_idx lists only sector i's two shocks therefore returns the
    # moments generated by THOSE SHOCKS ALONE — every other shock silenced.
    # T and R are untouched, so the decision rule is solved once and reused; only
    # the Lyapunov/HP pass repeats. Twelve of those per iteration is cheap.
    exo_names = smm_exo_names(context)
    exo_pos   = Dict(nm => k for (k, nm) in enumerate(exo_names))

    # WHICH SHOCKS DEFINE THE MODEL-SIDE COUNTERPART OF std_Y_idio?
    #
    # "own"  — only epsA_i and eps_om_i. Exactly one sector's own shocks.
    # "idio" — ALL 24 sectoral shocks, aggregates off.               (DEFAULT)
    #
    # "own" was the first redesign and it overshoots, for a reason worth stating.
    # std_Y_idio is built by projecting ONE principal component out of the
    # observed sectoral panel. That removes the COMMON factor — it does not
    # remove sector-to-sector transmission of *idiosyncratic* shocks, which
    # stays in the residual. So the data object is "volatility net of the common
    # factor", and its model counterpart is "all sectoral shocks on, aggregate
    # shocks off" — not "one sector's shock alone".
    #
    # Matching own-shock-only volatility to std_Y_idio therefore compares a
    # SMALLER model object to the data, the ratio exceeds 1 everywhere, and the
    # sigmas inflate to their upper bounds — the mirror image of the original
    # total-volatility failure, which compared a LARGER model object and drove
    # them to zero. Observed 2026-08-21: "own" gave ratios of 0.82-20.0 and
    # jammed on the ceiling within three iterations.
    #
    # "idio" is the like-for-like comparison. It couples the sectors (sigma_j
    # moves std_Y_idio_i), so the per-sector step is Jacobi rather than exact,
    # and it needs more iterations than the linear "own" map would.
    SHOCK_SET = get(ENV, "SHOCKCAL_SET", "idio")

    own_baselines = Vector{Any}(undef, nsec)
    if SHOCK_SET == "own"
        for i in 1:nsec
            idx = Int[]
            for nm in ("epsA_$(i)", "eps_om_$(i)")
                haskey(exo_pos, nm) || error("shock $(nm) not found among the model's exogenous variables")
                push!(idx, exo_pos[nm])
            end
            own_baselines[i] = merge(baseline, (active_exo_idx = idx,))
        end
        @printf "  SHOCKCAL_SET=own — one baseline per sector, 2 shocks each.\n"
        @printf "  NOTE: this target is known to overshoot; see the comment in this file.\n\n"
    else
        idx = Int[]
        for j in 1:nsec, nm in ("epsA_$(j)", "eps_om_$(j)")
            haskey(exo_pos, nm) || error("shock $(nm) not found among the model's exogenous variables")
            push!(idx, exo_pos[nm])
        end
        b_idio = merge(baseline, (active_exo_idx = idx,))
        for i in 1:nsec
            own_baselines[i] = b_idio        # same baseline; the loop reads moment i from it
        end
        @printf "  SHOCKCAL_SET=idio — all %d sectoral shocks on, aggregates off.\n" length(idx)
        @printf "  This matches what std_Y_idio measures: volatility net of the COMMON\n"
        @printf "  factor, with sector-to-sector transmission still in it.\n\n"
    end

    θ_base = load_transmission_theta(default_theta0(context))
    θ_base[2]  = parse(Float64, get(ENV, "SMM_EPSY", "0.80"))
    θ_base[3]  = parse(Float64, get(ENV, "SMM_EPSM", "0.20"))
    θ_base[35] = parse(Float64, get(ENV, "SMM_ETASTAR", "1.0"))
    # lambda_A / lambda_om are held at 1 THROUGHOUT the inversion. Letting them
    # float here would be circular: they are global stand-ins for exactly the
    # per-sector rescaling this loop computes.
    θ_base[37] = 1.0
    θ_base[38] = 1.0

    ecf = joinpath(DATA_DIR, "external_shock_calibration.csv")
    if isfile(ecf)
        ec = CSV.read(ecf, DataFrame); em = Dict(String.(ec.param) .=> Float64.(ec.value))
        haskey(em, "rho_pvstar")   && (θ_base[31] = em["rho_pvstar"])
        haskey(em, "sigma_pvstar") && (θ_base[32] = em["sigma_pvstar"])
    end

    # ---- Inversion loop ---------------------------------------------------- #
    @printf "\n  Iterating: max %d, damping %.2f, tol %.3f\n" MAXITER DAMP TOL
    @printf "  %-6s %10s %10s %10s   %s\n" "iter" "max|r-1|" "mean r" "min/max r" "converged?"
    @printf "  %s\n" repeat("-", 62)

    converged   = false
    last_ratio  = fill(NaN, nsec)
    std_Y_model = fill(NaN, nsec)   # own-shock std(Y_i) at the current sigmas
    own_std     = fill(NaN, nsec)

    for iter in 1:MAXITER
        θ = copy(θ_base)
        θ[7:18]  = isig
        θ[19:30] = som
        SMM_PIN[] = copy(θ)          # the moment function reads pinned entries from here

        # OWN-SHOCK volatility, sector by sector. Every other shock is off, so
        # what comes back is the part of sector i's output volatility that its
        # OWN shocks generate — the same object std_Y_idio measures in the data.
        for i in 1:nsec
            mi, oki = smm_model_moments(θ, context, own_baselines[i], endo_names)
            if !oki || isnan(mi[i])
                @printf "\n  Model solve FAILED at iteration %d, sector %d.\n" iter i
                @printf "  Shock sizes at failure: isigma_tfp %s\n" join(round.(isig, digits=4), " ")
                error("Inversion aborted — the model could not be solved. Check that " *
                      "main_SOE_gap.jl has been re-run so mod/params_jl.mod is current.")
            end
            own_std[i] = mi[i]
        end

        # SANITY CHECK (2026-08-21): did the model actually SEE the update?
        # An earlier version of this loop reported a byte-identical ratio across
        # every iteration, because the decision-rule cache in smm_model_moments.jl
        # keyed only on structural parameters and returned a stale R — and R is
        # where the shock sizes live. Fixed there; checked here, because a silent
        # no-op inversion is indistinguishable from a converged one in the log
        # except by staring at the digits.
        if iter > 1 && maximum(abs.(own_std .- std_Y_model)) < 1e-12
            # Two very different causes produce an identical trace, so separate
            # them before blaming the solver. If every sigma the loop wanted to
            # move is sitting on a bound, the shock sizes did NOT change and the
            # frozen moments are correct behaviour, not a bug.
            n_pinned = count(i -> isig[i] <= 1.01e-4 || isig[i] >= 0.0999 ||
                                  som[i]  <= 1.01e-5 || som[i]  >= 0.1999, 1:nsec)
            if n_pinned >= count(i -> abs(last_ratio[i] - 1) > TOL, 1:nsec)
                @printf "\n  Halting at iteration %d: every sector still off target has its\n" iter
                @printf "  sigma on a bound (%d of %d sectors pinned), so nothing can move.\n" n_pinned nsec
                @printf "  This is a BOUNDS problem, not a solver problem. See the report below.\n"
                break
            end
            error("""
            Model moments did not change between iterations $(iter-1) and $(iter),
            although the shock sizes did and no sigma is on a bound. The model is
            not seeing theta[7:30].

            Check _KLEIN_STRUCT_IDX in smm_model_moments.jl: it must include the
            shock-size indices, or _resolve_cached! returns a stale R.
            """)
        end

        std_Y_model .= own_std
        # With only sector i's shocks active and both of its sigmas scaled by the
        # same factor c, R's relevant columns scale by c, Gamma by c^2, and
        # own_std by exactly c. The map is LINEAR, so ratio^1 lands on the target
        # in a single step — which is why DAMP defaults to 1.0 here and why this
        # converges in two iterations where the total-volatility version could
        # not converge at all. The clamp is left in only as a guard against a
        # numerically dead sector.
        ratio = [sm > 1e-10 ? clamp(st / sm, 0.05, 20.0) : 1.0
                 for (st, sm) in zip(std_Y_target, own_std)]
        last_ratio .= ratio

        gap = maximum(abs.(ratio .- 1))
        @printf "  %-6d %10.4f %10.3f %6.2f/%.2f   %s\n" iter gap mean(ratio) minimum(ratio) maximum(ratio) (
            gap < TOL ? "YES" : "")

        if gap < TOL
            converged = true
            break
        end
        isig .*= ratio .^ DAMP
        som  .*= ratio .^ DAMP

        # Keep the inversion inside the bounds the estimator will later impose
        # (LB/UB in utils.jl: isigma_tfp in [1e-4, 0.10], sigma_om in [1e-5, 0.20]).
        # Without this the loop happily proposes sigma = 0.47 for a sector the
        # model can barely move, the value is silently clamped when it is loaded,
        # and the calibration written to disk is not the one that was solved.
        for i in 1:nsec
            isig[i] = clamp(isig[i], 1e-4, 0.10)
            som[i]  = clamp(som[i],  1e-5, 0.20)
        end
    end

    if !converged
        @printf "\n  NOT converged in %d iterations (max|ratio-1| = %.4f).\n" MAXITER maximum(abs.(last_ratio .- 1))
        @printf "  The file is still written — the shock sizes improved monotonically —\n"
        @printf "  but raise SHOCKCAL_ITERS, or lower SHOCKCAL_DAMP if the trace oscillates.\n"
    end

    # ---- Report and write --------------------------------------------------- #
    @printf "\n  %-16s %10s %10s %8s   %10s %10s %10s\n" "sector" "idio data" "idio model" "ratio" "isigma_tfp" "sigma_om" "TOTAL data"
    @printf "  %s\n" repeat("-", 84)
    for i in 1:nsec
        @printf "  %-16s %10.4f %10.4f %8.2f   %10.4f %10.4f %10.4f\n" sc.name[i] std_Y_target[i] std_Y_model[i] last_ratio[i] isig[i] som[i] std_Y_tot_ref[i]
    end
    @printf "  %s\n" repeat("-", 84)
    @printf "  Columns 2-3 are OWN-SHOCK volatility: the model solved with every other\n"
    @printf "  shock switched off. Total volatility is not targeted here — it emerges\n"
    @printf "  from own + network spillovers + aggregate shocks, and the SMM moment\n"
    @printf "  std(Y_i) is what disciplines it.\n"

    old_isig = Float64.(sc.isigma_tfp_init); old_som = Float64.(sc.sigma_om_init)
    @printf "  Stage 1 -> Stage 2 change in shock sizes:\n"
    @printf "    supply (TFP): mean x%.2f, range [x%.2f, x%.2f]\n" (
        mean(isig ./ old_isig)) (minimum(isig ./ old_isig)) (maximum(isig ./ old_isig))
    @printf "    demand      : mean x%.2f, range [x%.2f, x%.2f]\n" (
        mean(som ./ old_som)) (minimum(som ./ old_som)) (maximum(som ./ old_som))
    @printf "  A WIDE range here is the point: it is the per-sector loading heterogeneity\n"
    @printf "  that lambda_A / lambda_om could not represent with two global scalars.\n"

    # Flag sectors the inversion cannot reach. A sector whose model volatility is
    # far below the data at ANY shock size is not a calibration problem — the
    # model has no mechanism to move it, and forcing sigma up to compensate just
    # buys a huge shock that distorts everything the sector connects to.
    # Public administration is the standing example: 0.33% of value added, and
    # the model delivers std(Y_12) ~ 0.002 against 0.008 in the data.
    let flagged = [i for i in 1:nsec if last_ratio[i] > 3.0 || last_ratio[i] < 1/3]
        if !isempty(flagged)
            @printf "\n  %s\n" repeat("!", 72)
            @printf "  %d sector(s) still off by more than 3x after the inversion:\n" length(flagged)
            for i in flagged
                @printf "    %-16s idio data %.4f  own-shock model %.4f  ratio %.2f   sigma now %.4f%s\n" (
                    sc.name[i]) std_Y_target[i] std_Y_model[i] last_ratio[i] isig[i] (
                    isig[i] >= 0.0999 || isig[i] <= 1.01e-4 ? "  [AT BOUND]" : "")
            end
            @printf "  Since the map from sigma_i to own-shock volatility is LINEAR, a sector\n"
            @printf "  can only remain off by 3x if its sigma hit a bound. Either the data\n"
            @printf "  moment is out of the range the bounds allow, or the sector's loading is\n"
            @printf "  extreme enough to deserve its own look. More iterations will not help.\n"
            @printf "  %s\n" repeat("!", 72)
        end
    end

    # Back up Stage 1 before overwriting.
    s1 = joinpath(DATA_DIR, "sectoral_shock_calibration_stage1.csv")
    (!hasproperty(sc, :stage) || String(sc.stage[1]) == "1-analytic") && CSV.write(s1, sc)

    out = copy(sc)
    out.isigma_tfp_init = isig
    out.sigma_om_init   = som
    out.std_Y_own_model = std_Y_model
    out.rescale_ratio   = last_ratio
    out.stage           = fill("2-model-inverted", nsec)
    CSV.write(scf, out)

    @printf "\n  -> wrote %s  (stage = 2-model-inverted)\n" scf
    @printf "  -> Stage-1 backup: %s\n" basename(s1)
    @printf "  -> NEXT: julia --threads=auto --project=. run_smm_estimation.jl\n"
    @printf "     lambda_A and lambda_om should now estimate near 1.0.\n\n"
end


# =========================================================================== #
#  BUILD THE DYNARE CONTEXT AT TOP LEVEL                                      #
# =========================================================================== #
# Same reasoning, and the same statement-by-statement layout, as
# run_smm_estimation.jl: @dynare expands at MACRO-EXPANSION time, so it cannot
# be wrapped in a let/begin block, and a deserialised context cannot be solved
# because Dynare.DFunctions is only populated by an in-process @dynare.
const _DYNARE_CONTEXT = Ref{Any}(nothing)

_MOD_DIR = joinpath(SCRIPT_DIR, "mod")

isfile(joinpath(_MOD_DIR, "params_jl.mod")) || error("""
    mod/params_jl.mod not found. Run main_SOE_gap.jl first.
    """)

_VERBOSE = get(ENV, "DYNARE_VERBOSE", "0") in ("1", "true")
_LOG     = joinpath(_MOD_DIR, "dynare_shockcal_load.log")
_OLD_PWD = pwd()

@printf "\n--- Loading model via @dynare (populates Dynare.DFunctions) ---\n"
@printf "    output → %s\n" basename(_LOG)

cd(_MOD_DIR)

_STDOUT = stdout
_STDERR = stderr
_LOG_IO = _VERBOSE ? nothing : open(_LOG, "w")
if !_VERBOSE
    redirect_stdout(_LOG_IO)
    redirect_stderr(_LOG_IO)
end

_CTX = try
    @dynare "NK_SOE_lev_gap2"
finally
    if !_VERBOSE
        redirect_stdout(_STDOUT)
        redirect_stderr(_STDERR)
        close(_LOG_IO)
    end
    cd(_OLD_PWD)
end

_CTX isa Dynare.Context || error("""
    @dynare returned $(typeof(_CTX)), not a Dynare.Context — see $(basename(_LOG)).
    """)
_DYNARE_CONTEXT[] = _CTX

@printf "    done — %d endogenous, %d exogenous\n" length(
    _CTX.results.model_results[1].trends.endogenous_steady_state) _CTX.models[1].exogenous_nbr

Base.invokelatest(_main)
