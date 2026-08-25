"""
smm_estimation.jl
=================
SMM estimation for the NK-IOSOE Chile model.

STRATEGY
  First-order perturbation via Klein (2000) / Dynare.jl.
  For each θ, re-solve the linearised system without recompiling, then
  compute HP-filtered analytical moments via the spectral Lyapunov approach.

FIXES vs. original:
  1. HP filter applied to model moments (matches data HP-filtered log deviations)
  2. Monetary policy shock (eps_i, index 2) activated in Σe
  3. Lag-1 autocovariance computed via spectral approach (spurious term removed)
  4. Tighter, economically motivated bounds (fewer BK failures)
  5. Proportional weighting + 5× on rank correlations
  6. Per-thread SS cache with warm-start tolerance
  7. Pre-allocated scratch arrays, in-place Lyapunov, @views throughout

PREREQUISITES
  Run main_SOE_gap.jl first (EXERCISE=0).

Usage:
  julia --threads=auto --project=. run_smm_estimation.jl
"""

using LinearAlgebra, Statistics, StatsBase, Printf, Dates
using NLsolve, CSV, DataFrames, Dynare
using CMAEvolutionStrategy

# θ-REDUCTION: SMM_PIN holds a full 36-vector whose 28 non-transmission entries PIN
# the objective (elasticities, 24 measured sectoral shock sizes, external rho/sigma).
# Only the 7 transmission params are searched. This is now the ONLY estimation mode.
const SMM_PIN = Ref{Union{Nothing,Vector{Float64}}}(nothing)

# Thrown from the CMA-ES objective when the wall-clock self-limit is hit, so the
# optimiser unwinds and the caller can save the best θ before SLURM SIGKILLs the
# job. Detected by string match so it survives Task/Composite exception wrapping
# under multi-threaded evaluation.
struct SMMTimeout <: Exception end
Base.showerror(io::IO, ::SMMTimeout) = print(io, "SMMTimeout: wall-clock self-limit reached")
# Recognise the timeout whether bare or wrapped by the multi-threaded optimiser.
_is_timeout(err) = occursin("SMMTimeout", sprint(showerror, err))
_is_timeout(err::SMMTimeout) = true
_is_timeout(err::CompositeException) = any(_is_timeout, err.exceptions)
_is_timeout(err::TaskFailedException) = _is_timeout(err.task.exception)

SCRIPT_DIR  = @__DIR__
REPO_ROOT   = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR    = joinpath(REPO_ROOT, "Data")

# Dependency files are included by run_smm_estimation.jl BEFORE this file.

@printf "\n%s\n  SMM ESTIMATION: NK-IOSOE Chile Model\n%s\n\n" repeat("=",60) repeat("=",60)


# =========================================================================== #
#  1.  LOAD DATA MOMENTS                                                       #
# =========================================================================== #

NSEC     = 12
GOODS    = [1,2,3,4,5]
SERVICES = [6,7,8,9,10,11,12]

@printf "--- Loading data moments (CSV) ---\n"

sec_mom_file = joinpath(DATA_DIR, "sectoral_moments.csv")
agg_mom_file = joinpath(DATA_DIR, "aggregate_moments.csv")

if isfile(sec_mom_file) && isfile(agg_mom_file)
    sec = CSV.read(sec_mom_file, DataFrame)
    agg = CSV.read(agg_mom_file, DataFrame)
    agg_dict = Dict(String(r.moment) => Float64(r.value) for r in eachrow(agg))

    y_d          = Float64.(sec.std_Y)
    p_d          = Float64.(sec.std_PH)
    l_d          = Float64.(sec.std_L)
    corr_YPH_d   = hasproperty(sec, :corr_YPH) ? Float64.(sec.corr_YPH) : fill(0.0, NSEC)
    ac_y_d       = hasproperty(sec, :autocorr_Y) ? Float64.(sec.autocorr_Y) : fill(NaN, NSEC)
    d_std_GDP    = agg_dict["std_GDP"]
    d_std_pi     = agg_dict["std_pi"]
    d_corr_GDPpi = agg_dict["corr_GDPpi"]
    d_omG        = agg_dict["omG"]
    d_std_Q      = agg_dict["std_Q"]
    d_autocorr_Q = agg_dict["autocorr_Q"]
    d_corr_GDPQ  = agg_dict["corr_GDPQ"]
    d_TBGDP      = agg_dict["TBGDP"]
    # std(TB/GDP): HP-filtered std dev of trade-balance-to-GDP ratio.
    # Computed by compute_data_moments.jl (section 7b) and saved to aggregate_moments.csv.
    # NO fallback (2026-07-08): stale csv must be regenerated, not papered over.
    haskey(agg_dict, "std_TBGDP") || error(
        "aggregate_moments.csv is STALE: std_TBGDP not found. " *
        "Regenerate with: julia --project=. compute_data_moments.jl")
    d_std_TBGDP  = agg_dict["std_TBGDP"]
    # Employment comovement moments (added 2026-07-08): corr(N,GDP) ≈ +0.70,
    # corr(N,GDP/N) ≈ +0.07 in Chilean data (HP-1600, 2009Q1–2023Q4). They
    # discipline the demand/supply shock mix and kappaw. NO fallback: a stale
    # aggregate_moments.csv (pre-2026-07-08, also missing the sector-8/10
    # valid-window fix) must not silently produce results.
    (haskey(agg_dict, "corr_NGDP") && haskey(agg_dict, "corr_NAPL")) || error("""
        aggregate_moments.csv is STALE: corr_NGDP / corr_NAPL not found.
        It predates the 2026-07-08 moment additions (and the sector-8/10
        valid-window fix). Regenerate it before estimating:
            julia --project=. compute_data_moments.jl
        """)
    haskey(agg_dict, "autocorr_GDP") || error("""
        aggregate_moments.csv is STALE: autocorr_GDP not found.
        The persistence moments (64-77) were added 2026-08-21. Regenerate:
          rm Data/sectoral_moments.csv Data/aggregate_moments.csv
          julia --project=. compute_data_moments.jl
        """)
    d_ac_GDP     = agg_dict["autocorr_GDP"]
    d_ac_pi      = agg_dict["autocorr_pi"]
    haskey(agg_dict, "rbar_YY") || error("""
        aggregate_moments.csv is STALE: rbar_YY not found.
        Moments 78-84 were added 2026-08-21. Regenerate:
          rm Data/sectoral_moments.csv Data/aggregate_moments.csv
          julia --project=. compute_data_moments.jl
          python3 Data/build_reallocation_calibration.py
        """)
    d_rbar_YY    = agg_dict["rbar_YY"]
    d_ratio_stdC = agg_dict["ratio_stdC"]
    d_ratio_stdI = agg_dict["ratio_stdI"]
    d_corr_CGDP  = agg_dict["corr_CGDP"]
    d_corr_IGDP  = agg_dict["corr_IGDP"]
    # 79-80: goods/services average sectoral employment volatility, already in
    # the sectoral CSV (constant down each group's rows).
    d_std_Lg     = hasproperty(sec, :std_Lg) ? first(skipmissing(filter(!isnan, Float64.(sec.std_Lg)))) : NaN
    d_std_Ls     = hasproperty(sec, :std_Ls) ? first(skipmissing(filter(!isnan, Float64.(sec.std_Ls)))) : NaN
    d_corr_NGDP  = agg_dict["corr_NGDP"]
    d_corr_NAPL  = agg_dict["corr_NAPL"]
    # std(omG), moment 61 (added 2026-08-19): HP-filtered std of the nominal
    # goods expenditure share. Under Cobb-Douglas this IS the FGI (2023) omega_t,
    # so it identifies sigma_omg directly. NO fallback — a stale
    # aggregate_moments.csv would silently mis-score the new reallocation shock.
    haskey(agg_dict, "std_omG") || error("""
        aggregate_moments.csv is STALE: std_omG not found.
        It predates the 2026-08-19 goods-share moment. Regenerate with:
            python3 Data/build_reallocation_calibration.py
        (compute_data_moments.jl does not yet compute this moment, and it
        hardcodes omG = 0.57 — re-running it will clobber both entries.)
        """)
    d_std_omG    = agg_dict["std_omG"]
    # Moments 62-63 (2026-08-19): the relative-price channel that identifies cl.
    (haskey(agg_dict, "std_pigap") && haskey(agg_dict, "corr_pigap_om")) || error("""
        aggregate_moments.csv is STALE: std_pigap / corr_pigap_om not found.
        These are the goods-services relative price moments that identify the
        labour adjustment cost cl; without them cl is not interior-identified
        and a free search returns the upper bound. Regenerate with:
            python3 Data/build_reallocation_calibration.py
        """)
    d_std_pigap     = agg_dict["std_pigap"]
    d_corr_pigap_om = agg_dict["corr_pigap_om"]
    @printf "  Loaded sectoral_moments.csv + aggregate_moments.csv\n\n"
elseif get(ENV, "SMM_SMOKE", "0") == "1"
    # SMOKE-TEST MODE ONLY (set by smoke_test.jl): allow include-time syntax/
    # contract checks on a fresh bundle where the moment CSVs are built later
    # in the pipeline. These placeholders can NEVER produce estimation results:
    # run_smm_estimation.jl does not set SMM_SMOKE, so a real run against
    # missing/stale CSVs still hard-errors above.
    @printf "  [SMOKE MODE] moment CSVs absent — placeholder values for syntax check ONLY.\n"
    y_d = fill(0.04, NSEC); p_d = fill(0.02, NSEC); l_d = fill(0.03, NSEC)
    corr_YPH_d = fill(0.0, NSEC)
    d_std_GDP=0.0215; d_std_pi=0.0041; d_corr_GDPpi=-0.15
    d_omG=0.5304; d_std_Q=0.0520; d_autocorr_Q=0.75
    d_corr_GDPQ=-0.15; d_TBGDP=-0.02; d_std_TBGDP=0.025
    d_corr_NGDP=0.698; d_corr_NAPL=0.065; d_std_omG=0.0172
    d_std_pigap=0.0154; d_corr_pigap_om=0.376
else
    # NO placeholder fallback (removed 2026-07-08): estimating against
    # invented moments silently produces meaningless results.
    error("""
        Data moment files not found in $(DATA_DIR):
            sectoral_moments.csv / aggregate_moments.csv
        Generate them first:
            julia --project=. compute_data_moments.jl
        """)
end

# omG (goods expenditure share) is pinned by SS calibration → contributes 0 to the
# loss and carries no information for estimation.  Replaced by std(TB/GDP), which
# directly identifies the export and import elasticities.
# corr(Y_i,PH_i): negative under TFP shocks, positive under demand shocks — key identifier
# for supply vs demand decomposition per sector (12 new moments, Option-A).
# MOMENT_NAMES / N_MOMENTS / MOMENT_BLOCKS now live in utils.jl (shared with
# main_SOE_gap.jl, 2026-07-10) — utils.jl must be included before this file.

data_moments = [y_d; p_d; l_d; d_std_GDP; d_std_pi; d_corr_GDPpi;
                d_std_TBGDP; d_std_Q; d_autocorr_Q; d_corr_GDPQ; 1.0; 1.0; 1.0; corr_YPH_d;
                d_corr_NGDP; d_corr_NAPL; d_std_omG; d_std_pigap; d_corr_pigap_om;
                # 64-77: persistence (2026-08-21)
                d_ac_GDP; d_ac_pi; ac_y_d;
                # 78-84: literature-standard block
                d_rbar_YY; d_std_Lg; d_std_Ls;
                d_ratio_stdC; d_ratio_stdI; d_corr_CGDP; d_corr_IGDP]
@assert length(data_moments) == length(MOMENT_NAMES) "data_moments ($(length(data_moments))) ≠ MOMENT_NAMES ($(length(MOMENT_NAMES)))"

# GUARD: NaN anywhere in data_moments poisons the objective for ALL evaluations
# (best_obj initialises to NaN and is never updated). Catch this immediately.
let nan_dm = findall(isnan, data_moments)
    if !isempty(nan_dm)
        nan_names = [MOMENT_NAMES[k] for k in nan_dm]
        error("""
        NaN detected in data_moments at positions $(nan_dm):
          $(join(nan_names, ", "))

        The most common cause is a missing row in aggregate_moments.csv.
        Run  julia --project=. compute_data_moments.jl  to regenerate it.
        """)
    end
end


# =========================================================================== #
#  2.  PARAMETER BOUNDS (tighter, economically motivated)                     #
# =========================================================================== #

# PARAM_LABELS / N_THETA / CSV_PARAM_NAMES / LB / UB, the weighting matrix
# (build_weighting_matrix), and the report printers (print_param_table,
# print_fit_table, MOMENT_BLOCKS) now live in utils.jl — shared with
# main_SOE_gap.jl so the printed fit is the SAME objective the estimator
# minimizes (moved 2026-07-10).


# =========================================================================== #
#  3b. SETUP CONTRACT CHECKS                                                   #
# =========================================================================== #
# The dimensions 36 (params) and 60 (moments) are hard-coded in many places
# (LB/UB, layout comments, Σe indexing, slicing in smm_run, the decomposition
# printouts). Historically, adding one parameter or moment broke the run deep
# inside CMA-ES with an opaque BoundsError. This validates every cross-cutting
# invariant ONCE, up front, and fails with a message that says exactly what to
# update. Call it at the very top of smm_run.

# N_MOMENTS (= 60) is defined in utils.jl.

# Weighted contribution of one moment block, ψ' W ψ restricted to `rng`.
# Used by every decomposition printer so they all follow MOMENT_BLOCKS and can
# never again drift out of sync with N_MOMENTS (2026-08-19).
_blk0(ψ, W, rng) = dot(view(ψ, rng), view(W, rng, rng) * view(ψ, rng))

function validate_smm_setup(data_moments)
    errs = String[]
    # The decomposition printers sum MOMENT_BLOCKS; if that does not tile
    # 1:N_MOMENTS the printed blocks silently stop adding up to the objective.
    _covered = sort(vcat([collect(rng) for (rng, _) in MOMENT_BLOCKS]...))
    _covered == collect(1:N_MOMENTS) ||
        push!(errs, "MOMENT_BLOCKS does not tile 1:$(N_MOMENTS) exactly " *
                    "(covers $(length(_covered)) indices) — update MOMENT_BLOCKS in utils.jl")

    length(LB) == N_THETA ||
        push!(errs, "length(LB)=$(length(LB)) ≠ N_THETA=$N_THETA")
    length(UB) == N_THETA ||
        push!(errs, "length(UB)=$(length(UB)) ≠ N_THETA=$N_THETA")
    length(PARAM_LABELS) == N_THETA ||
        push!(errs, "length(PARAM_LABELS)=$(length(PARAM_LABELS)) ≠ N_THETA=$N_THETA")
    length(CSV_PARAM_NAMES) == N_THETA ||
        push!(errs, "length(CSV_PARAM_NAMES)=$(length(CSV_PARAM_NAMES)) ≠ N_THETA=$N_THETA — update the CSV name list for the new θ layout")
    N_THETA >= 36 ||
        push!(errs, "N_THETA=$N_THETA but the moment fn reads θ[31:36]; need ≥36")

    all(LB .< UB) ||
        push!(errs, "LB ≥ UB at indices $(findall(LB .>= UB)) " *
                    "(params: $(PARAM_LABELS[findall(LB .>= UB)]))")

    length(data_moments) == N_MOMENTS ||
        push!(errs, "length(data_moments)=$(length(data_moments)) ≠ N_MOMENTS=$N_MOMENTS")
    length(MOMENT_NAMES) == N_MOMENTS ||
        push!(errs, "length(MOMENT_NAMES)=$(length(MOMENT_NAMES)) ≠ N_MOMENTS=$N_MOMENTS")

    W = build_weighting_matrix(data_moments)
    size(W) == (N_MOMENTS, N_MOMENTS) ||
        push!(errs, "build_weighting_matrix returned $(size(W)), expected ($N_MOMENTS,$N_MOMENTS)")

    # Klein structural-param indices must be valid θ positions.
    bad_klein = filter(i -> i < 1 || i > N_THETA, _KLEIN_STRUCT_IDX)
    isempty(bad_klein) ||
        push!(errs, "_KLEIN_STRUCT_IDX has out-of-range entries $bad_klein for N_THETA=$N_THETA")

    # --- capital block (LPR semi-fixed capital, 2026-08-20) ----------------- #
    # The baseline NamedTuple reads alphaK_i / Kbar_i / chiI_i / nuK out of
    # params_jl.mod. A params file written BEFORE the capital block simply has no
    # such lines, and pvec() returns NaN for a missing name. Those NaNs would
    # propagate silently into the steady-state solve and produce a converged-
    # looking but meaningless model, so check here and fail loudly instead.
    # Fix: re-run  julia --project=. main_SOE_gap.jl  to regenerate params_jl.mod.
    let pf = joinpath(@__DIR__, "mod", "params_jl.mod")
        if isfile(pf)
            txt = read(pf, String)
            for nm in ["alphaK_1", "Kbar_1", "chiI_1", "RKss_1", "nuK", "PIinv_ss"]
                occursin(Regex("(^|\\n)\\s*$(nm)\\s*="), txt) ||
                    push!(errs, "params_jl.mod has no `$nm` — it predates the capital " *
                                "block. Re-run main_SOE_gap.jl before estimating.")
            end
        else
            push!(errs, "mod/params_jl.mod not found — run main_SOE_gap.jl first.")
        end
    end

    # The calibration CSVs the capital block reads must exist and be coherent.
    let cf = joinpath(@__DIR__, "..", "..", "Data", "sector_calibration.csv"),
        kf = joinpath(@__DIR__, "..", "..", "Data", "capital_calibration.csv")
        if isfile(cf)
            d = CSV.read(cf, DataFrame)
            if !hasproperty(d, :alpha_K) || !hasproperty(d, :chi_I)
                push!(errs, "sector_calibration.csv lacks alpha_K / chi_I — " *
                            "run python3 Data/build_sector_calibration.py")
            else
                abs(sum(d.chi_I) - 1) < 1e-8 ||
                    push!(errs, "chi_I sums to $(sum(d.chi_I)), expected 1")
                lab = 1 .- d.alpha .- d.alpha_V .- d.alpha_K
                all(>(0), lab) ||
                    push!(errs, "non-positive labour weight in sectors $(findall(<=(0), lab))")
            end
        else
            push!(errs, "Data/sector_calibration.csv not found")
        end
        isfile(kf) || push!(errs, "Data/capital_calibration.csv not found — " *
                                  "run python3 Data/build_sector_calibration.py")
    end

    if !isempty(errs)
        error("""
        SMM setup is inconsistent — fix before estimating:

          - $(join(errs, "\n          - "))

        These dimensions must stay in sync whenever you add a parameter or a
        moment. The usual edit points are:
          params  → LB, UB, PARAM_LABELS, N_THETA, the θ[...] slices in
                    smm_model_moments, default_theta0, _KLEIN_STRUCT_IDX
          moments → MOMENT_NAMES, N_MOMENTS, build_weighting_matrix, the return
                    vector of smm_model_moments, the ψ-decomposition printouts
        """)
    end
    @printf "  Setup contract checks passed: %d params, %d moments.\n" N_THETA N_MOMENTS
    return nothing
end


# Atomic CSV write: write to a temp file in the same directory, then rename.
# A `mv` is atomic on the same filesystem, so a job killed mid-write never
# leaves a half-written checkpoint that would poison the next warm start.
function atomic_write_csv(path::AbstractString, df)
    mkpath(dirname(path))
    tmp = path * ".tmp_$(getpid())_$(Threads.threadid())"
    CSV.write(tmp, df)
    mv(tmp, path; force=true)
    return path
end

# Persist the current best θ to the checkpoint file (warm-start source). Called
# live, on every improvement, so a 1–2 day cluster run that dies loses at most
# the current generation rather than everything since the last manual save.
function save_checkpoint(θ::AbstractVector, obj::Real)
    # Use CSV_PARAM_NAMES (single source of truth for CSV output) — a
    # hardcoded name list here silently went stale when kappaw (θ[36]) was
    # added on 2026-07-08 and killed the run at the first checkpoint write.
    length(θ) == length(CSV_PARAM_NAMES) ||
        error("save_checkpoint: length(θ)=$(length(θ)) ≠ length(CSV_PARAM_NAMES)=$(length(CSV_PARAM_NAMES))")
    df = DataFrame(
        param = CSV_PARAM_NAMES,
        value = collect(Float64, θ),
        obj_hat = vcat([Float64(obj)], fill(NaN, length(θ)-1)))
    try
        atomic_write_csv(joinpath(ESTIMATION_DIR, "smm_checkpoint.csv"), df)
        # Plain-text twins (house convention, cf. third-year paper best_sol.txt):
        # best_sol.txt  — current best θ, one "name  value" line per parameter
        # min_loss.txt  — loss at the best θ + timestamp
        _ts = Libc.strftime("%Y-%m-%d %H:%M:%S", time())
        open(joinpath(ESTIMATION_DIR, "best_sol.txt"), "w") do io
            @printf(io, "# NK-IOSOE SMM best solution  |  obj = %.8f  |  %s\n", Float64(obj), _ts)
            for (nm, v) in zip(CSV_PARAM_NAMES, θ)
                @printf(io, "%-16s  %.10g\n", nm, Float64(v))
            end
        end
        open(joinpath(ESTIMATION_DIR, "min_loss.txt"), "w") do io
            @printf(io, "%.8f\n# loss at best_sol.txt  |  %s\n", Float64(obj), _ts)
        end
        # Stamp the objective definition this obj was computed under, so a
        # later main_SOE_gap.jl run can tell whether the stored number is still
        # comparable (2026-08-19; see OBJECTIVE PROVENANCE in utils.jl).
        write_objective_provenance(ESTIMATION_DIR,
            objective_dep_files(SCRIPT_DIR, joinpath(SCRIPT_DIR, "mod"), DATA_DIR))
    catch err
        @printf "  [warn] checkpoint write failed: %s\n" sprint(showerror, err)
    end
    return nothing
end


# =========================================================================== #
#  4.  BASELINE CALIBRATION STRUCT                                            #
# =========================================================================== #

function build_baseline(context::Dynare.Context,
                         endo_names::Vector{String},
                         d_std_Y, d_std_PH, d_std_L, d_TBGDP, d_omG)
    function pvec(nm)
        idx = param_idx(context, nm)
        idx === nothing && return NaN
        tid = _tid()
        p = _SMM_PARAMS_READY[tid][] ? _SMM_PARAMS[tid][] : _load_smm_params!(context)
        (isempty(p) || idx > length(p)) ? NaN : p[idx]
    end
    nsec = NSEC

    # ---- SMM_KAPPA_SCALE: diagnostic scale on the sectoral Rotemberg costs ---
    #
    # WHY. std(PH_i) is 9.1% of the objective and corr(Y_i,PH_i) is 38.8% — 48%
    # between them — and they may share one mechanism: the model's sectoral
    # PRICES barely move, so the output-price correlation collapses toward zero.
    # Model vs data std(PH_i) at the 2026-08-25 estimate: Agriculture 3.6x too
    # smooth, Manufacturing 2.9x, Finance 2.4x, Transport 2.3x, Utilities 2.1x.
    # The tell is MINING — the one sector with NO Phillips curve (law of one
    # price, no kappa applied) — and it is the one sector that matches (0.132
    # data vs 0.145 model).
    #
    # The kappa_i are CALIBRATED from price-change frequencies (8.0 Agriculture
    # to 508.6 Public Admin), never estimated, so nothing has ever tested whether
    # their LEVEL is right. This scale multiplies all twelve at once.
    #
    # DIAGNOSTIC ONLY, and read it sceptically: the errors go BOTH ways
    # (Construction and Business Services are ~2x too VOLATILE), so a single
    # global scale is a blunt instrument and may not help. The question it
    # answers is narrow: do std(PH_i) and corr(Y_i,PH_i) improve TOGETHER when
    # prices are made more flexible? If yes, they share a mechanism and a free
    # scale is worth building. If only one moves, they do not, and this is dead.
    #
    # Applied HERE, before the pvec reads below, so the context parameters and
    # baseline.modkappa (which feeds the steady state) can never disagree.
    let ks = parse(Float64, get(ENV, "SMM_KAPPA_SCALE", "1.0"))
        if ks != 1.0
            ks > 0 || error("SMM_KAPPA_SCALE must be > 0, got $(ks)")
            # THREADS=1 IS REQUIRED HERE, not merely preferred. set_param! writes
            # into _SMM_PARAMS[tid], a PER-THREAD vector, and build_baseline runs
            # on thread 1 only. Under multithreading the scale would apply to
            # thread 1's copy and every other thread would silently solve at the
            # ORIGINAL kappas — a partial, non-reproducible calibration that no
            # output would reveal. (Separately, threading already causes ~36% of
            # evaluations to return a stale decision rule, so --threads=1 is the
            # standing recommendation anyway.)
            Threads.nthreads() == 1 || error("""
                SMM_KAPPA_SCALE=$(ks) requires --threads=1 (running with $(Threads.nthreads())).
                set_param! writes to the per-thread _SMM_PARAMS vector and only
                thread 1 runs build_baseline, so other threads would keep the
                ORIGINAL kappas and the run would be silently inconsistent.
                Re-run with: julia --threads=1 --project=. ...
                """)
            base = [pvec("kappa_$(i)") for i in 1:nsec]
            for i in 1:nsec
                isnan(base[i]) && continue
                set_param!(context, "kappa_$(i)", ks * base[i])
            end
            # DO NOT invalidate and reload the cache here. `set_param!` writes
            # into _SMM_PARAMS itself — see the note at smm_model_moments.jl:636
            # ("context.work.params was unreliable on ARM Macs, so this file
            # keeps its own vector"), and resolve_first_order! pushes _SMM_PARAMS
            # into the context before any solver runs. `_load_smm_params!` refills
            # _SMM_PARAMS *from context.work.params*, i.e. from the ORIGINAL
            # params_jl.mod values, so calling it after set_param! silently
            # DISCARDS the writes.
            #
            # That is exactly what happened on the first attempt: the 2026-08-25
            # kappa sweep returned byte-identical moments at 0.1x and 4x because
            # a reload added "for safety" wiped every set_param! a few lines
            # earlier. The read-back below is what caught it; keep it.
            # READ BACK, do not trust the write. Printing ks*base would only
            # prove I can multiply; it says nothing about whether set_param!
            # reached the array the solver reads. A byte-identical sweep is the
            # classic symptom of a parameter that never arrived.
            rb = [pvec("kappa_$(i)") for i in 1:nsec]
            @printf "  [kappa scale] SMM_KAPPA_SCALE=%.3f applied to all %d sectoral Rotemberg costs\n" ks nsec
            @printf "                kappa_1  %8.3f -> requested %8.3f, READ BACK %8.3f\n" base[1] ks*base[1] rb[1]
            @printf "                kappa_12 %8.3f -> requested %8.3f, READ BACK %8.3f\n" base[nsec] ks*base[nsec] rb[nsec]
            bad = count(i -> !isnan(rb[i]) && abs(rb[i] - ks*base[i]) > 1e-9*max(1.0, abs(ks*base[i])), 1:nsec)
            if bad > 0
                error("""
                    SMM_KAPPA_SCALE did not take: $(bad) of $(nsec) kappa_i read back
                    at their OLD values after set_param!. The scale would be silently
                    ignored and the sweep would return identical moments at every
                    grid point. Do not interpret such a sweep as "kappa does not
                    matter" — it means the parameter never reached the solver.
                    """)
            end
        end
    end

    ss_vec = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let i = findfirst(==(nm), endo_names); i === nothing ? 1.0 : ss_vec[i]; end
    (
        nsec=nsec, goods=GOODS, services=SERVICES,
        Y_ss=[get_ss("Y_$(i)") for i in 1:nsec],
        GDP_ss=get_ss("GDP"),
        ombar_val=let v=pvec("ombar"); isnan(v) ? d_omG : v; end,
        tb_target=d_TBGDP,
        modalpha=[pvec("alpha_$(i)") for i in 1:nsec],
        modalphaV=[pvec("alphaV_$(i)") for i in 1:nsec],
        # LPR semi-fixed capital, carried through from params_jl.mod. No NaN
        # fallbacks: a params file without these was written by a pre-2026-08-20
        # build and must not be silently reinterpreted as a no-capital model.
        modalphaK=[pvec("alphaK_$(i)") for i in 1:nsec],
        modKbar=[pvec("Kbar_$(i)") for i in 1:nsec],
        modchiI=[pvec("chiI_$(i)") for i in 1:nsec],
        nuK_val=pvec("nuK"),
        # Exogenous government demand (2026-08-20). modGi is the solved SS level
        # vector written by main_SOE_gap.jl; modgammaG and gshare_target let the
        # estimation re-solve the steady state with the same government block.
        modGi=[pvec("Gi_$(i)") for i in 1:nsec],
        modgammaG=let g=[pvec("Gi_$(i)") for i in 1:nsec]
            s = sum(g); s > 0 ? g ./ s : zeros(nsec)
        end,
        gshare_target=pvec("gshare_target"),
        subsMC_val=let v=pvec("subsMC"); (isnan(v) || v <= 0) ? 1.0 : v end,
        modbeta=[pvec("beta_$(i)_$(j)") for i in 1:nsec, j in 1:nsec],
        modgammag=[pvec("gammag_$(i)") for i in 1:nsec],
        modgammas=[pvec("gammas_$(i)") for i in 1:nsec],
        modvarrho=[pvec("varrho_$(i)") for i in 1:nsec],
        modchiX=[pvec("chiX_$(i)") for i in 1:nsec],
        modkappa=[pvec("kappa_$(i)") for i in 1:nsec],
        gamma_val=pvec("gamma"), psi_val=pvec("psi"), chi_val=1.0,
        epsilon_val=pvec("epsilon"), beta_val=pvec("beta"),
        PVstar_ss=pvec("PVstar_ss"), sigmaH_val=pvec("sigmaH"),
        etastar_val=pvec("etastar"), omegaX_val=pvec("omegaX"),
        ystar_ss_val=pvec("ystar_ss"),
        data_std_Y=d_std_Y, data_std_PH=d_std_PH, data_std_L=d_std_L,
        # Diagonal Σe positions to activate, resolved by shock NAME (fixes C1).
        active_exo_idx=active_shock_indices(context, nsec),
    )
end


# =========================================================================== #
#  5.  INITIAL θ FROM CONTEXT                                                 #
# =========================================================================== #

function default_theta0(context::Dynare.Context)
    pv(nm) = let idx = param_idx(context, nm)
        if idx === nothing; 0.0
        else
            tid = _tid()
            p = _SMM_PARAMS_READY[tid][] ? _SMM_PARAMS[tid][] : _load_smm_params!(context)
            (isempty(p) || idx > length(p)) ? 0.0 : p[idx]
        end
    end
    # Option-A layout: 38 params.
    #   θ[36] = kappaw                            (2026-07-08)
    #   θ[37:38] = lambda_A, lambda_om            (2026-08-20, shock-mix scales)
    # The λ entries have no counterpart in params_jl.mod — they multiply the
    # measured shock vectors inside the objective rather than being written to
    # the model — so they are seeded at 1.0 here, which reproduces the measured
    # sizes exactly.
    θ = [pv("ilabcosts"); pv("epsY_1"); pv("epsM_1"); log(pv("kappaV"));
         pv("rho_om1"); pv("rho_tfp1");
         [pv("isigma_tfp_$(i)") for i in 1:12];
         [max(pv("sigma_om_$(i)"), 0.01) for i in 1:12];
         pv("rho_pvstar"); pv("sigma_pvstar"); pv("rho_zeta"); pv("sigma_zeta");
         pv("etastar");
         let k = pv("kappaw"); k > 0 ? k : 115.0 end;
         1.0; 1.0]
    # Guard: this constructor and the θ layout in utils.jl must stay in sync.
    # They are edited in different files, so a mismatch shows up as a
    # BoundsError several hundred lines away from the cause.
    length(θ) == N_THETA || error(
        "default_theta0 built a $(length(θ))-vector but N_THETA=$(N_THETA). " *
        "Update BOTH default_theta0 (smm_estimation.jl) and PARAM_LABELS / " *
        "CSV_PARAM_NAMES / LB / UB (utils.jl) when changing the θ layout.")
    return θ
end




# =========================================================================== #
#  7.  PER-THREAD SS CACHE                                                    #
# =========================================================================== #

mutable struct SSCache
    epsY::Float64; epsM::Float64; etastar::Float64
    ss_vec::Vector{Float64}; valid::Bool
end
SSCache() = SSCache(NaN, NaN, NaN, Float64[], false)

const _SS_CACHES = [SSCache() for _ in 1:max(1, Threads.nthreads())]
const _SS_TOL    = 1e-4

function recompute_ss_cached!(context, epsY, epsM, baseline, endo_names; etastar=nothing)
    tid   = min(Threads.threadid(), length(_SS_CACHES))
    cache = _SS_CACHES[tid]
    eta_eff = something(etastar, baseline.etastar_val)
    eta_cached = isnan(cache.etastar) ? baseline.etastar_val : cache.etastar
    if cache.valid && abs(epsY-cache.epsY) < _SS_TOL && abs(epsM-cache.epsM) < _SS_TOL &&
       abs(eta_eff - eta_cached) < _SS_TOL
        for i in 1:baseline.nsec
            set_param!(context, "epsY_$(i)", epsY)
            set_param!(context, "epsM_$(i)", epsM)
        end
        ss_mut = context.results.model_results[1].trends.endogenous_steady_state
        length(cache.ss_vec) == length(ss_mut) && copyto!(ss_mut, cache.ss_vec)
        return true
    end
    ok = recompute_ss!(context, epsY, epsM, baseline, endo_names; etastar=etastar)
    if ok
        cache.epsY    = epsY; cache.epsM = epsM; cache.etastar = eta_eff
        cache.ss_vec  = copy(context.results.model_results[1].trends.endogenous_steady_state)
        cache.valid   = true
    end
    return ok
end


# =========================================================================== #
#  8.  IN-PLACE LYAPUNOV SOLVER                                               #
# =========================================================================== #

function dlyap_inplace!(P::Matrix{Float64}, A::AbstractMatrix{Float64},
                         Q::Matrix{Float64};
                         tmp1::Matrix{Float64}, tmp2::Matrix{Float64},
                         Ak::Matrix{Float64}, maxiter::Int=100, tol::Float64=1e-12)
    copyto!(P, Q); copyto!(Ak, A)
    n = size(P, 1)
    for _ in 1:maxiter
        mul!(tmp1, Ak, P); mul!(tmp2, tmp1, Ak')
        err = 0.0
        @inbounds for j in 1:n, i in 1:n
            v = tmp2[i,j] + P[i,j]; err = max(err, abs(v - P[i,j])); P[i,j] = v
        end
        @inbounds for j in 1:n, i in 1:j-1
            avg = 0.5*(P[i,j]+P[j,i]); P[i,j]=avg; P[j,i]=avg
        end
        copyto!(tmp1, Ak); mul!(Ak, tmp1, tmp1)
        err / (norm(P, Inf) + 1e-14) < tol && break
    end
    @inbounds for j in 1:n, i in 1:j-1
        avg = 0.5*(P[i,j]+P[j,i]); P[i,j]=avg; P[j,i]=avg
    end
    return P
end


# =========================================================================== #
#  9.  PER-THREAD SCRATCH ARRAYS                                              #
# =========================================================================== #

mutable struct SMMScratch
    Σe::Matrix{Float64}; Q_lyap::Matrix{Float64}
    P::Matrix{Float64};  tmp1::Matrix{Float64}
    tmp2::Matrix{Float64}; Ak::Matrix{Float64}
    endo_idx::Dict{String,Int}
    hp_w_var::Vector{Float64}; hp_w_lag1::Vector{Float64}
    initialized::Bool
end
SMMScratch() = SMMScratch(zeros(0,0),zeros(0,0),zeros(0,0),zeros(0,0),zeros(0,0),zeros(0,0),
                           Dict{String,Int}(),Float64[],Float64[],false)

const _SCRATCH = [SMMScratch() for _ in 1:max(1, Threads.nthreads())]

function _get_scratch!(n_exo::Int, n_state::Int, endo_names::Vector{String})
    tid = min(Threads.threadid(), length(_SCRATCH))
    s = _SCRATCH[tid]
    if !s.initialized || size(s.Σe,1) != n_exo || size(s.P,1) != n_state
        s.Σe      = zeros(n_exo, n_exo)
        s.Q_lyap  = zeros(n_state, n_state)
        s.P       = zeros(n_state, n_state)
        s.tmp1    = zeros(n_state, n_state)
        s.tmp2    = zeros(n_state, n_state)
        s.Ak      = zeros(n_state, n_state)
        s.endo_idx = Dict(nm => i for (i,nm) in enumerate(endo_names))
        s.hp_w_var, s.hp_w_lag1 = build_hp_weights(1600.0, 256)
        s.initialized = true
    end
    return s
end


# =========================================================================== #
#  10. OPTIMIZED MOMENT FUNCTION                                              #
# =========================================================================== #

"""
Opt-in diagnostic sink (2026-08-24). When set to `Ref`-held `true`,
`smm_model_moments` stashes the HP-filtered covariance matrix it already
computed, together with the variable names it is indexed by, so a diagnostic
script can look at objects the 84-moment vector does not expose — cross-sector
correlation structure, principal components, loadings on aggregates.

`nothing` (the default) means the branch is never taken, so the estimation pays
literally nothing for this. NEVER leave it armed inside a CMA-ES loop: it copies
a ~50x50 matrix per evaluation and would be retained across millions of calls.
"""
const MOMENT_DIAG_SINK = Ref{Any}(nothing)

"""
The twelve MEASURED sectoral TFP persistences rho_hat_{A,i}, read from
`rho_A_init` in the sectoral shock calibration CSV at setup (2026-08-24).

The model's `rho_tfp1_i` are then set to `clamp(lambda_rho * rho_hat_i, 0, 0.98)`
with `lambda_rho = θ[6]`, exactly mirroring how `lambda_A` scales the twelve
measured shock SIZES: the cross-sectional PATTERN is measured from the data and
pinned, one free scalar absorbs the level. No new free parameters.

`nothing` means the CSV predates this change; the objective then falls back to a
single common persistence, `θ[6]` reverts to its old meaning, and setup prints a
loud warning. That fallback exists so an old calibration file still runs — it is
NOT the intended configuration, because a common rho cannot fit sectoral
autocorrelations that run from 0.11 to 0.89 in the data.
"""
const RHO_A_MEASURED = Ref{Union{Nothing,Vector{Float64}}}(nothing)

function smm_model_moments(θ::AbstractVector{<:Real}, context, baseline, endo_names)
    nsec = baseline.nsec; NAN58 = fill(NaN, N_MOMENTS)   # name kept; length = N_MOMENTS (60)
    # Option-A layout (36 params):
    #   θ[1:4]  = ilabcosts, epsY, epsM, log(kappaV)
    #   θ[5]    = rho_om (common persistence for all 12 sectoral demand shocks)
    #   θ[6]    = lambda_rho, scale on the 12 MEASURED sectoral persistences
    #             (was the single common rho_A until 2026-08-24)
    #   θ[7:18] = isigma_tfp_1:12
    #   θ[19:30]= sigma_om_1:12
    #   θ[31:35]= rho_pvstar, sigma_pvstar, rho_zeta, sigma_zeta, etastar
    #   θ[36]   = kappaw (Rotemberg wage stickiness; 0 = flexible)
    ilabcosts=θ[1]; epsY=θ[2]; epsM=θ[3]; kappaV=exp(θ[4])
    # θ[6]: lambda_rho, the SCALE on the measured sectoral persistences when
    # RHO_A_MEASURED is populated; the single common rho_A only in the legacy
    # fallback. See RHO_A_MEASURED above.
    rho_om=θ[5]; lambda_rho=θ[6]
    _rhoA   = RHO_A_MEASURED[]
    # POWER transform, NOT multiplicative (fixed 2026-08-24 after the first
    # sector-specific run scored 8.18 against the common-rho 6.74).
    #
    # rho_i = rho_hat_i ^ (1 / lambda_rho)
    #
    # The first design was clamp(lambda_rho * rho_hat_i, 0, 0.98). Multiplying a
    # parameter bounded on [0,1) and then clamping COLLAPSES the top of the
    # distribution: at the estimate lambda_rho went to its bound of 2.0 and SEVEN
    # of twelve sectors clamped to 0.98, flattening exactly the cross-sectional
    # dispersion the change exists to introduce. That run tested a mangled
    # specification, not this hypothesis.
    #
    # The power map fixes it: it is monotone, it PRESERVES THE ORDERING of the
    # measured persistences, it maps (0,1) into (0,1) so no clamping ever binds,
    # and lambda_rho = 1 is the identity, so the measured values are nested. Above
    # 1 it pushes every sector toward unity while keeping them ordered and
    # distinct; below 1, toward zero.
    rho_A_i = _rhoA === nothing ? nothing :
              [clamp(r <= 0 ? 0.0 : r^(1 / lambda_rho), 0.0, 0.98) for r in _rhoA]
    # Scalar used for the legacy path, the (unused) global rho_tfp1, and the
    # validity screen below. With measured values this is their mean.
    rho_A = _rhoA === nothing ? lambda_rho : mean(rho_A_i)
    isigma_tfp  = @view θ[7:18]
    sigma_om_vec= @view θ[19:30]
    rho_pvstar=θ[31]; sigma_pvstar=θ[32]; rho_zeta=θ[33]; sigma_zeta=θ[34]
    etastar = length(θ) >= 35 ? θ[35] : baseline.etastar_val
    # kappaw (θ[36], added 2026-07-08): Rotemberg wage stickiness. Fallback to
    # the params_jl.mod value only for legacy 35-length θ vectors.
    kappaw = length(θ) >= 36 ? θ[36] :
             let v = get_param_val(context,"kappaw"); isnan(v) ? 115.0 : v end

    # SHOCK-MIX SCALES (θ[37:38], added 2026-08-20). These are FREE and multiply
    # the PINNED, data-measured shock vectors. lambda = 1 reproduces the measured
    # calibration exactly, so the pre-2026-08-20 model is nested at (1,1).
    #
    # Why they exist: with all 24 sectoral shock sizes pinned, nothing in the
    # free parameter set could change the supply/demand mix, and the
    # corr(Y_i,PH_i) block — 40% of the objective — was structurally
    # unimprovable. Scaling preserves the measured cross-sectional pattern (what
    # the rank-correlation moments discipline) and frees only the level and mix.
    lambda_A  = length(θ) >= 37 ? θ[37] : 1.0
    lambda_om = length(θ) >= 38 ? θ[38] : 1.0

    # θ-REDUCTION: override the pinned params with their measured/fixed values.
    # NOTE these are COPIES, not views: the λ scaling must not write back into
    # SMM_PIN, or every subsequent evaluation would compound the last one's λ.
    if SMM_PIN[] !== nothing
        _p = SMM_PIN[]
        # Take a value from SMM_PIN only if it is actually pinned. epsY is free
        # as of 2026-08-21; reading it from _p regardless would have let CMA-ES
        # search a parameter the model never received.
        epsY = _FREE_EPSY ? θ[2] : _p[2]
        epsM = _FREE_EPSM ? θ[3] : _p[3]
        isigma_tfp   = _p[7:18]  .* lambda_A
        sigma_om_vec = _p[19:30] .* lambda_om
        rho_pvstar   = _p[31]; sigma_pvstar = _p[32]
    else
        isigma_tfp   = isigma_tfp   .* lambda_A
        sigma_om_vec = sigma_om_vec .* lambda_om
    end
    (lambda_A <= 0 || lambda_om <= 0) && return NAN58, false

    (!(0<epsY<5)||!(0<epsM<2)||ilabcosts<=0||kappaV<=0||rho_om<0||rho_om>=1||
     any(<(0),sigma_om_vec)||rho_A<0||rho_A>=1||any(<(0),isigma_tfp)||rho_pvstar<0||rho_pvstar>=1||
     sigma_pvstar<0||rho_zeta<0||rho_zeta>=1||sigma_zeta<0||
     !(0.1<etastar<8.0)||kappaw<0||kappaw>1e4) && return NAN58, false

    set_param!(context,"kappaw",kappaw)
    set_param!(context,"ilabcosts",ilabcosts); set_param!(context,"kappaV",kappaV)
    set_param!(context,"rho_om1",rho_om)
    set_param!(context,"etastar",etastar)
    set_param!(context,"rho_tfp1",rho_A);      set_param!(context,"rho_pvstar",rho_pvstar)
    set_param!(context,"sigma_pvstar",sigma_pvstar); set_param!(context,"rho_zeta",rho_zeta)
    set_param!(context,"sigma_zeta",sigma_zeta)
    for i in 1:nsec
        set_param!(context,"isigma_tfp_$(i)",isigma_tfp[i])
        set_param!(context,"sigma_om_$(i)",sigma_om_vec[i])   # now EXISTS in the model (fixes C2)
        # SECTOR-SPECIFIC TFP PERSISTENCE (2026-08-24). The .mod always had
        # rho_tfp1_i per sector; this line used to overwrite all twelve with one
        # common value, which is what made the autocorr(Y_i) block (19% of the
        # objective) unfittable — the model was 4x too persistent for Agriculture
        # and 3x too transitory for Finance at the same time.
        set_param!(context,"rho_tfp1_$(i)", rho_A_i === nothing ? rho_A : rho_A_i[i])
        # θ[1] DRIVES THE LABOUR ADJUSTMENT COST (2026-08-19).
        #
        # It used to set only `ilabcosts`, which is a DEAD parameter: in
        # NK_SOE_lev_gap2.mod `ilabcosts` appears in exactly one place, the
        # definition of the reporting variable `Lab_costs` (line 319), and
        # `Lab_costs` is used in no other equation. The live parameter is
        # `cl_i`, which enters BOTH the labour-agency FOC (line 557) and
        # labour-market clearing (line 302) — and it was pinned to 0 by
        # main_SOE_gap.jl's `SMM_CL` default, so the friction was switched off
        # and θ[1] had no effect on anything at any value.
        #
        # `cl_i` IS Ferrante, Graves & Iacoviello's (2023 JME) hiring cost c,
        # one-for-one: their FOC + envelope condition is reproduced term for
        # term at .mod:557-560, and their aggregate labour-market clearing
        # (their eq. 24) at .mod:302-307. FGI estimate c = 19.1 (s.e. 12.6).
        # θ[1]'s bounds [0.001, 50] bracket that comfortably.
        #
        # The steady state does not depend on cl (at L_i/L_i(-1) = 1 every
        # adjustment term vanishes), so this cannot disturb the calibration.
        set_param!(context,"cl_$(i)",ilabcosts)
    end

    epsY_prev    = get_param_val(context,"epsY_1"); epsM_prev = get_param_val(context,"epsM_1")
    etastar_prev = get_param_val(context,"etastar")
    need_ss = abs(epsY-epsY_prev)>1e-8 || abs(epsM-epsM_prev)>1e-8 ||
              abs(etastar - (isnan(etastar_prev) ? baseline.etastar_val : etastar_prev)) > 1e-8
    for i in 1:nsec; set_param!(context,"epsY_$(i)",epsY); set_param!(context,"epsM_$(i)",epsM); end
    if need_ss
        ok = recompute_ss_cached!(context, epsY, epsM, baseline, endo_names; etastar=etastar)
        !ok && return NAN58, false
    end

    success, T, R = _resolve_cached!(context, θ)
    !success && return NAN58, false

    n_exo   = size(R, 2)   # should be 28 after Option-A mod recompile
    sr      = context.models[1].i_bkwrd_b
    n_state = length(sr)
    sc      = _get_scratch!(n_exo, n_state, endo_names)

    # Σe: activate shocks by NAME (fixes C1). The diagonal positions were
    # resolved once in build_baseline from the model's own exogenous ordering
    # (eps_i, eps_pvstar, eps_zeta, epsA_1:12, eps_om_1:12; oil & labour-supply off),
    # so a change to the .mod shock list can never silently scramble Σe again.
    fill!(sc.Σe, 0.0)
    @inbounds for k in baseline.active_exo_idx
        k <= n_exo && (sc.Σe[k,k] = 1.0)
    end

    Tsr = T[sr,:]; Rsr = R[sr,:]
    RΣ  = Rsr * sc.Σe   # 78×17
    mul!(sc.Q_lyap, RΣ, Rsr')
    @inbounds for j in 1:n_state, i in 1:j-1
        avg=0.5*(sc.Q_lyap[i,j]+sc.Q_lyap[j,i]); sc.Q_lyap[i,j]=avg; sc.Q_lyap[j,i]=avg
    end

    dlyap_inplace!(sc.P, Tsr, sc.Q_lyap; tmp1=sc.tmp1, tmp2=sc.tmp2, Ak=sc.Ak)
    # Reject: negative variances, NaN, OR non-finite / explosive entries. The
    # doubling iteration returns a finite-but-huge P when the state transition is
    # near-unit-root (BK borderline); 1e12 on a normalised state is unphysical.
    (any(diag(sc.P).<-1e-10) || !all(isfinite, sc.P) ||
     maximum(abs, diag(sc.P)) > 1e12) && return NAN58, false

    # HP-filtered covariance on the ~40 needed variables only (25× faster)
    needed_names = vcat(["Y_$(i)" for i in 1:nsec], ["PH_$(i)" for i in 1:nsec],
                        ["L_$(i)" for i in 1:nsec],
                        ["GDP","GDP_vol","pi","Q","TB","N","om_g","pi_g","pi_s",
                         # 2026-08-21, for the great-ratio block 81-84. EInv is
                         # NOMINAL investment expenditure, so real investment is
                         # EInv/PI_inv and both are needed.
                         "C","EInv","PI_inv"])
    ei          = sc.endo_idx
    needed_idx  = [get(ei, nm, 0) for nm in needed_names]
    valid_mask  = needed_idx .> 0
    nidx_valid  = needed_idx[valid_mask]

    T_sub = Matrix{Float64}(T[nidx_valid, :])
    R_sub = Matrix{Float64}(R[nidx_valid, :])
    RsubΣRsub = R_sub * sc.Σe * R_sub'
    RsubΣRsub = (RsubΣRsub + RsubΣRsub') / 2

    Γ_sub, Γ1_sub = hp_filtered_cov_fast(
        Matrix{Float64}(Tsr), sc.Q_lyap, T_sub, RsubΣRsub,
        sc.hp_w_var, sc.hp_w_lag1)

    # Map back to full needed_names ordering
    n_needed = length(needed_names)
    Γ_v  = zeros(n_needed, n_needed)
    Γ1_v = zeros(n_needed, n_needed)
    sub_pos = findall(valid_mask)
    for (si,pi) in enumerate(sub_pos), (sj,pj) in enumerate(sub_pos)
        Γ_v[pi,pj]=Γ_sub[si,sj]; Γ1_v[pi,pj]=Γ1_sub[si,sj]
    end

    # Local index lookup for needed_names
    ei_sub = Dict(nm=>i for (i,nm) in enumerate(needed_names))
    ys     = context.results.model_results[1].trends.endogenous_steady_state

    # Diagnostic sink — see MOMENT_DIAG_SINK above. Disarmed by default.
    if MOMENT_DIAG_SINK[] !== nothing
        MOMENT_DIAG_SINK[] = (Gamma = copy(Γ_v), names = copy(needed_names))
    end

    @inline pstd(vn) = let k=get(ei_sub,vn,0)
        k==0 ? 0.0 : sqrt(max(Γ_v[k,k],0.0)) / max(abs(ys[needed_idx[k]]),1e-12)
    end
    @inline xcorr(v1,v2) = let i1=get(ei_sub,v1,0), i2=get(ei_sub,v2,0)
        (i1==0||i2==0) ? NaN :
        let d=sqrt(max(Γ_v[i1,i1],0.0)*max(Γ_v[i2,i2],0.0))
            d<1e-15 ? 0.0 : clamp(Γ_v[i1,i2]/d,-1.0,1.0)
        end
    end

    # First-order autocorrelation of any needed variable. Gamma1_v is the
    # LAG-1 HP-filtered covariance, already computed alongside Gamma_v, so the
    # persistence moments cost nothing extra.
    @inline pac(vn) = let k=get(ei_sub,vn,0)
        (k==0 || Γ_v[k,k] <= 1e-15) ? NaN : clamp(Γ1_v[k,k]/Γ_v[k,k], -1.0, 1.0)
    end

    i_Q  = get(ei_sub,"Q",0)
    i_TB = get(ei_sub,"TB",0)
    acQ  = (i_Q>0 && Γ_v[i_Q,i_Q]>1e-15) ? Γ1_v[i_Q,i_Q]/Γ_v[i_Q,i_Q] : NaN

    # std(TB/GDP): normalize by steady-state GDP (TB is a level variable in the model)
    GDP_ss = baseline.GDP_ss > 0 ? baseline.GDP_ss : 1.0
    std_TBGDP = (i_TB>0) ? sqrt(max(Γ_v[i_TB,i_TB],0.0)) / GDP_ss : 0.0

    std_Y  = [pstd("Y_$(i)")  for i in 1:nsec]
    std_PH = [pstd("PH_$(i)") for i in 1:nsec]
    std_L  = [pstd("L_$(i)")  for i in 1:nsec]

    # corr(Y_i, PH_i): negative under TFP shocks, positive under demand shocks
    corr_YPH = [xcorr("Y_$(i)", "PH_$(i)") for i in 1:nsec]

    dY=baseline.data_std_Y; dPH=baseline.data_std_PH; dL=baseline.data_std_L
    vy=isfinite.(std_Y).&isfinite.(dY)
    vp=isfinite.(std_PH).&isfinite.(dPH)
    vl=isfinite.(std_L).&isfinite.(dL)
    rY = sum(vy)>=3 ? safe_spearman(std_Y[vy],dY[vy]) : 0.0
    rP = sum(vp)>=3 ? safe_spearman(std_PH[vp],dPH[vp]) : 0.0
    rL = sum(vl)>=3 ? safe_spearman(std_L[vl],dL[vl]) : 0.0

    # Employment comovement moments (59–60).
    # corr(N,GDP) directly from the level-deviation covariance (scale-free).
    # corr(N, GDP/N): labor productivity apl = gdp − n in log-deviations. With
    # Γ_v in LEVEL deviations, convert to log-dev (co)variances by dividing by
    # steady states:  v_n = Γ_nn/N̄²,  v_g = Γ_gg/Ḡ²,  c_ng = Γ_ng/(N̄Ḡ).
    # Then corr(n, g−n) = (c_ng − v_n)/√(v_n·(v_g + v_n − 2c_ng)).
    corr_NGDP_m = xcorr("N","GDP_vol")
    corr_NAPL_m = let iN=get(ei_sub,"N",0), iG=get(ei_sub,"GDP_vol",0)
        if iN==0 || iG==0
            NaN
        else
            Nbar = max(abs(ys[needed_idx[iN]]), 1e-12)
            Gbar = max(abs(ys[needed_idx[iG]]), 1e-12)
            v_n  = max(Γ_v[iN,iN],0.0)/Nbar^2
            v_g  = max(Γ_v[iG,iG],0.0)/Gbar^2
            c_ng = Γ_v[iN,iG]/(Nbar*Gbar)
            den  = sqrt(max(v_n,0.0)*max(v_g + v_n - 2c_ng, 0.0))
            den < 1e-15 ? 0.0 : clamp((c_ng - v_n)/den, -1.0, 1.0)
        end
    end

    # Moment 61 (2026-08-19): std of the goods expenditure share.
    # In the model omega_t = exp(om_g), and Dynare gives Γ in deviations of om_g
    # from its SS log(ombar), so sqrt(Γ) = std(log omega) and the LEVEL std is
    # sqrt(Γ)·ombar. That matches how std_omG is built from the data in
    # Data/build_reallocation_calibration.py (std of the HP cycle of the level).
    # Under Cobb-Douglas omega_t IS the goods expenditure share (FGI 2023 eq. 12),
    # so this moment maps one-for-one onto sigma_omg.
    _ombar = get_param_val(context, "ombar")
    isnan(_ombar) && (_ombar = baseline.d_omG)
    std_omG = let k = get(ei_sub, "om_g", 0)
        k == 0 ? 0.0 : sqrt(max(Γ_v[k,k], 0.0)) * _ombar
    end

    # Moments 62-63 (2026-08-19): the goods-services relative price channel that
    # identifies cl. gap = pi_g - pi_s, so var(gap) = var(pi_g)+var(pi_s)-2cov.
    # Γ_v holds LEVEL deviations of gross inflation, and gross-inflation
    # deviations equal net-inflation deviations to first order, which is the same
    # object as the data's HP cycle of log gross inflation
    # (compute_data_moments.jl:539-542). om_g = log(omega), matching hp(log omega)
    # in the data.
    i_pg = get(ei_sub, "pi_g", 0); i_ps = get(ei_sub, "pi_s", 0); i_og = get(ei_sub, "om_g", 0)
    var_gap = (i_pg == 0 || i_ps == 0) ? 0.0 :
              max(Γ_v[i_pg,i_pg] + Γ_v[i_ps,i_ps] - 2Γ_v[i_pg,i_ps], 0.0)
    std_pigap = sqrt(var_gap)
    corr_pigap_om = (i_pg == 0 || i_ps == 0 || i_og == 0 || var_gap < 1e-20) ? 0.0 :
        let cov_go = Γ_v[i_pg,i_og] - Γ_v[i_ps,i_og],
            d = sqrt(var_gap * max(Γ_v[i_og,i_og], 0.0))
            d < 1e-15 ? 0.0 : clamp(cov_go / d, -1.0, 1.0)
        end

    # ---- 78: average pairwise cross-sectoral output correlation ----------- #
    # The network literature's cross-sectional summary statistic (FSW 2011 Table
    # 7; Atalay 2017 eq. 18). Smooth in theta, unlike the rank correlations it
    # replaces, and governed by exactly the input elasticities this paper is about.
    rbar_YY = let acc = 0.0, np = 0
        for i in 1:nsec, j in (i+1):nsec
            c = xcorr("Y_$(i)", "Y_$(j)")
            isnan(c) && continue
            acc += c; np += 1
        end
        np > 0 ? acc / np : NaN
    end

    # ---- 79-80: goods / services average sectoral employment volatility ---- #
    # Output-share weighted within each group, matching how d_std_Lg and
    # d_std_Ls are built in compute_data_moments.jl section 9.
    std_Lg_m, std_Ls_m = let Yss = baseline.Y_ss
        gi = baseline.goods; si = baseline.services
        wg = sum(Yss[gi]) > 0 ? Yss[gi] ./ sum(Yss[gi]) : fill(1/length(gi), length(gi))
        ws = sum(Yss[si]) > 0 ? Yss[si] ./ sum(Yss[si]) : fill(1/length(si), length(si))
        (sum(wg .* std_L[gi]), sum(ws .* std_L[si]))
    end

    # ---- 81-84: SOE great ratios ------------------------------------------ #
    # Gamma_v holds LEVEL deviations, so convert to log-deviation (co)variances
    # by dividing by steady states — the same conversion corr(N,GDP/N) uses.
    # Real investment is EInv/PI_inv, so var(log I) = v_e + v_p - 2c_ep.
    ratio_stdC_m, ratio_stdI_m, corr_CGDP_m, corr_IGDP_m =
        let iC = get(ei_sub,"C",0), iE = get(ei_sub,"EInv",0),
            iP = get(ei_sub,"PI_inv",0), iG = get(ei_sub,"GDP_vol",0)
            if iC == 0 || iE == 0 || iP == 0 || iG == 0
                (NaN, NaN, NaN, NaN)
            else
                sb(k) = max(abs(ys[needed_idx[k]]), 1e-12)
                v_c = max(Γ_v[iC,iC],0.0)/sb(iC)^2
                v_g = max(Γ_v[iG,iG],0.0)/sb(iG)^2
                v_e = max(Γ_v[iE,iE],0.0)/sb(iE)^2
                v_p = max(Γ_v[iP,iP],0.0)/sb(iP)^2
                c_ep = Γ_v[iE,iP]/(sb(iE)*sb(iP))
                v_i  = max(v_e + v_p - 2c_ep, 0.0)          # real investment
                c_cg = Γ_v[iC,iG]/(sb(iC)*sb(iG))
                c_ig = (Γ_v[iE,iG]/(sb(iE)*sb(iG))) - (Γ_v[iP,iG]/(sb(iP)*sb(iG)))
                sg = sqrt(v_g)
                (sg > 1e-12 ? sqrt(v_c)/sg : NaN,
                 sg > 1e-12 ? sqrt(v_i)/sg : NaN,
                 (v_c*v_g) > 1e-24 ? clamp(c_cg/sqrt(v_c*v_g), -1.0, 1.0) : 0.0,
                 (v_i*v_g) > 1e-24 ? clamp(c_ig/sqrt(v_i*v_g), -1.0, 1.0) : 0.0)
            end
        end

    # Return 77 moments: 12×std_Y + 12×std_PH + 12×std_L + 10×aggregate +
    # 12×corr(Y_i,PH_i) + corr(N,GDP) + corr(N,GDP/N) + std(omG)
    # + std(pi_g-pi_s) + corr(pi_g-pi_s, om_g)
    # + autocorr(GDP) + autocorr(pi) + 12×autocorr(Y_i)      [64-77, 2026-08-21]
    return [std_Y;std_PH;std_L;pstd("GDP_vol");pstd("pi");xcorr("GDP_vol","pi");
            std_TBGDP;pstd("Q");acQ;xcorr("GDP_vol","Q");rY;rP;rL;corr_YPH;
            corr_NGDP_m;corr_NAPL_m;std_omG;std_pigap;corr_pigap_om;
            pac("GDP_vol");pac("pi");[pac("Y_$(i)") for i in 1:nsec];
            rbar_YY; std_Lg_m; std_Ls_m;
            ratio_stdC_m; ratio_stdI_m; corr_CGDP_m; corr_IGDP_m], true
end


# =========================================================================== #
#  11. MAIN ESTIMATION FUNCTION                                               #
# =========================================================================== #

function smm_run(context::Dynare.Context; endo_names_override=nothing)

    endo_names = if endo_names_override !== nothing && length(endo_names_override) >= 400
        endo_names_override
    else
        Dynare.get_endogenous(context.symboltable)
    end
    @printf "  endo_names: %d variables\n" length(endo_names)

    # Fail loudly NOW if params/moments/bounds are out of sync (see §3b).
    validate_smm_setup(data_moments)

    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)

    # Verify the name-based shock mapping resolved (fixes C1). The expected
    # count is DERIVED from active_shock_indices' own name list rather than
    # hardcoded (2026-08-19): the literal "4 + 2*NSEC = 28" went stale the
    # moment eps_omg was added, and the error text then blamed a stale context
    # when the context was in fact correct (591 endogenous, 31 shocks, 29
    # active). A count guard whose expectation is a magic number just moves the
    # staleness from the data to the check.
    let na = length(baseline.active_exo_idx), exo = smm_exo_names(context),
        want = length(active_shock_indices(context, NSEC))
        @printf "  Active shocks (by name): %d of %d exogenous\n" na length(exo)
        if na != want
            error("""
            Active-shock mismatch: build_baseline resolved $(na) but
            active_shock_indices resolves $(want) from the same context
            ($(length(exo)) exogenous). These two must agree — they are separate
            code paths over the same name list. Check that build_baseline and
            active_shock_indices in smm_model_moments.jl use the same set.
            """)
        elseif na == 0
            error("""
            No shocks matched by name against the loaded context's exogenous
            list ($(length(exo)) shocks). The context is almost certainly the
            OLD model. Rebuild it:
              julia --project=. main_SOE_gap.jl     # recompiles the unified .mod
            then re-run the estimation.
            """)
        end
    end

    n_threads_active = Threads.nthreads()
    contexts_th = n_threads_active > 1 ?
        ((@printf "  Spawning %d per-thread contexts...\n" n_threads_active);
         [deepcopy(context) for _ in 1:n_threads_active]) : [context]
    baselines_th = [deepcopy(baseline) for _ in 1:length(contexts_th)]

    W  = build_weighting_matrix(data_moments)
    # The old text here claimed "5× rank correlations"; the actual weight in
    # build_weighting_matrix has been 2.0 since 2026-07-09. Report what the code
    # does, and flag the labour rank moment being switched off (2026-08-21).
    @printf "  Weighting: VA-share sectoral blocks + 2× rank correlations + 3× correlations.\n"
    @printf "  Labour rank correlation (moment 46): weight=%.1f%s\n\n" W[46,46] (
        W[46,46] == 0 ? "  [UNTARGETED via SMM_W_LABOR_RANK]" : "")

    # Initial θ — θ-REDUCTION (the only mode): pin 28 params via SMM_PIN (elasticities +
    # 24 measured sectoral shocks + external rho/sigma + etastar); estimate only the 7 transmission
    # params (θ 1,4,5,6,33,34,35,36), seeded at interior values.
    θ0 = default_theta0(context)
    θ0[2] = parse(Float64, get(ENV, "SMM_EPSY", "0.80"))   # epsY (Atalay eps_Q) — pinned
    θ0[3] = parse(Float64, get(ENV, "SMM_EPSM", "0.20"))   # epsM (Atalay eps_m) — pinned
    # SHOCK_STAGE selects which sectoral shock calibration to load (2026-08-21):
    #   1 (DEFAULT) — sectoral_shock_calibration_stage1.csv, the analytic split
    #   2            — sectoral_shock_calibration.csv, model-inverted
    #
    # DEFAULT CHANGED TO 1 ON 2026-08-21. Stage 2 is retired from the baseline.
    # Head to head on identical data, moments and weights:
    #     Stage 2  obj = 14.752   8 sigmas at LB, 2 at UB — degenerate
    #     Stage 1  obj = 11.462   all sigmas at measured values, 0.0055-0.0387
    # Stage 1 is 22% better AND is the only one of the two that produces a shock
    # vector defensible in print. lambda_A came back at 0.836, near 1, which is
    # the sign the measured sizes were already close in scale.
    # SHOCK_STAGE=2 still works, for the robustness appendix.
    _stage_sel = get(ENV, "SHOCK_STAGE", "1")
    _scf = _stage_sel == "1" ?
        joinpath(DATA_DIR, "sectoral_shock_calibration_stage1.csv") :
        joinpath(DATA_DIR, "sectoral_shock_calibration.csv")
    if _stage_sel == "1" && !isfile(_scf)
        error("""
            SHOCK_STAGE=1 but $(basename(_scf)) does not exist.
            It is written by calibrate_sectoral_shocks.jl as a backup of the
            Stage-1 input, or you can regenerate it directly:
              julia --project=. compute_sectoral_shocks.jl
              cp Data/sectoral_shock_calibration.csv Data/sectoral_shock_calibration_stage1.csv
            """)
    end
    if isfile(_scf)
        _sc = CSV.read(_scf, DataFrame)
        θ0[7:18]  = Float64.(_sc.isigma_tfp_init)   # measured (may exceed old bounds — pinned, not searched)
        θ0[19:30] = Float64.(_sc.sigma_om_init)
        # Measured sectoral TFP persistences (2026-08-24). Stored OUTSIDE θ, like
        # the shock sizes' cross-sectional pattern, and scaled by θ[6].
        if !SECTORAL_RHO
            RHO_A_MEASURED[] = nothing
            @printf "  sectoral TFP persistence: COMMON rho_A = theta[6] (baseline; SMM_SECTORAL_RHO=1 to use the measured vector — tested and rejected, see utils.jl)\n"
        elseif hasproperty(_sc, :rho_A_init) && !all(isnan, Float64.(_sc.rho_A_init))
            _r = Float64.(_sc.rho_A_init)
            any(isnan, _r) && error("""
                rho_A_init in $(basename(_scf)) contains NaN for some sectors.
                A partly-measured persistence vector would silently mix measured
                and default values. Re-run compute_sectoral_shocks.jl.
                """)
            RHO_A_MEASURED[] = _r
            @printf "  sectoral TFP persistence: MEASURED, %.3f-%.3f (mean %.3f), scaled by lambda_rho = theta[6]\n" minimum(_r) maximum(_r) mean(_r)
        else
            RHO_A_MEASURED[] = nothing
            @printf "\n  %s\n" repeat("!", 70)
            @printf "  WARNING: no rho_A_init column in %s.\n" basename(_scf)
            @printf "  Falling back to a SINGLE COMMON TFP persistence, and theta[6] reverts\n"
            @printf "  to meaning rho_A rather than a scale — but its bounds are now [%.2f, %.2f],\n" LB[6] UB[6]
            @printf "  which are SCALE bounds, so the fallback is not a supported configuration.\n"
            @printf "  Fix, in order:\n"
            @printf "    FORCE_RECOMPUTE=1 julia --project=. compute_data_moments.jl\n"
            @printf "    julia --project=. compute_sectoral_shocks.jl\n"
            @printf "  %s\n\n" repeat("!", 70)
        end
    end
    _ecf = joinpath(DATA_DIR, "external_shock_calibration.csv")
    if isfile(_ecf)
        _ec = CSV.read(_ecf, DataFrame); _em = Dict(String.(_ec.param) .=> Float64.(_ec.value))
        haskey(_em, "rho_pvstar")   && (θ0[31] = _em["rho_pvstar"])
        haskey(_em, "sigma_pvstar") && (θ0[32] = _em["sigma_pvstar"])
    end
    θ0[35] = parse(Float64, get(ENV, "SMM_ETASTAR", "1.0"))  # η* export elasticity — now PINNED (paper Table 3 = 1; Feenstra ~1-1.5), was estimated & drifted to ~4.5
    # λ = 1 means "use the measured shock sizes as they are". Seeding at 1 makes
    # the first evaluation identical to the pre-2026-08-20 model, so the run
    # starts from the known objective and can only improve from there.
    θ0[37] = parse(Float64, get(ENV, "SMM_LAMBDA_A",  "1.0"))
    θ0[38] = parse(Float64, get(ENV, "SMM_LAMBDA_OM", "1.0"))

    # NEWLY PINNED (2026-08-21): log(kappaV) and rho_zeta. These MUST be set
    # BEFORE the SMM_PIN snapshot below. For any index NOT in FREE_THETA the
    # objective reads its value out of SMM_PIN — it never sees the later θ0
    # assignments — so a pinned parameter seeded after the snapshot would
    # silently take whatever default_theta0 happened to return. That ordering
    # trap is why these two lines live here and not with the free seeds.
    #
    # The objective is flat in both (|d obj| = 1.8e-4 and 3.4e-4 in the
    # pre-flight scan), so the chosen values are immaterial to the fit; they are
    # interior and round so the paper can report them as calibrated, not
    # estimated.
    #
    # kappaV is not merely weakly identified — it is STRUCTURALLY INERT, and it
    # is worth knowing why before anyone tries to "fix" its identification. The
    # import-pricing FOC (NK_SOE_lev_gap2.mod ~line 924) is
    #     1 - eV + eV*Q*PVstar/PV - kV*(.) + beta*kV*(.) = 0
    # Divide through by eV: the Rotemberg terms enter with weight kV/eV.
    # params_jl.mod sets epsilonV = 1e13, so kV/eV = 1.3e-7 and the equation
    # collapses to PV = Q*PVstar, the law of one price. Any kappaV in any
    # plausible range is multiplied by ~1e-7 and cannot move a single moment.
    # To carry the weight the DOMESTIC block does (epsilon = 10, kappa_i in
    # 8..509) kappaV would have to be ~1e14-1e15.
    #
    # TWO THINGS FOR THE PAPER, both verified against the code 2026-08-21:
    #  (1) The 1e13 in Table 3 is epsilonV, the import-variety elasticity, NOT
    #      kappa_v. eV -> infinity is what delivers the competitive/flexible
    #      import price, so the DESCRIPTION is right and the SYMBOL is wrong.
    #      kappa_v's actual value in params_jl.mod is ~1.3e6, and it is inert.
    #  (2) With eV = 1e13 the model has COMPLETE, INSTANTANEOUS exchange-rate
    #      pass-through. The text describes incomplete pass-through and cites
    #      Romero (2025 JIE) for it. Text and calibration contradict each other;
    #      pick one. If incomplete pass-through is wanted, eV has to come down
    #      to a normal markup (~10) and kappaV then becomes a real parameter.
    θ0[4]  = log(parse(Float64, get(ENV, "SMM_KAPPAV", "1e6")))  # log(kappaV) — PINNED

    SMM_PIN[] = copy(θ0)                  # objective reads the pinned entries from here
    # FREE transmission params — interior seeds (never at a bound):
    θ0[1]  = 1.0        # ilabcosts
    θ0[5]  = 0.5        # rho_om
    # θ[6] = lambda_rho, the scale on the measured sectoral persistences.
    # Seed at 1.0 — the measured values reproduced exactly, which is the nesting
    # point and the natural null. (Was 0.5, when θ[6] was rho_A itself.)
    θ0[6]  = SECTORAL_RHO ? 1.0 : 0.5
    θ0[33] = 0.7        # rho_zeta (free again — see the FREE_THETA note in utils.jl)
    θ0[34] = 0.030      # sigma_zeta (interior under the raised UB of 0.15)
    # θ0[35] (etastar) and θ0[4] (kappaV) pinned above
    θ0[36] = 100.0      # kappaw

    # ---- WARM START (added 2026-08-24) ------------------------------------- #
    #
    # THE BUG THIS FIXES. save_checkpoint's docstring calls smm_checkpoint.csv
    # "the warm-start source" and utils.jl line 457 says it is deliberately left
    # untagged so it can stay that source — but NOTHING EVER READ IT outside
    # report-only mode. Every run, including every point of the epsY and epsM
    # sweeps, cold-started from the hardcoded seeds directly above. Prior work
    # was written to disk and thrown away.
    #
    # That is not cosmetic. CMA-ES on this objective lands in different basins
    # from the same seed, so a cold restart can return a WORSE point than one
    # already on disk — which is exactly what happened on 2026-08-24: a fresh
    # run finished at 8.7043 while the stored theta-hat scored 8.5484 on the
    # same data. It also means the sweep profiles are noisier than they look:
    # each grid point was an independent cold start, so part of the curvature
    # is basin luck rather than the parameter being profiled.
    #
    # ONLY the FREE dimensions are taken from the file. Pinned entries must keep
    # the values SMM_PIN was just snapshotted with, or a stale checkpoint would
    # silently reintroduce old pins (the ordering trap documented above).
    #
    #   SMM_WARM_START=0            cold start (the old behaviour)
    #   SMM_WARM_START=1            read estimation_results/smm_checkpoint.csv
    #   SMM_WARM_START=<path>       read that file instead
    let ws = get(ENV, "SMM_WARM_START", "0")
        if ws != "0"
            path = ws == "1" ? joinpath(ESTIMATION_DIR, "smm_checkpoint.csv") :
                   (isabspath(ws) ? ws : joinpath(ESTIMATION_DIR, ws))
            if !isfile(path)
                @printf "  [warm start] %s not found — cold-starting from the seeds.\n" path
            else
                try
                    _w = CSV.read(path, DataFrame)
                    _d = Dict(String(r.param) => Float64(r.value) for r in eachrow(_w)
                              if !ismissing(r.value) && isfinite(Float64(r.value)))
                    taken = String[]
                    for k in FREE_THETA
                        nm = CSV_PARAM_NAMES[k]
                        if haskey(_d, nm)
                            θ0[k] = clamp(_d[nm], LB[k], UB[k])
                            push!(taken, nm)
                        end
                    end
                    @printf "  [warm start] %d of %d free params seeded from %s\n" length(taken) length(FREE_THETA) basename(path)
                    @printf "               %s\n" join(taken, ", ")
                    @printf "               (pinned entries NOT taken from the file — SMM_PIN already snapshotted.)\n"
                catch e
                    @printf "  [warm start] could not read %s (%s) — cold-starting.\n" basename(path) sprint(showerror, e)
                end
            end
        end
    end

    θ0 = clamp.(θ0, LB, UB)              # pinned dims clamped so scaled space stays [0,1]; objective overrides exactly
    @printf "  θ-reduction: %d transmission params estimated, %d pinned (etastar calibrated to 1).\n" length(FREE_THETA) (N_THETA - length(FREE_THETA))
    @printf "  shock-mix scales free: lambda_A, lambda_om (1.0 = measured sizes; the pinned-mix model is nested at 1,1).\n"
    @printf "  pins: epsY=%s epsM=%.2f | 24 measured sectoral shocks | rho_pvstar=%.3f sigma_pvstar=%.3f\n" (
        _FREE_EPSY ? @sprintf("FREE [%.2f,%.2f], seed %.2f", LB[2], UB[2], θ0[2]) : @sprintf("%.2f", SMM_PIN[][2])
    ) SMM_PIN[][3] SMM_PIN[][31] SMM_PIN[][32]
    _FREE_EPSY && @printf "  NOTE: epsY is free, so every evaluation recomputes the steady state — slower per eval.\n"
    @printf "  pin (structurally inert given epsilonV=1e13): kappaV=%.1e\n" exp(SMM_PIN[][4])
    @printf "  free seeds: sigma_zeta=%.3f ilabcosts=%.2f kappaw=%.0f | pinned etastar=%.2f\n" θ0[34] θ0[1] θ0[36] θ0[35]
    # Shock-calibration provenance: an SMM run on Stage-1 shock sizes is not
    # comparable to one on Stage-2 sizes, and the objective level differs. Say
    # which one this is, in the log, every run.
    let _f = _scf   # the file actually loaded above, not a second hardcoded path
        if isfile(_f)
            _sc = CSV.read(_f, DataFrame)
            _stage = hasproperty(_sc, :stage) ? String(_sc.stage[1]) : "UNKNOWN (pre-2026-08-21 file)"
            _src   = hasproperty(_sc, :moment_source) ? String(_sc.moment_source[1]) : "UNKNOWN"
            @printf "  sectoral shocks: stage=%s | inverted from %s | SHOCK_STAGE=%s (%s)\n" _stage _src _stage_sel basename(_scf)
            if _stage_sel == "2"
                @printf "  NOTE: running on Stage-2 (model-inverted) shock sizes. The baseline is\n"
                @printf "  Stage 1 (obj 11.46 vs 14.75); Stage 2 is kept for the robustness appendix.\n"
            end
        end
    end

    # Pre-flight
    @printf "=== PRE-FLIGHT ===\n"
    ok_resolve, g_test, _, _ = resolve_first_order!(context)
    if !ok_resolve
        error("Model re-solve failed at pre-flight. Run main_SOE_gap.jl first.")
    end
    @printf "  resolve_first_order! OK, g1_1 size=%s\n" string(size(g_test))

    # Safe wrapper: a throw in the moment fn becomes (NaNs, false) so pre-flight
    # reports it cleanly instead of dumping a stack trace.
    _safe_moments(θv) = try
        smm_model_moments(θv, context, baseline, endo_names)
    catch e
        @printf "  pre-flight eval threw: %s\n" sprint(showerror, e)
        (fill(NaN, N_MOMENTS), false)
    end

    m_test, ok_test = _safe_moments(θ0)
    if !ok_test || any(isnan, m_test)
        # Interior free seeds didn't solve — retry at the values main_SOE_gap.jl used
        # (known to solve), keeping the 28 pins.
        @printf "  Interior seeds failed pre-flight; retrying free seeds at known-solving values.\n"
        θ0[1]  = clamp(0.18,       LB[1],  UB[1])    # ilabcosts
        θ0[4]  = clamp(log(4.6e6), LB[4],  UB[4])    # log(kappaV)
        θ0[34] = clamp(0.020,      LB[34], UB[34])   # sigma_zeta (was 0.049; keep the fallback modest post-discipline)
        # θ0[35] (etastar) stays at its pinned value (1.0) — no longer a free seed
        m_test, ok_test = _safe_moments(θ0)
    end
    if !ok_test || any(isnan, m_test)
        @printf "  FAILED: NaN at %s\n" string(findall(isnan, m_test))
        error("Pre-flight failed even from default θ₀. Re-run main_SOE_gap.jl (EXERCISE=0) to rebuild the context, then check the model compiles.")
    end
    obj_test    = dot(data_moments .- m_test, W * (data_moments .- m_test))
    m_test_copy = copy(m_test)   # used as fallback for best_moments below
    @printf "  obj(θ₀) = %.6f\n" obj_test

    ψ0 = data_moments .- m_test
    # Driven by MOMENT_BLOCKS (2026-08-19) rather than hardcoded 1:12 … 59:60.
    # Adding moment 61 (std(omG)) left it inside the TOTAL but outside every
    # printed slice, so the blocks silently stopped summing to the objective.
    # Now any future moment addition is picked up automatically.
    @printf "  Decomp: %s  (Σ=%.3f, obj=%.3f)\n" join(
        [@sprintf("%s=%.3f", first(split(lbl)), _blk0(ψ0, W, rng))
         for (rng, lbl) in MOMENT_BLOCKS], " ") sum(
        _blk0(ψ0, W, rng) for (rng, _) in MOMENT_BLOCKS) obj_test

    print_param_table(θ0)
    print_fit_table(data_moments, m_test, W)

    # ---- SENSITIVITY GUARD (2026-08-19) ---------------------------------- #
    # Refuse to start if the objective does not respond to the free parameters.
    #
    # From 2026-07 until 2026-08-19 it did not: resolve_first_order!'s fast path
    # called Dynare's solver, which reads context.work.params, while set_param!
    # writes to _SMM_PARAMS. Every evaluation therefore returned the decision
    # rule at the ORIGINAL parameters and the objective was EXACTLY constant in
    # all 7 free dims. CMA-ES duly "converged" after ~126 evaluations with θ
    # unchanged, every Jacobian column was zero (std err 0.0, t = Inf), and two
    # separate estimation runs produced nothing — with no error raised.
    #
    # A silent no-op optimiser is the worst possible failure mode, so this is
    # now a hard stop rather than a warning.
    let bad = String[]
        # FREE_THETA, not _FREE: the local alias is not bound until the CMA-ES
        # section further down.
        for p in FREE_THETA
            θt = copy(θ0)
            θt[p] = θ0[p] + 0.25*(UB[p] - θ0[p]) + 1e-6   # a quarter of the way up
            mt, okt = smm_model_moments(θt, context, baseline, endo_names)
            if !okt || any(isnan, mt)
                @printf "  [sensitivity] θ[%d] %-18s : model FAILED at the test point (skipped)\n" p PARAM_LABELS[p]
                continue
            end
            ψt = data_moments .- mt
            Δ  = abs(dot(ψt, W*ψt) - obj_test)
            @printf "  [sensitivity] θ[%d] %-18s : |Δobj| = %.3e%s\n" p PARAM_LABELS[p] Δ (Δ < 1e-10 ? "   <<<< FLAT" : "")
            Δ < 1e-10 && push!(bad, PARAM_LABELS[p])
        end
        # A flat direction is normally a hard error — it means CMA-ES would
        # "converge" instantly having estimated nothing. But there is a second
        # cause that is NOT a modelling problem: resolve_first_order! silently
        # falling back to the stale decision rule, which makes every parameter
        # look flat at that theta. Look for "[WARN] ... STALE decision rule"
        # just above; if it is there, the model is fragile at this point in the
        # parameter space, not the parameter dead.
        #
        # In a SWEEP that distinction matters operationally: one fragile grid
        # point should not abort the remaining points. SMM_STRICT_SENS=0
        # downgrades this to a warning so the sweep records the failure and
        # moves on. Keep it STRICT for a single production run.
        _strict_sens = get(ENV, "SMM_STRICT_SENS", "1") != "0"
        if !isempty(bad) && !_strict_sens
            @printf "\n  %s\n" repeat("!", 70)
            @printf "  FLAT DIRECTIONS: %s\n" join(bad, ", ")
            @printf "  SMM_STRICT_SENS=0, so continuing anyway. If a STALE-decision-rule\n"
            @printf "  warning appears above, this is solver fragility at this theta, not a\n"
            @printf "  dead parameter — the estimate from this run is NOT trustworthy.\n"
            @printf "  %s\n\n" repeat("!", 70)
            empty!(bad)
        end
        isempty(bad) || error("""

            OBJECTIVE DOES NOT RESPOND TO: $(join(bad, ", "))

            These parameters are in the free set but moving them a quarter of the
            way to their upper bound leaves the objective unchanged to 1e-10.
            Estimating them is meaningless — CMA-ES will "converge" immediately
            at θ₀ and report success.

            Most likely causes:
              1. θ is not reaching the solver. Check that set_param! writes
                 somewhere resolve_first_order! actually reads
                 (_SMM_PARAMS vs context.work.params).
              2. resolve_first_order! is returning the stale decision rule —
                 look for "[WARN] ... STALE decision rule" above.
              3. The parameter is genuinely dead in the .mod (as ilabcosts was
                 until 2026-08-19: it appeared only in a reporting variable).

            Run probe_free_dims.jl to see the full-range response of each.
            """)
    end
    @printf "=== PRE-FLIGHT PASSED ===\n\n"

    # Report-only mode: print the fit of the current θ (warm start / checkpoint)
    # and exit without optimizing. Usage: SMM_REPORT_ONLY=1 julia run_smm_estimation.jl
    if get(ENV, "SMM_REPORT_ONLY", "0") in ("1", "true")
        @printf "  SMM_REPORT_ONLY set — fit report printed above, skipping CMA-ES.\n"
        @printf "  (smm_estimates.csv / checkpoint NOT modified.)\n\n"

        # ── ENGINE CROSS-CHECK (2026-07) ─────────────────────────────────────
        # Score the IN-PROCESS moment engine at the CHECKPOINT θ, so its
        # per-block objective can be compared block-by-block against the
        # Dynare-subprocess table that main_SOE_gap.jl prints at the SAME θ.
        # The two solve paths share the .mod and Dynare's perturbation solver but
        # NOT the steady state: smm_model_moments normalises by the SS that
        # recompute_ss_cached! produces, whereas main normalises by the Dynare
        # subprocess SS. If the checkpoint's recorded obj (~95) and main's
        # recomputed obj (~40) disagree, this localises the gap:
        #   • all 7 blocks off by a common ratio ⇒ global SS-scale difference;
        #   • one block exploding ⇒ a specific variable's SS diverges (watch the
        #     tiny sectors: Public Admin SS(Y)=0.0022 amplifies any SS error).
        _ckpt = joinpath(ESTIMATION_DIR, "smm_checkpoint.csv")
        if isfile(_ckpt)
            _cdf = CSV.read(_ckpt, DataFrame)
            θ_ck = Float64.(_cdf.value)
            if length(θ_ck) == N_THETA
                m_ck, ok_ck = smm_model_moments(θ_ck, context, baseline, endo_names)
                if ok_ck && !any(isnan, m_ck)
                    ψc   = data_moments .- m_ck
                    objc = dot(ψc, W * ψc)
                    _blk(a,b) = dot(ψc[a:b], W[a:b,a:b] * ψc[a:b])
                    @printf "\n  ── ENGINE CROSS-CHECK: in-process engine at CHECKPOINT θ ──\n"
                    @printf "  obj(in-process, checkpoint θ) = %.4f   [checkpoint recorded %.4f]\n" objc _cdf.obj_hat[1]
                    @printf "  Per-block: Y=%.3f PH=%.3f L=%.3f Agg=%.3f Rank=%.3f CorrYP=%.3f NLab=%.3f\n" _blk(1,12) _blk(13,24) _blk(25,36) _blk(37,43) _blk(44,46) _blk(47,58) _blk(59,60)
                    @printf "  Compare these 7 blocks to main_SOE_gap.jl's Dynare table (same θ):\n"
                    @printf "  uniform ratio ⇒ global SS scale; one block exploding ⇒ that variable's SS diverges.\n\n"
                else
                    @printf "  [cross-check] in-process eval failed at checkpoint θ (ok=%s, NaN=%s)\n" ok_ck string(any(isnan, m_ck))
                end
            end
        end
        return θ0, obj_test, m_test
    end

    # CMA-ES
    max_evals = 300_000
    _FREE = FREE_THETA   # single source of truth in utils.jl — shared with smm_inference.jl
    @printf "--- CMA-ES ---\n"
    @printf "  %d FREE params (of %d; %d pinned) | %d moments | max %d evals | %d threads\n" length(_FREE) N_THETA (N_THETA-length(_FREE)) N_MOMENTS max_evals n_threads_active
    @printf "  free: %s\n" join(PARAM_LABELS[_FREE], ", ")
    @printf "  %-6s  %-10s  %-70s  %-7s  %-8s  %-8s\n" "eval" "best_obj" "[Y  PH  L  Agg  Rank  CorrYP  NLab  Om  ACagg  ACsec]" "fail" "ms/eval" "Klein%"
    @printf "  %s\n" repeat("-",90)

    best_θ        = Ref(clamp.(θ0, LB, UB))   # stored in original parameter space
    best_obj      = Ref(obj_test)              # initialise with pre-flight result
    best_moments  = Ref(m_test_copy)           # moments at best θ — avoids re-evaluation
    best_lock     = ReentrantLock()
    print_lock    = ReentrantLock()
    fail_count    = Threads.Atomic{Int}(0)
    eval_count    = Threads.Atomic{Int}(0)
    last_printed  = Threads.Atomic{Int}(0)   # last eval_count at which we printed
    t_start       = Ref(time())
    t_last_print  = Ref(time())
    PRINT_EVERY   = 50    # print every N evaluations (any thread can trigger)

    # Machine-readable progress log (2026-07-08): one row per PRINT_EVERY block.
    # Plottable objective trajectory + fit decomposition; survives node failure
    # alongside smm_checkpoint.csv. Header written fresh at every run start.
    mkpath(ESTIMATION_DIR)
    progress_log = joinpath(ESTIMATION_DIR, "smm_progress_log.csv")
    open(progress_log, "w") do io
        println(io, "timestamp,elapsed_s,evals,best_obj,decomp_Y,decomp_PH,decomp_L,decomp_Agg,decomp_Rank,decomp_CorrYP,decomp_NLab,decomp_Om,decomp_ACagg,decomp_ACsec,fails,ms_per_eval,klein_hit_pct")
    end

    # Wall-clock self-limit: set SMM_MAX_HOURS a bit under the SLURM --time so the
    # run stops itself and saves, rather than being SIGKILLed mid-write.
    max_seconds = let h = tryparse(Float64, get(ENV, "SMM_MAX_HOURS", ""))
        (h === nothing || h <= 0) ? Inf : h * 3600
    end
    max_seconds < Inf && @printf "  Wall-clock self-limit: %.2f h (SMM_MAX_HOURS)\n" (max_seconds/3600)

    _tid_cma() = min(Threads.threadid(), length(contexts_th))

    # ---- STAGNATION STOP (2026-08-24) ------------------------------------- #
    # ftol and xtol above are deliberately near-disabled, because premature
    # convergence used to be the main failure mode. The cost is the opposite
    # failure: once CMA-ES reaches a flat basin it keeps drawing candidates for
    # thousands of evaluations that change nothing. Observed in the epsM sweep —
    # the objective was frozen at 9.4462 from evaluation 1800 to 3000 and still
    # running, and the epsY sweep's runtimes ranged 356s to 2616s for the same
    # work, entirely because some points happened to trip xtol and others did not.
    #
    # A spread-based criterion (ftol/xtol) asks "are the candidates close
    # together?". What actually matters is "has the BEST value improved
    # lately?". Track that directly: stop when the best objective has not
    # improved by more than SMM_STAG_TOL (relative) over the last
    # SMM_STAG_EVALS evaluations. Throwing SMMTimeout reuses the existing
    # graceful-unwind path, so best_θ is saved exactly as on a wall-clock stop.
    stag_evals = parse(Int,     get(ENV, "SMM_STAG_EVALS", "600"))
    stag_tol   = parse(Float64, get(ENV, "SMM_STAG_TOL",   "1e-6"))
    stag_ref   = Ref(Inf)      # best objective at the last improvement
    stag_at    = Ref(0)        # eval count at the last improvement
    stagnated  = Ref(false)

    obj_fn = θ_sc -> begin
        # Graceful wall-clock stop — throw so CMA-ES unwinds; caller saves best θ.
        (max_seconds < Inf && (time() - t_start[]) > max_seconds) && throw(SMMTimeout())

        # Graceful stagnation stop, same unwind path.
        if stag_evals > 0
            b = best_obj[]
            if isfinite(b) && (stag_ref[] - b) > stag_tol * max(abs(stag_ref[]), 1.0)
                stag_ref[] = b; stag_at[] = eval_count[]
            elseif eval_count[] - stag_at[] > stag_evals && stag_at[] > 0
                if !stagnated[]
                    stagnated[] = true
                    @printf "\n  STAGNATION STOP: no improvement > %.0e in %d evaluations (best %.6f).\n" stag_tol stag_evals b
                    @printf "  SMM_STAG_EVALS=0 disables; raise SMM_STAG_EVALS for a longer patience.\n\n"
                    flush(stdout)
                end
                throw(SMMTimeout())
            end
        end

        tid   = _tid_cma()
        θ = copy(SMM_PIN[])           # 36-vec: 29 pinned entries; 7 free set below
        θ[_FREE] = LB[_FREE] .+ θ_sc .* span   # unscale the 7 free dims [0,1] → original
        # A single bad evaluation must NEVER take down a multi-hour run. Any
        # exception inside the moment computation is converted to the failure
        # penalty (1e8); only SMMTimeout is allowed to propagate.
        local m
        obj = try
            mm, ok = smm_model_moments(θ, contexts_th[tid], baselines_th[tid], endo_names)
            m = mm
            (!ok || any(isnan, mm)) ? 1e8 : dot(data_moments.-mm, W*(data_moments.-mm))
        catch err
            err isa SMMTimeout && rethrow(err)
            1e8
        end

        n = Threads.atomic_add!(eval_count, 1) + 1   # new count (1-based)
        obj >= 1e7 && Threads.atomic_add!(fail_count, 1)

        if isfinite(obj) && obj < best_obj[]
            lock(best_lock) do
                if isfinite(obj) && obj < best_obj[]
                    best_obj[]      = obj
                    best_θ[]        = copy(θ)   # θ in original parameter space
                    best_moments[]  = copy(m)   # save moments — avoids re-evaluation bug
                    # Live checkpoint on EVERY improvement (2026-07-10). The old
                    # 60s throttle lost the final best of an interrupted run:
                    # an improvement inside the throttle window was never saved
                    # and Ctrl-C bypasses the end-of-run save (local run 2026-07-09:
                    # best 9839 lost, checkpoint kept the earlier 13083 θ).
                    # Improvements are rare and the atomic 1KB write costs <1ms,
                    # so throttling buys nothing.
                    save_checkpoint(best_θ[], best_obj[])
                end
            end
        end

        # Print when we cross the next multiple of PRINT_EVERY.
        # Any thread can trigger; trylock prevents duplicate/garbled output.
        if n - last_printed[] >= PRINT_EVERY
            if trylock(print_lock)
                try
                    if n - last_printed[] >= PRINT_EVERY   # re-check inside lock
                        t_now  = time()
                        dt     = t_now - t_last_print[]
                        ms     = dt > 0 ? 1000.0 * dt / max(n - last_printed[], 1) : 0.0
                        tot    = _KLEIN_HITS[] + _KLEIN_MISSES[]
                        kpct   = tot > 0 ? round(Int, 100*_KLEIN_HITS[]/tot) : 0
                        b_obj  = best_obj[]
                        # Snapshot best moments for decomposition (brief lock, no alloc in hot path)
                        b_mom  = lock(best_lock) do; copy(best_moments[]); end
                        ψ_now  = data_moments .- b_mom
                        # MOMENT_BLOCKS-driven (2026-08-19) — see the pre-flight
                        # note above. dOm is the goods expenditure share.
                        dY  = _blk0(ψ_now, W, MOMENT_BLOCKS[1][1])
                        dPH = _blk0(ψ_now, W, MOMENT_BLOCKS[2][1])
                        dL  = _blk0(ψ_now, W, MOMENT_BLOCKS[3][1])
                        dAg = _blk0(ψ_now, W, MOMENT_BLOCKS[4][1])
                        dRk = _blk0(ψ_now, W, MOMENT_BLOCKS[5][1])
                        dCY = _blk0(ψ_now, W, MOMENT_BLOCKS[6][1])
                        dNL = _blk0(ψ_now, W, MOMENT_BLOCKS[7][1])
                        dOm = length(MOMENT_BLOCKS) >= 8 ?
                              _blk0(ψ_now, W, MOMENT_BLOCKS[8][1]) : 0.0
                        # 10-11: the persistence blocks (2026-08-21). Guarded on
                        # length so an older MOMENT_BLOCKS still prints.
                        dACa = length(MOMENT_BLOCKS) >= 10 ?
                               _blk0(ψ_now, W, MOMENT_BLOCKS[10][1]) : 0.0
                        dACs = length(MOMENT_BLOCKS) >= 11 ?
                               _blk0(ψ_now, W, MOMENT_BLOCKS[11][1]) : 0.0
                        @printf "  %-6d  %-10.4f  [Y=%.2f PH=%.2f L=%.2f Agg=%.2f Rk=%.2f CY=%.2f NL=%.2f Om=%.2f ACa=%.2f ACs=%.2f]  fail=%-5d  %.1fms  Klein=%d%%\n" n b_obj dY dPH dL dAg dRk dCY dNL dOm dACa dACs fail_count[] ms kpct
                        flush(stdout)
                        # Append to the machine-readable progress log (best effort:
                        # a full disk or NFS hiccup must never kill the run).
                        try
                            open(progress_log, "a") do io
                                @printf io "%s,%.1f,%d,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.1f,%d\n" Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS") (t_now - t_start[]) n b_obj dY dPH dL dAg dRk dCY dNL dOm dACa dACs fail_count[] ms kpct
                            end
                        catch; end
                        last_printed[] = n
                        t_last_print[] = t_now
                    end
                finally
                    unlock(print_lock)
                end
            end
        end
        obj
    end

    # Per-parameter sigma: (UB-LB)/6 so each parameter gets its own initial
    # search radius proportional to its feasible range.
    # This is critical: without it, CMA-ES uses the same step size for
    # log(kappaV) (range 12) and isigma_tfp_i (range 0.10), wasting
    # hundreds of evaluations on infeasible or uninformative candidates.
    # CMAEvolutionStrategy.jl takes scalar sigma — pass the mean, but
    # pre-scale the parameter space so all dimensions have unit range.
    insigma = 0.25   # initial step for the 7-dim reduced search (larger avoids early stall)

    # Rescale θ to [0,1] so CMA-ES works in a unit hypercube.
    # The objective wrapper maps back to the original scale.
    # TRUE 7-dim search: scale only the 7 free transmission dims to [0,1].
    span  = UB[_FREE] .- LB[_FREE]
    θ0_sc = (clamp.(θ0[_FREE], LB[_FREE], UB[_FREE]) .- LB[_FREE]) ./ span
    LB_sc = zeros(length(_FREE))
    UB_sc = ones(length(_FREE))

    # The optimiser is wrapped so that a wall-clock timeout OR any unexpected
    # internal failure does not discard hours of search — best_θ is always the
    # running best, and is saved below regardless of how minimize() exits.
    try
        CMAEvolutionStrategy.minimize(
            obj_fn, θ0_sc, insigma;
            lower = LB_sc,
            upper = UB_sc,
            maxiter = max_evals,
            ftol  = 1e-10,   # effectively disabled: HP-filtered obj is smooth,
                              # premature convergence was the main stagnation cause
            xtol  = 1e-8,    # stop only when parameter changes are genuinely tiny
            seed  = 42,
            verbosity = 0,
            multi_threading = n_threads_active > 1)
    catch err
        if _is_timeout(err)
            @printf "\n  %s — stopping CMA-ES; best θ retained.\n" (
                stagnated[] ? "Stagnation stop" : "Wall-clock self-limit reached")
        else
            @printf "\n  [warn] CMA-ES ended early: %s\n  Proceeding with best θ found so far.\n" sprint(showerror, err)
        end
    end
    # Ensure the very latest best is on disk before the results section.
    save_checkpoint(best_θ[], best_obj[])

    θ_hat = clamp.(best_θ[], LB, UB)
    total_t = time() - t_start[]
    tot_kl  = _KLEIN_HITS[] + _KLEIN_MISSES[]
    @printf "\nDone: %d evals | %.1fs | %.1fms/eval | Klein cache=%d%%\n" eval_count[] total_t 1000*total_t/max(1,eval_count[]) round(Int,100*_KLEIN_HITS[]/max(1,tot_kl))

    # Solver-path accounting (2026-08-19). Diagnostics are rate-limited during
    # the run so the log stays readable; the totals must still be reported,
    # because a silent solver fallback is exactly what let the flat-objective
    # bug survive two full estimation runs.
    @printf "Solver: %d Klein aborts | %d STALE-decision-rule evaluations" _KLEIN_ABORTS[] _STALE_DR_USES[]
    if _STALE_DR_USES[] > 0
        @printf "  (%.0f%% of evals wasted)\n" 100*_STALE_DR_USES[]/max(1, eval_count[])
        @printf "  A stale evaluation returns the θ₀ decision rule, so its objective is\n"
        @printf "  ~obj(θ₀) regardless of θ. Those candidates are REJECTED by the optimiser,\n"
        @printf "  so best_θ is always from a genuine evaluation — but the work is wasted.\n"
        @printf "  Cause: Dynare's solver does not survive deepcopy / is not thread-safe, so\n"
        @printf "  worker threads fall through. Re-run with --threads=1 to confirm the result.\n"
    else
        @printf "  — all evaluations used the real solver.\n"
    end
    @printf "\n"

    # Use moments cached at best θ — re-evaluating on the main context would give a
    # different result because context has cold SS state (frozen at θ₀ from pre-flight),
    # while the best was found on a per-thread context with warm SS cache.
    obj_hat = best_obj[]
    m_hat   = best_moments[]
    ψ_hat   = data_moments .- m_hat

    # Results table
    @printf "\n%s\n  SMM RESULTS\n%s\n\n" repeat("=",60) repeat("=",60)
    @printf "  Objective: %.6f\n\n" obj_hat
    @printf "  %-6s  %-18s  %10s  %10s  %s\n" "Idx" "Parameter" "Initial" "Estimate" "Bound?"
    @printf "  %s\n" repeat("-",60)
    for k in 1:N_THETA
        span_k = UB[k] - LB[k]
        flag = θ_hat[k] <= LB[k] + 0.02*span_k ? "<< at LB" :
               θ_hat[k] >= UB[k] - 0.02*span_k ? ">> at UB" : ""
        if k==4
            @printf "  %3d  %-18s  %10.4f  %10.4f  [log] %s\n" k PARAM_LABELS[k] θ0[k] θ_hat[k] flag
            @printf "  %3s  %-18s  %10.2e  %10.2e  [level]\n" "--" "kappaV" exp(θ0[k]) exp(θ_hat[k])
        else
            @printf "  %3d  %-18s  %10.4f  %10.4f  %s\n" k PARAM_LABELS[k] θ0[k] θ_hat[k] flag
        end
    end
    print_fit_table(data_moments, m_hat, W)

    # Save (atomic writes — a kill mid-write can't corrupt these files).
    # All estimation outputs go to ESTIMATION_DIR (julia_dynare/estimation_results).
    df_res = DataFrame(param=vcat(PARAM_LABELS,fill("",N_MOMENTS-N_THETA)),
                        theta=vcat(θ_hat,fill(NaN,N_MOMENTS-N_THETA)),
                        moment=MOMENT_NAMES, data=data_moments, model=m_hat, diff=ψ_hat)
    atomic_write_csv(joinpath(ESTIMATION_DIR, tagged("smm_results.csv")), df_res)

    df_est = DataFrame(
        param=CSV_PARAM_NAMES,
        value=θ_hat, obj_hat=vcat([obj_hat],fill(NaN,N_THETA-1)))
    atomic_write_csv(joinpath(ESTIMATION_DIR, tagged("smm_estimates.csv")), df_est)
    save_checkpoint(θ_hat, obj_hat)   # checkpoint CSV + best_sol.txt + min_loss.txt

    @printf "  Results: %s\n  Estimates: %s\n\n" joinpath(ESTIMATION_DIR, tagged("smm_results.csv")) joinpath(ESTIMATION_DIR, tagged("smm_estimates.csv"))

    # INFERENCE IS DELIBERATELY NOT RUN (2026-08-21, Agustín's decision).
    #
    # The standard errors this used to produce were not reportable and were
    # actively misleading:
    #   * the moment Jacobian G at θ̂ is rank-deficient, so diag(V) came out at
    #     ~0 and every t-statistic was Inf (see probe_free_dims.jl);
    #   * without Data/moment_cov.csv the code fell back on delta-method
    #     DIAGONAL variances for S, which ignores the (large) correlation
    #     between 60 moments built from the same 60 quarters of data.
    # Numbers that look like publication SEs but are neither are worse than no
    # numbers at all, so the call site is removed rather than silenced.
    #
    # smm_inference.jl is kept on disk, unwired, so the machinery survives. To
    # bring it back: supply a moving-block-bootstrap Data/moment_cov.csv, fix
    # the rank deficiency (free parameters that actually move the moments), and
    # restore the call here.

    @printf "  Re-run main_SOE_gap.jl (EXERCISE=0) to apply estimates.\n\n"

    return θ_hat, obj_hat, m_hat
end
