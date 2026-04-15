"""
smm_estimation.jl
=================
Simulated Method of Moments (SMM) estimation for the NK-IOSOE Chile model.

STRATEGY
  The model is solved at first order by Dynare.jl (NK_SOE_lev_gap2.mod).
  For each candidate θ, smm_model_moments() re-solves the linearised system
  (via Dynare.compute_first_order_solution!) WITHOUT recompiling the model.
  Theoretical unconditional second moments are derived analytically from the
  state-space representation via the discrete Lyapunov equation.

PREREQUISITES
  1. Run compute_data_moments.jl at least ONCE to produce data_moments_chile.mat.
  2. main_SOE_gap.jl handles the initial Dynare compile; we import the context.

ESTIMATED PARAMETERS  θ (23-element vector)
  [1]     ilabcosts          aggregate labour adjustment cost (inverse)
  [2]     epsY               elast. of subst. in production
  [3]     epsM               elast. of subst. between materials
  [4]     log(kappaV)        log of import price adj. cost
  [5]     rho_om             AR persistence, goods-services shock
  [6]     sigma_om           std dev, goods-services shock
  [7]     rho_A              AR persistence, sectoral TFP shocks (common)
  [8-19]  isigma_tfp_1,...12 std dev of each sector's TFP shock
  [20]    rho_pvstar         AR persistence, import price shock
  [21]    sigma_pvstar       std dev of import price shock
  [22]    rho_xi             AR persistence, preference shock
  [23]    sigma_xi           std dev of preference shock

Usage:
  julia --project=. smm_estimation.jl
or from REPL:
  include("smm_estimation.jl")
  # Then call:  smm_run(context)   where context comes from main_SOE_gap.jl
"""

using LinearAlgebra, Statistics, StatsBase, Printf
using NLsolve, CSV, DataFrames, Dynare
using CMAEvolutionStrategy

SCRIPT_DIR  = @__DIR__
REPO_ROOT   = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR    = joinpath(REPO_ROOT, "Data")

# NOTE: steady_ntwsoe_system.jl, steady_ntwsoe.jl, utils.jl, and
# smm_model_moments.jl are included by run_smm_estimation.jl BEFORE this
# file is included, to avoid double-include const-redefinition issues.
# When using smm_estimation.jl standalone, include those files first.

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
    d_std_GDP    = agg_dict["std_GDP"]
    d_std_pi     = agg_dict["std_pi"]
    d_corr_GDPpi = agg_dict["corr_GDPpi"]
    d_omG        = agg_dict["omG"]
    d_std_Q      = agg_dict["std_Q"]
    d_autocorr_Q = agg_dict["autocorr_Q"]
    d_corr_GDPQ  = agg_dict["corr_GDPQ"]
    d_TBGDP      = agg_dict["TBGDP"]
    @printf "  Loaded sectoral_moments.csv + aggregate_moments.csv\n\n"
else
    @printf "  CSV moment files not found in %s\n" DATA_DIR
    @printf "  Run compute_data_moments.jl first to generate them.\n\n"
    y_d          = fill(0.04, NSEC)
    p_d          = fill(0.02, NSEC)
    l_d          = fill(0.03, NSEC)
    d_std_GDP    = 0.0215; d_std_pi     = 0.0041;  d_corr_GDPpi = -0.15
    d_omG        = 0.57;   d_std_Q      = 0.0520;  d_autocorr_Q = 0.75
    d_corr_GDPQ  = -0.15;  d_TBGDP      = -0.02
end

# Full 46-element data moment vector
data_moments = [
    y_d;           # 1-12   std(Y_i)
    p_d;           # 13-24  std(PH_i)
    l_d;           # 25-36  std(L_i)
    d_std_GDP;     # 37
    d_std_pi;      # 38
    d_corr_GDPpi;  # 39
    d_omG;         # 40   SS target (passive)
    d_std_Q;       # 41
    d_autocorr_Q;  # 42   identifies rho_pvstar
    d_corr_GDPQ;   # 43   over-identifies pvstar shock
    1.0;           # 44   rank corr output: perfect match target
    1.0;           # 45   rank corr prices
    1.0;           # 46   rank corr labor
]
@assert length(data_moments) == 46 "Expected 46 data moments, got $(length(data_moments))"


# =========================================================================== #
#  2.  LOAD INITIAL MODEL CONTEXT AND BUILD BASELINE                          #
# =========================================================================== #

"""
    build_baseline(context, endo_names, data_moments) -> NamedTuple

Assemble the fixed calibration objects passed into smm_model_moments on
every iteration.  Mirrors the `baseline` struct in MATLAB smm_estimation.m.
"""
function build_baseline(context::Dynare.Context,
                         endo_names::Vector{String},
                         d_std_Y, d_std_PH, d_std_L,
                         d_TBGDP, d_omG)
    # Load calibration vectors from _SMM_PARAMS (robust against uninitialized work.params)
    function pvec(nm)
        idx = param_idx(context, nm)
        idx === nothing && return NaN
        p = _SMM_PARAMS_READY[] ? _SMM_PARAMS[] : _load_smm_params!(context)
        (isempty(p) || idx > length(p)) ? NaN : p[idx]
    end

    nsec  = NSEC
    # Read sectoral parameter vectors from Dynare context
    modalpha   = [pvec("alpha_$(i)")   for i in 1:nsec]
    modalphaV  = [pvec("alphaV_$(i)")  for i in 1:nsec]
    modvarrho  = [pvec("varrho_$(i)")  for i in 1:nsec]
    modgammag  = [pvec("gammag_$(i)")  for i in 1:nsec]
    modgammas  = [pvec("gammas_$(i)")  for i in 1:nsec]
    modchiX    = [pvec("chiX_$(i)")    for i in 1:nsec]
    modkappa   = [pvec("kappa_$(i)")   for i in 1:nsec]
    modbeta    = [pvec("beta_$(i)_$(j)") for i in 1:nsec, j in 1:nsec]

    ss_vec = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let idx = findfirst(==(nm), endo_names)
                     idx === nothing ? 1.0 : ss_vec[idx]
                 end

    Y_ss_vec = [get_ss("Y_$(i)") for i in 1:nsec]

    return (
        nsec        = nsec,
        goods       = GOODS,
        services    = SERVICES,
        Y_ss        = Y_ss_vec,
        GDP_ss      = get_ss("GDP"),
        ombar_val   = let v = pvec("ombar"); isnan(v) ? d_omG : v; end,
        tb_target   = d_TBGDP,
        modalpha    = modalpha,
        modalphaV   = modalphaV,
        modbeta     = modbeta,
        modgammag   = modgammag,
        modgammas   = modgammas,
        modvarrho   = modvarrho,
        modchiX     = modchiX,
        modkappa    = modkappa,
        gamma_val   = pvec("gamma"),
        psi_val     = pvec("psi"),
        chi_val     = 1.0,
        epsilon_val = pvec("epsilon"),
        beta_val    = pvec("beta"),
        PVstar_ss   = pvec("PVstar_ss"),
        sigmaH_val  = pvec("sigmaH"),
        etastar_val = pvec("etastar"),
        omegaX_val  = pvec("omegaX"),
        ystar_ss_val= pvec("ystar_ss"),
        data_std_Y  = d_std_Y,
        data_std_PH = d_std_PH,
        data_std_L  = d_std_L,
    )
end


# =========================================================================== #
#  3.  WEIGHTING MATRIX                                                        #
#                                                                              #
#  Diagonal relative-scale: W[i,i] = 1 / max(|d_i|, floor)^2                #
#  Makes every moment contribute equally in percentage terms.                 #
#  Downweights structurally-mismatched moments (corr(GDP,pi)) and             #
#  near-zero moments that would otherwise dominate mechanically.              #
# =========================================================================== #

function build_weighting_matrix(data_moments::Vector{<:Real})
    floor_w = 0.01
    w_diag  = 1.0 ./ max.(abs.(data_moments), floor_w).^2
    W = Diagonal(w_diag)   # use Diagonal for efficiency; convert to Matrix below

    # Cap weights for small-magnitude moments to prevent mechanical dominance
    # (equivalent to MATLAB Fix 3: |d_i| < 0.02 → 0.25× natural weight)
    W_mat = collect(Matrix(W))
    for i in eachindex(data_moments)
        if abs(data_moments[i]) < 0.02
            W_mat[i, i] *= 0.25
        end
    end

    # Fix 1: corr(GDP,pi) [39] — structurally cannot match with supply shocks only
    W_mat[39, 39] *= 0.02
    # Fix passive moments
    W_mat[40, 40] *= 0.10   # mean omG: imposed by SS calibration
    W_mat[43, 43] *= 0.50   # corr(GDP,Q): over-identified
    W_mat[44, 44] *= 0.50   # rank corr output
    W_mat[45, 45] *= 0.50   # rank corr prices
    W_mat[46, 46] *= 0.50   # rank corr labor

    return W_mat
end


# =========================================================================== #
#  4.  PARAMETER BOUNDS AND INITIAL VALUES                                     #
# =========================================================================== #

const PARAM_LABELS = vcat(
    ["ilabcosts", "epsY", "epsM", "log(kappaV)", "rho_om", "sigma_om", "rho_A"],
    ["isigma_tfp_$(i)" for i in 1:12],
    ["rho_pvstar", "sigma_pvstar", "rho_xi", "sigma_xi"],
)
const N_THETA = length(PARAM_LABELS)  # 23

const LB = [1e-3; 0.10; 0.01; log(1e3);  -0.99; 1e-5; -0.99;
            fill(1e-5, 12); 0.50;  0.005; 0.00;  0.0  ]
const UB = [100.0; 3.00; 1.50; log(1e16); 0.99;  0.50; 0.99;
            fill(0.50, 12); 0.99;  0.50;  0.99;  0.20 ]


"""
    default_theta0(context) -> Vector{Float64}

Build the initial parameter vector from the current Dynare context.
"""
function default_theta0(context::Dynare.Context)
    pv(nm) = let idx = param_idx(context, nm)
                 if idx === nothing; 0.0
                 else
                     p = _SMM_PARAMS_READY[] ? _SMM_PARAMS[] : _load_smm_params!(context)
                     (isempty(p) || idx > length(p)) ? 0.0 : p[idx]
                 end
             end
    [
        pv("ilabcosts");
        pv("epsY_1");
        pv("epsM_1");
        log(pv("kappaV"));
        pv("rho_om1");
        pv("sigma_om");
        pv("rho_tfp1");
        [pv("isigma_tfp_$(i)") for i in 1:12];
        pv("rho_pvstar");
        pv("sigma_pvstar");
        pv("rho_xi");
        pv("sigma_xi");
    ]
end


# =========================================================================== #
#  5.  SMM OBJECTIVE FUNCTION                                                  #
# =========================================================================== #

"""
    smm_objective(θ, data_moments, W, context, baseline, endo_names)
      -> (obj::Float64, moments::Vector{Float64})

Compute the SMM loss (θ - data_moments)' W (θ - data_moments).
Returns 1e8 if the model fails to solve.
"""
function smm_objective(θ::AbstractVector{<:Real},
                        data_moments::Vector{<:Real},
                        W::Matrix{<:Real},
                        context::Dynare.Context,
                        baseline::NamedTuple,
                        endo_names::Vector{String})
    moments, ok = smm_model_moments(θ, context, baseline, endo_names)
    if !ok || any(isnan, moments)
        return 1e8, fill(NaN, 46)
    end
    ψ   = data_moments .- moments
    obj = dot(ψ, W * ψ)
    return obj, moments
end


# =========================================================================== #
#  6.  WARM START                                                              #
# =========================================================================== #

function load_warm_start(n_theta::Int)
    # Warm-start from CSV checkpoint (written by previous SMM runs)
    ckpt_file = joinpath(DATA_DIR, "smm_checkpoint.csv")

    if isfile(ckpt_file)
        try
            df = CSV.read(ckpt_file, DataFrame)
            θ_prev = Float64.(df.value)
            if length(θ_prev) == n_theta && all(θ_prev .>= LB) && all(θ_prev .<= UB)
                obj_prev = parse(Float64, String(df[1, :obj_hat]))  # stored in first row
                @printf "  Warm start: smm_checkpoint.csv  (obj=%.6f)\n\n" obj_prev
                return θ_prev
            end
        catch
        end
    end

    @printf "  No valid checkpoint found — using default θ₀.\n\n"
    return nothing
end


# =========================================================================== #
#  7.  MAIN ESTIMATION FUNCTION                                                #
# =========================================================================== #

"""
    smm_run(context) -> (theta_hat, obj_hat, moments_hat)

Run the full SMM estimation.  Pass the Dynare.jl `context` returned by
`@dynare "NK_SOE_lev_gap2"` (from main_SOE_gap.jl or equivalent).

Example:
    include("main_SOE_gap.jl")          # solves model, returns context
    include("smm_estimation.jl")
    theta_hat, obj, moments = smm_run(context)
"""
function smm_run(context::Dynare.Context; endo_names_override=nothing)

    # context.symboltable may have inconsistent variable count (Dynare.jl issue).
    # Use the override (from dynare_endo_names.csv) when provided.
    endo_names = if endo_names_override !== nothing && length(endo_names_override) >= 400
        endo_names_override
    else
        Dynare.get_endogenous(context.symboltable)
    end
    @printf "  endo_names: %d variables\n" length(endo_names)

    # --- Build baseline calibration struct ---
    baseline = build_baseline(context, endo_names,
                               y_d, p_d, l_d, d_TBGDP, d_omG)

    # --- Weighting matrix ---
    W = build_weighting_matrix(data_moments)
    @printf "  Weighting matrix built (%dx%d diagonal).\n\n" size(W,1) size(W,2)

    # --- Initial θ ---
    # default_theta0 reads from _SMM_PARAMS (params_jl.mod).  The baseline
    # has isigma_tfp_i = 0 for all sectors (no TFP shocks in calibration),
    # so the initial model moments are near-zero — very flat landscape.
    # Override the zero shock parameters with data-informed starting values:
    #   isigma_tfp_i ≈ std(Y_i) * 0.4  (rough approximation)
    #   sigma_om     ≈ 0.05
    θ0 = default_theta0(context)
    # Override zero/tiny shock starting values with data-informed guesses
    θ0[6]  = max(θ0[6],  0.05)                          # sigma_om
    for i in 1:12
        θ0[7+i] = max(θ0[7+i], clamp(y_d[i] * 0.4, LB[7+i], UB[7+i]))
    end                                                   # isigma_tfp_i
    θ0[21] = max(θ0[21], 0.05)                           # sigma_pvstar
    θ0[23] = max(θ0[23], 0.05)                           # sigma_xi
    # Try warm start
    θ_warm = load_warm_start(N_THETA)
    θ0     = something(θ_warm, θ0)
    θ0     = clamp.(θ0, LB, UB)

    # --- Diagnostics: check param_idx and resolve_first_order! work --------
    @printf "=== PRE-FLIGHT DIAGNOSTICS ===\n"

    # 1. Check param_idx for a known parameter
    test_params = ["ilabcosts", "kappaV", "rho_om1", "sigma_om",
                   "rho_tfp1", "epsY_1", "epsM_1"]
    @printf "  param_idx results:\n"
    all_idx_ok = true
    for nm in test_params
        idx = param_idx(context, nm)
        @printf "    %-20s → %s\n" nm (idx === nothing ? "MISSING — param update will skip" : string(idx))
        idx === nothing && (all_idx_ok = false)
    end
    all_idx_ok || @printf "  WARNING: some params not found → set_param! will skip them.\n"
    @printf "\n"

    # 2. Check resolve_first_order! works
    @printf "  resolve_first_order! test:\n"
    ok_resolve, g_test, _, _ = resolve_first_order!(context)
    @printf "    success=%s  g1_1 size=%s\n\n" ok_resolve (ok_resolve ? string(size(g_test)) : "N/A")

    if !ok_resolve
        @printf "  PRE-FLIGHT: model re-solve failed.\n"
        @printf "  The Klein (2000) pure-Julia solver uses:\n"
        @printf "    - SparseDynamicG1!.jl (compiled Jacobian)\n"
        @printf "    - GenericSchur.jl (QZ decomposition, no LAPACK)\n"
        @printf "  Check that mod/NK_SOE_lev_gap2/model/julia/ exists\n"
        @printf "  and that main_SOE_gap.jl ran successfully first.\n\n"
        error("Model re-solve failed at pre-flight. Run main_SOE_gap.jl first.")
    end

    # 3. Full pre-flight
    @printf "=== PRE-FLIGHT CHECK ===\n"
    m_test, ok_test = smm_model_moments(θ0, context, baseline, endo_names)
    if !ok_test || any(isnan, m_test)
        @printf "  PRE-FLIGHT FAILED — model did not solve at θ₀.\n"
        @printf "  NaN moment indices: %s\n" string(findall(isnan, m_test))
        @printf "  → Check that main_SOE_gap.jl ran EXERCISE=0 (Baseline) first.\n"
        error("SMM pre-flight failed. Fix model setup before estimation.")
    end
    obj_test = smm_objective(θ0, data_moments, W, context, baseline, endo_names)[1]

    @printf "  %-34s  %9s  %9s\n" "Moment" "Data" "Model"
    @printf "  %s\n" repeat("-", 56)
    for (i, nm) in enumerate(MOMENT_NAMES)
        @printf "  %-34s  %9.5f  %9.5f\n" nm data_moments[i] m_test[i]
    end
    @printf "\n  obj(θ₀) = %.6f\n" obj_test
    @printf "=== PRE-FLIGHT PASSED — launching optimizer ===\n\n"

    # --- CMA-ES optimisation ---
    # CMA-ES (Covariance Matrix Adaptation Evolution Strategy) is a state-of-the-art
    # black-box derivative-free optimizer — the standard choice for SMM estimation.
    max_evals = 30_000   # ~1300 generations for n=23; enough for convergence
    @printf "--- CMA-ES (black-box, derivative-free) ---\n"
    @printf "  %d parameters  |  %d moments  |  max %d evaluations\n" N_THETA 46 max_evals
    @printf "  %-6s  %-14s  %-10s\n" "eval" "objective" "BK/NaN"
    @printf "  %s\n" repeat("-", 36)

    # Track best solution manually — robust across CMAEvolutionStrategy API versions.
    best_θ   = Ref(clamp.(θ0, LB, UB))
    best_obj = Ref(Inf)
    fail_count = Ref(0)
    eval_count = Ref(0)

    obj_fn = θ -> begin
        obj, _ = smm_objective(θ, data_moments, W, context, baseline, endo_names)
        eval_count[] += 1
        if isfinite(obj) && obj < best_obj[]
            best_obj[] = obj
            best_θ[]   = copy(θ)
        end
        obj >= 1e7 && (fail_count[] += 1)
        # Print progress every 200 evaluations
        if eval_count[] % 200 == 0
            @printf "  %-6d  %-14.6f  %-10d\n" eval_count[] best_obj[] fail_count[]
            flush(stdout)
        end
        obj
    end

    # Scalar σ0: mean of (UB-LB)/6 spans ~1/6 of the feasible range.
    insigma = mean((UB .- LB) ./ 6)

    CMAEvolutionStrategy.minimize(
        obj_fn,
        clamp.(θ0, LB, UB),
        insigma;
        lower     = LB,
        upper     = UB,
        maxiter   = max_evals,
        ftol      = 1e-6,
        xtol      = 1e-6,
        seed      = 42,
        verbosity = 0,   # suppress CMA-ES internal output; we print our own
    )

    θ_hat = clamp.(best_θ[], LB, UB)
    @printf "  %s\n" repeat("-", 36)
    @printf "  %-6d  %-14.6f  %-10d\n\n" eval_count[] best_obj[] fail_count[]
    @printf "CMA-ES done after %d evaluations.  Best obj = %.6f\n\n" eval_count[] best_obj[]

    # --- Final evaluation at best θ ---
    obj_hat, moments_hat = smm_objective(θ_hat, data_moments, W, context, baseline, endo_names)
    ψ_hat = data_moments .- moments_hat

    # =========================================================================== #
    #  8.  RESULTS TABLE                                                           #
    # =========================================================================== #

    @printf "\n%s\n  SMM RESULTS\n%s\n\n" repeat("=",60) repeat("=",60)
    @printf "  Objective at θ̂: %.6f\n\n" obj_hat

    @printf "  %-6s  %-18s  %10s  %10s\n" "Idx" "Parameter" "Initial" "Estimate"
    @printf "  %s\n" repeat("-", 50)
    for k in 1:N_THETA
        if k == 4
            @printf "  %3d  %-18s  %10.4f  %10.4f  [log]\n" k PARAM_LABELS[k] θ0[k] θ_hat[k]
            @printf "  %3s  %-18s  %10.2e  %10.2e  [level]\n" "--" "kappaV" exp(θ0[k]) exp(θ_hat[k])
        else
            @printf "  %3d  %-18s  %10.4f  %10.4f\n" k PARAM_LABELS[k] θ0[k] θ_hat[k]
        end
    end

    @printf "\n  Moment fit:\n"
    @printf "  %-36s  %9s  %9s  %9s\n" "Moment" "Data" "Model" "Diff"
    @printf "  %s\n" repeat("-", 68)
    for (i, nm) in enumerate(MOMENT_NAMES)
        @printf "  %-36s  %9.5f  %9.5f  %+9.5f\n" nm data_moments[i] moments_hat[i] ψ_hat[i]
    end

    # =========================================================================== #
    #  9.  SAVE RESULTS                                                            #
    # =========================================================================== #

    # smm_results.mat — full results (compatible with MATLAB smm_estimation.m format)
    smm_results_dict = Dict{String,Any}(
        "theta_hat"     => θ_hat,
        "obj_hat"       => obj_hat,
        "psi_hat"       => ψ_hat,
        "moments_hat"   => moments_hat,
        "data_moments"  => data_moments,
        "param_labels"  => PARAM_LABELS,
        "moment_names"  => MOMENT_NAMES,
        "W"             => W,
    )
    # ---- smm_results.csv — full results table --------------------------------
    df_results = DataFrame(
        param   = vcat(PARAM_LABELS, fill("", 46 - N_THETA)),
        theta   = vcat(θ_hat,        fill(NaN, 46 - N_THETA)),
        moment  = MOMENT_NAMES,
        data    = data_moments,
        model   = moments_hat,
        diff    = ψ_hat,
    )
    results_path = joinpath(DATA_DIR, "smm_results.csv")
    CSV.write(results_path, df_results)
    @printf "\nResults saved to:\n  %s\n" results_path

    # ---- smm_estimates.csv — parameter estimates (read by main_SOE_gap.jl) --
    sector_labels = ["isigma_tfp_$(i)" for i in 1:NSEC]
    df_est = DataFrame(
        param = vcat(["ilabcosts","epsY","epsM","log_kappaV",
                      "rho_om","sigma_om","rho_A"],
                     sector_labels,
                     ["rho_pvstar","sigma_pvstar","rho_xi","sigma_xi"]),
        value = θ_hat,
        obj_hat = vcat([obj_hat], fill(NaN, N_THETA-1)),
    )
    estimates_path = joinpath(DATA_DIR, "smm_estimates.csv")
    CSV.write(estimates_path, df_est)
    @printf "  Estimates saved to:\n  %s\n" estimates_path

    # ---- checkpoint for warm starts ----------------------------------------
    CSV.write(joinpath(DATA_DIR, "smm_checkpoint.csv"), df_est)

    @printf "\n  Re-run main_SOE_gap.jl (EXERCISE=0) to apply estimates.\n\n"

    return θ_hat, obj_hat, moments_hat
end


# =========================================================================== #
#  ENTRY POINT                                                                 #
#                                                                              #
#  When run as a script, this tries to load the Dynare context from the       #
#  compiled model output (written by main_SOE_gap.jl).  If not found, it      #
#  tells the user to run main_SOE_gap.jl first.                               #
# =========================================================================== #

if abspath(PROGRAM_FILE) == @__FILE__
    # Primary: context saved by run_dynare_subprocess.jl
    context_file = joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls")
    # Fallback: Dynare.jl's own cache location
    context_file_alt = joinpath(SCRIPT_DIR, "mod", "NK_SOE_lev_gap2",
                                "output", "NK_SOE_lev_gap2.jls")
    ctx_path = isfile(context_file) ? context_file :
               isfile(context_file_alt) ? context_file_alt : ""

    if !isempty(ctx_path)
        using Serialization
        @printf "Loading compiled Dynare context from:\n  %s\n\n" ctx_path
        context = deserialize(ctx_path)
        smm_run(context)
    else
        @printf """
        Dynare context not found at:
          %s

        Run main_SOE_gap.jl first (EXERCISE=0) to compile the model:
          julia --project=. main_SOE_gap.jl

        Then re-run SMM estimation:
          julia --project=. smm_estimation.jl

        Or from a Julia REPL:
          include("main_SOE_gap.jl")   # returns `context`
          include("smm_estimation.jl")
          smm_run(context)
        """ context_file
    end
end
