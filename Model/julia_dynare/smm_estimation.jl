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

using LinearAlgebra, Statistics, StatsBase, Printf
using NLsolve, CSV, DataFrames, Dynare
using CMAEvolutionStrategy

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
    # Fallback to 0.025 (≈ Chilean historical estimate) so the optimiser never sees NaN.
    d_std_TBGDP  = get(agg_dict, "std_TBGDP", 0.025)
    @printf "  Loaded sectoral_moments.csv + aggregate_moments.csv\n\n"
else
    @printf "  CSV not found in %s — using placeholder values.\n\n" DATA_DIR
    y_d = fill(0.04, NSEC); p_d = fill(0.02, NSEC); l_d = fill(0.03, NSEC)
    d_std_GDP=0.0215; d_std_pi=0.0041; d_corr_GDPpi=-0.15
    d_omG=0.57; d_std_Q=0.0520; d_autocorr_Q=0.75
    d_corr_GDPQ=-0.15; d_TBGDP=-0.02; d_std_TBGDP=0.025
end

# omG (goods expenditure share) is pinned by SS calibration → contributes 0 to the
# loss and carries no information for estimation.  Replaced by std(TB/GDP), which
# directly identifies the export and import elasticities.
data_moments = [y_d; p_d; l_d; d_std_GDP; d_std_pi; d_corr_GDPpi;
                d_std_TBGDP; d_std_Q; d_autocorr_Q; d_corr_GDPQ; 1.0; 1.0; 1.0]
@assert length(data_moments) == 46 "Expected 46 data moments, got $(length(data_moments))"

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

const MOMENT_NAMES = vcat(
    ["std(Y_$(i))"  for i in 1:12], ["std(PH_$(i))" for i in 1:12],
    ["std(L_$(i))"  for i in 1:12],
    ["std(GDP)", "std(pi)", "corr(GDP,pi)", "std(TB/GDP)",
     "std(Q)", "autocorr(Q)", "corr(GDP,Q)",
     "rank corr: output (model vs data)", "rank corr: prices (model vs data)",
     "rank corr: labor  (model vs data)"])


# =========================================================================== #
#  2.  PARAMETER BOUNDS (tighter, economically motivated)                     #
# =========================================================================== #

const PARAM_LABELS = vcat(
    ["ilabcosts", "epsY", "epsM", "log(kappaV)", "rho_om", "sigma_om", "rho_A"],
    ["isigma_tfp_$(i)" for i in 1:12],
    ["rho_pvstar", "sigma_pvstar", "rho_xi", "sigma_xi"],
    ["etastar"],   # export demand elasticity η*; identifies std(TB/GDP)
)
const N_THETA = length(PARAM_LABELS)   # 24

const LB = [1e-3; 0.30; 0.05; log(1e3);  -0.95; 1e-5;  0.10;
            fill(1e-4, 12); 0.50;  0.005; 0.00; 0.0;  0.50]
const UB = [50.0; 1.50; 0.50; log(1e8);   0.95; 0.20;  0.95;
            fill(0.10, 12); 0.99;  0.20;  0.95; 0.05; 6.00]


# =========================================================================== #
#  3.  WEIGHTING MATRIX (proportional + high rank-corr weight)                #
# =========================================================================== #

function build_weighting_matrix(dm::Vector{<:Real})
    w = ones(46)
    # Standard inverse-variance weight for: sectoral stds (1-36), std(GDP) (37),
    # std(pi) (38), std(TB/GDP) (40), std(Q) (41).
    # Indices 39,42,43,44,45,46 start at w=1.0 and are overridden below.
    for k in vcat(1:36, [37, 38, 40, 41])
        d = abs(dm[k]); w[k] = d > 1e-4 ? 1.0/d^2 : 1.0/0.01^2
    end
    w[39] *= 0.20     # corr(GDP,π): structurally hard for supply-shock model
    # w[40]: std(TB/GDP) — now gets standard inverse-variance weight (no override)
    w[42]  = 2.0      # autocorr(Q): identifies rho_pvstar
    w[44]  = 2.0; w[45] = 2.0; w[46] = 2.0   # rank correlations: moderate (reduced from 5× to avoid dominating loss)
    return Diagonal(w) |> Matrix
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
    [pv("ilabcosts"); pv("epsY_1"); pv("epsM_1"); log(pv("kappaV"));
     pv("rho_om1"); pv("sigma_om"); pv("rho_tfp1");
     [pv("isigma_tfp_$(i)") for i in 1:12];
     pv("rho_pvstar"); pv("sigma_pvstar"); pv("rho_xi"); pv("sigma_xi");
     pv("etastar")]
end


# =========================================================================== #
#  6.  WARM START                                                              #
# =========================================================================== #

function load_warm_start(n_theta::Int)
    ckpt = joinpath(DATA_DIR, "smm_checkpoint.csv")
    if isfile(ckpt)
        try
            df = CSV.read(ckpt, DataFrame)
            θ_prev = Float64.(df.value)
            if length(θ_prev) == n_theta
                θ_clamped = clamp.(θ_prev, LB, UB)
                obj_prev = try parse(Float64, String(df[1, :obj_hat])) catch; NaN; end
                @printf "  Warm start: smm_checkpoint.csv (obj=%.6f)\n\n" obj_prev
                return θ_clamped
            end
        catch; end
    end
    @printf "  No valid checkpoint — using default θ₀.\n\n"
    return nothing
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

function smm_model_moments(θ::AbstractVector{<:Real}, context, baseline, endo_names)
    nsec = baseline.nsec; NAN46 = fill(NaN, 46)
    ilabcosts=θ[1]; epsY=θ[2]; epsM=θ[3]; kappaV=exp(θ[4])
    rho_om=θ[5]; sigma_om=θ[6]; rho_A=θ[7]
    isigma_tfp = @view θ[8:19]
    rho_pvstar=θ[20]; sigma_pvstar=θ[21]; rho_xi=θ[22]; sigma_xi=θ[23]
    etastar = length(θ) >= 24 ? θ[24] : baseline.etastar_val   # export demand elasticity η*

    (!(0<epsY<5)||!(0<epsM<2)||ilabcosts<=0||kappaV<=0||abs(rho_om)>=1||
     sigma_om<0||abs(rho_A)>=1||any(<(0),isigma_tfp)||abs(rho_pvstar)>=1||
     sigma_pvstar<0||abs(rho_xi)>=1||sigma_xi<0||
     !(0.1<etastar<8.0)) && return NAN46, false

    set_param!(context,"ilabcosts",ilabcosts); set_param!(context,"kappaV",kappaV)
    set_param!(context,"rho_om1",rho_om);      set_param!(context,"sigma_om",sigma_om)
    set_param!(context,"etastar",etastar)
    set_param!(context,"rho_tfp1",rho_A);      set_param!(context,"rho_pvstar",rho_pvstar)
    set_param!(context,"sigma_pvstar",sigma_pvstar); set_param!(context,"rho_xi",rho_xi)
    set_param!(context,"sigma_xi",sigma_xi)
    for i in 1:nsec; set_param!(context,"isigma_tfp_$(i)",isigma_tfp[i]); end

    epsY_prev    = get_param_val(context,"epsY_1"); epsM_prev = get_param_val(context,"epsM_1")
    etastar_prev = get_param_val(context,"etastar")
    need_ss = abs(epsY-epsY_prev)>1e-8 || abs(epsM-epsM_prev)>1e-8 ||
              abs(etastar - (isnan(etastar_prev) ? baseline.etastar_val : etastar_prev)) > 1e-8
    for i in 1:nsec; set_param!(context,"epsY_$(i)",epsY); set_param!(context,"epsM_$(i)",epsM); end
    if need_ss
        ok = recompute_ss_cached!(context, epsY, epsM, baseline, endo_names; etastar=etastar)
        !ok && return NAN46, false
    end

    success, T, R = _resolve_cached!(context, θ)
    !success && return NAN46, false

    n_exo   = size(R, 2)
    sr      = context.models[1].i_bkwrd_b
    n_state = length(sr)
    sc      = _get_scratch!(n_exo, n_state, endo_names)

    # Σe: activate all shocks except epschi (unused labor supply shock)
    fill!(sc.Σe, 0.0)
    sc.Σe[1,1]=1.0; sc.Σe[2,2]=1.0; sc.Σe[4,4]=1.0   # eps_om, eps_i, eps_pvstar
    for i in 5:min(16,n_exo); sc.Σe[i,i]=1.0; end     # epsA_1:12
    n_exo>=17 && (sc.Σe[17,17]=1.0)                    # eps_xi

    Tsr = T[sr,:]; Rsr = R[sr,:]
    RΣ  = Rsr * sc.Σe   # 78×17
    mul!(sc.Q_lyap, RΣ, Rsr')
    @inbounds for j in 1:n_state, i in 1:j-1
        avg=0.5*(sc.Q_lyap[i,j]+sc.Q_lyap[j,i]); sc.Q_lyap[i,j]=avg; sc.Q_lyap[j,i]=avg
    end

    dlyap_inplace!(sc.P, Tsr, sc.Q_lyap; tmp1=sc.tmp1, tmp2=sc.tmp2, Ak=sc.Ak)
    (any(diag(sc.P).<-1e-10) || any(isnan.(sc.P))) && return NAN46, false

    # HP-filtered covariance on the ~40 needed variables only (25× faster)
    needed_names = vcat(["Y_$(i)" for i in 1:nsec], ["PH_$(i)" for i in 1:nsec],
                        ["L_$(i)" for i in 1:nsec], ["GDP","pi","Q","TB"])
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

    @inline pstd(vn) = let k=get(ei_sub,vn,0)
        k==0 ? 0.0 : sqrt(max(Γ_v[k,k],0.0)) / max(abs(ys[needed_idx[k]]),1e-12)
    end
    @inline xcorr(v1,v2) = let i1=get(ei_sub,v1,0), i2=get(ei_sub,v2,0)
        (i1==0||i2==0) ? NaN :
        let d=sqrt(max(Γ_v[i1,i1],0.0)*max(Γ_v[i2,i2],0.0))
            d<1e-15 ? 0.0 : clamp(Γ_v[i1,i2]/d,-1.0,1.0)
        end
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

    dY=baseline.data_std_Y; dPH=baseline.data_std_PH; dL=baseline.data_std_L
    vy=isfinite.(std_Y).&isfinite.(dY)
    vp=isfinite.(std_PH).&isfinite.(dPH)
    vl=isfinite.(std_L).&isfinite.(dL)
    rY = sum(vy)>=3 ? safe_spearman(std_Y[vy],dY[vy]) : 0.0
    rP = sum(vp)>=3 ? safe_spearman(std_PH[vp],dPH[vp]) : 0.0
    rL = sum(vl)>=3 ? safe_spearman(std_L[vl],dL[vl]) : 0.0

    return [std_Y;std_PH;std_L;pstd("GDP");pstd("pi");xcorr("GDP","pi");
            std_TBGDP;pstd("Q");acQ;xcorr("GDP","Q");rY;rP;rL], true
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

    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)

    n_threads_active = Threads.nthreads()
    contexts_th = n_threads_active > 1 ?
        ((@printf "  Spawning %d per-thread contexts...\n" n_threads_active);
         [deepcopy(context) for _ in 1:n_threads_active]) : [context]
    baselines_th = [deepcopy(baseline) for _ in 1:length(contexts_th)]

    W  = build_weighting_matrix(data_moments)
    @printf "  Weighting: proportional std devs + 5× rank correlations.\n\n"

    # Initial θ: data-informed starting values for shock parameters
    θ0 = default_theta0(context)
    θ0[6]  = max(θ0[6],  0.05)
    for i in 1:12; θ0[7+i] = max(θ0[7+i], clamp(y_d[i]*0.4, LB[7+i], UB[7+i])); end
    θ0[21] = max(θ0[21], 0.05)
    θ0[23] = max(θ0[23], 0.02)
    θ_warm = load_warm_start(N_THETA)
    θ0 = clamp.(something(θ_warm, θ0), LB, UB)

    # Pre-flight
    @printf "=== PRE-FLIGHT ===\n"
    ok_resolve, g_test, _, _ = resolve_first_order!(context)
    if !ok_resolve
        error("Model re-solve failed at pre-flight. Run main_SOE_gap.jl first.")
    end
    @printf "  resolve_first_order! OK, g1_1 size=%s\n" string(size(g_test))

    m_test, ok_test = smm_model_moments(θ0, context, baseline, endo_names)
    if !ok_test || any(isnan, m_test)
        @printf "  FAILED: NaN at %s\n" string(findall(isnan, m_test))
        error("Pre-flight failed. Fix model setup before estimation.")
    end
    obj_test    = dot(data_moments .- m_test, W * (data_moments .- m_test))
    m_test_copy = copy(m_test)   # used as fallback for best_moments below
    @printf "  obj(θ₀) = %.6f\n" obj_test

    ψ0 = data_moments .- m_test
    @printf "  Decomp: Y=%.3f PH=%.3f L=%.3f Agg=%.3f Rank=%.3f\n" dot(ψ0[1:12],W[1:12,1:12]*ψ0[1:12]) dot(ψ0[13:24],W[13:24,13:24]*ψ0[13:24]) dot(ψ0[25:36],W[25:36,25:36]*ψ0[25:36]) dot(ψ0[37:43],W[37:43,37:43]*ψ0[37:43]) dot(ψ0[44:46],W[44:46,44:46]*ψ0[44:46])

    @printf "\n  %-34s  %9s  %9s\n" "Moment" "Data" "Model"
    @printf "  %s\n" repeat("-",56)
    for (i,nm) in enumerate(MOMENT_NAMES)
        @printf "  %-34s  %9.5f  %9.5f\n" nm data_moments[i] m_test[i]
    end
    @printf "=== PRE-FLIGHT PASSED ===\n\n"

    # CMA-ES
    max_evals = 300_000
    @printf "--- CMA-ES ---\n"
    @printf "  %d params | %d moments | max %d evals | %d threads\n" N_THETA 46 max_evals n_threads_active
    @printf "  %-6s  %-10s  %-44s  %-7s  %-8s  %-8s\n" "eval" "best_obj" "[Y    PH    L     Agg   Rank ]" "fail" "ms/eval" "Klein%"
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

    _tid_cma() = min(Threads.threadid(), length(contexts_th))

    obj_fn = θ_sc -> begin
        tid   = _tid_cma()
        θ     = LB .+ θ_sc .* span    # unscale [0,1] → original parameter space
        m, ok = smm_model_moments(θ, contexts_th[tid], baselines_th[tid], endo_names)
        obj   = (!ok || any(isnan, m)) ? 1e8 : dot(data_moments.-m, W*(data_moments.-m))

        n = Threads.atomic_add!(eval_count, 1) + 1   # new count (1-based)
        obj >= 1e7 && Threads.atomic_add!(fail_count, 1)

        if isfinite(obj) && obj < best_obj[]
            lock(best_lock) do
                if isfinite(obj) && obj < best_obj[]
                    best_obj[]      = obj
                    best_θ[]        = copy(θ)   # θ in original parameter space
                    best_moments[]  = copy(m)   # save moments — avoids re-evaluation bug
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
                        dY  = dot(ψ_now[1:12],  W[1:12,1:12]   * ψ_now[1:12])
                        dPH = dot(ψ_now[13:24], W[13:24,13:24] * ψ_now[13:24])
                        dL  = dot(ψ_now[25:36], W[25:36,25:36] * ψ_now[25:36])
                        dAg = dot(ψ_now[37:43], W[37:43,37:43] * ψ_now[37:43])
                        dRk = dot(ψ_now[44:46], W[44:46,44:46] * ψ_now[44:46])
                        @printf "  %-6d  %-10.4f  [Y=%.2f PH=%.2f L=%.2f Agg=%.2f Rk=%.2f]  fail=%-5d  %.1fms  Klein=%d%%\n" n b_obj dY dPH dL dAg dRk fail_count[] ms kpct
                        flush(stdout)
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
    insigma = 0.08   # 8% of [0,1] after rescaling — tighter start prevents premature step-size collapse

    # Rescale θ to [0,1] so CMA-ES works in a unit hypercube.
    # The objective wrapper maps back to the original scale.
    span  = UB .- LB
    θ0_sc = (clamp.(θ0, LB, UB) .- LB) ./ span     # scaled θ₀ ∈ [0,1]
    LB_sc = zeros(N_THETA)
    UB_sc = ones(N_THETA)

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

    θ_hat = clamp.(best_θ[], LB, UB)
    total_t = time() - t_start[]
    tot_kl  = _KLEIN_HITS[] + _KLEIN_MISSES[]
    @printf "\nDone: %d evals | %.1fs | %.1fms/eval | Klein cache=%d%%\n\n" eval_count[] total_t 1000*total_t/max(1,eval_count[]) round(Int,100*_KLEIN_HITS[]/max(1,tot_kl))

    # Use moments cached at best θ — re-evaluating on the main context would give a
    # different result because context has cold SS state (frozen at θ₀ from pre-flight),
    # while the best was found on a per-thread context with warm SS cache.
    obj_hat = best_obj[]
    m_hat   = best_moments[]
    ψ_hat   = data_moments .- m_hat

    # Results table
    @printf "\n%s\n  SMM RESULTS\n%s\n\n" repeat("=",60) repeat("=",60)
    @printf "  Objective: %.6f\n\n" obj_hat
    @printf "  %-6s  %-18s  %10s  %10s\n" "Idx" "Parameter" "Initial" "Estimate"
    @printf "  %s\n" repeat("-",50)
    for k in 1:N_THETA
        if k==4
            @printf "  %3d  %-18s  %10.4f  %10.4f  [log]\n" k PARAM_LABELS[k] θ0[k] θ_hat[k]
            @printf "  %3s  %-18s  %10.2e  %10.2e  [level]\n" "--" "kappaV" exp(θ0[k]) exp(θ_hat[k])
        else
            @printf "  %3d  %-18s  %10.4f  %10.4f\n" k PARAM_LABELS[k] θ0[k] θ_hat[k]
        end
    end
    @printf "\n  %-36s  %9s  %9s  %9s  %6s\n" "Moment" "Data" "Model" "Diff" "Rel%"
    @printf "  %s\n" repeat("-",75)
    for (i,nm) in enumerate(MOMENT_NAMES)
        re = abs(data_moments[i])>1e-6 ? abs(ψ_hat[i])/abs(data_moments[i])*100 : 0.0
        @printf "  %-36s  %9.5f  %9.5f  %+9.5f  %5.1f%%\n" nm data_moments[i] m_hat[i] ψ_hat[i] re
    end
    @printf "\n  Decomp: Y=%.3f PH=%.3f L=%.3f Agg=%.3f Rank=%.3f\n\n" dot(ψ_hat[1:12],W[1:12,1:12]*ψ_hat[1:12]) dot(ψ_hat[13:24],W[13:24,13:24]*ψ_hat[13:24]) dot(ψ_hat[25:36],W[25:36,25:36]*ψ_hat[25:36]) dot(ψ_hat[37:43],W[37:43,37:43]*ψ_hat[37:43]) dot(ψ_hat[44:46],W[44:46,44:46]*ψ_hat[44:46])

    # Save
    df_res = DataFrame(param=vcat(PARAM_LABELS,fill("",46-N_THETA)),
                        theta=vcat(θ_hat,fill(NaN,46-N_THETA)),
                        moment=MOMENT_NAMES, data=data_moments, model=m_hat, diff=ψ_hat)
    CSV.write(joinpath(DATA_DIR,"smm_results.csv"), df_res)

    df_est = DataFrame(
        param=vcat(["ilabcosts","epsY","epsM","log_kappaV","rho_om","sigma_om","rho_A"],
                   ["isigma_tfp_$(i)" for i in 1:NSEC],
                   ["rho_pvstar","sigma_pvstar","rho_xi","sigma_xi","etastar"]),
        value=θ_hat, obj_hat=vcat([obj_hat],fill(NaN,N_THETA-1)))
    CSV.write(joinpath(DATA_DIR,"smm_estimates.csv"), df_est)
    CSV.write(joinpath(DATA_DIR,"smm_checkpoint.csv"), df_est)

    @printf "  Results: %s\n  Estimates: %s\n\n" joinpath(DATA_DIR,"smm_results.csv") joinpath(DATA_DIR,"smm_estimates.csv")
    @printf "  Re-run main_SOE_gap.jl (EXERCISE=0) to apply estimates.\n\n"

    return θ_hat, obj_hat, m_hat
end
