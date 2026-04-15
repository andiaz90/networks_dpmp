"""
smm_estimation_v2.jl
====================
Optimized SMM estimation for the NK-IOSOE Chile model.

KEY PERFORMANCE CHANGES:
  1. SS caching with warm-start tolerance (avoids redundant nlsolve)
  2. Pre-allocated scratch arrays for Σe, moment vector, string→index dict
  3. In-place Lyapunov solver (no temporaries in hot loop)
  4. @views for θ slicing (avoids allocation)
  5. Pre-computed endo_name → index dictionary (no string lookup per moment)
  6. Tighter bounds, proportional weighting, higher rank-corr weights

PREREQUISITES
  Same as smm_estimation.jl — run main_SOE_gap.jl first.
"""

using LinearAlgebra, Statistics, StatsBase, Printf
using NLsolve, CSV, DataFrames, Dynare
using CMAEvolutionStrategy

SCRIPT_DIR  = @__DIR__
REPO_ROOT   = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR    = joinpath(REPO_ROOT, "Data")

@printf "\n%s\n  SMM ESTIMATION v2: NK-IOSOE Chile Model\n%s\n\n" repeat("=",60) repeat("=",60)


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
    @printf "  CSV not found — using placeholder values.\n\n"
    y_d = fill(0.04, NSEC); p_d = fill(0.02, NSEC); l_d = fill(0.03, NSEC)
    d_std_GDP=0.0215; d_std_pi=0.0041; d_corr_GDPpi=-0.15
    d_omG=0.57; d_std_Q=0.0520; d_autocorr_Q=0.75
    d_corr_GDPQ=-0.15; d_TBGDP=-0.02
end

data_moments = [y_d; p_d; l_d; d_std_GDP; d_std_pi; d_corr_GDPpi;
                d_omG; d_std_Q; d_autocorr_Q; d_corr_GDPQ; 1.0; 1.0; 1.0]
@assert length(data_moments) == 46


# =========================================================================== #
#  2.  TIGHTER PARAMETER BOUNDS                                                #
# =========================================================================== #

const PARAM_LABELS_V2 = vcat(
    ["ilabcosts", "epsY", "epsM", "log(kappaV)", "rho_om", "sigma_om", "rho_A"],
    ["isigma_tfp_$(i)" for i in 1:12],
    ["rho_pvstar", "sigma_pvstar", "rho_xi", "sigma_xi"],
)
const N_THETA_V2 = length(PARAM_LABELS_V2)

const LB_V2 = [1e-3; 0.30; 0.05; log(1e3);  -0.95; 1e-5;  0.10;
               fill(1e-4, 12); 0.50;  0.005; 0.00;  0.0  ]
const UB_V2 = [50.0; 1.50; 0.50; log(1e8);   0.95; 0.20;  0.95;
               fill(0.10, 12); 0.99;  0.20;  0.95;  0.05 ]


# =========================================================================== #
#  3.  PROPORTIONAL WEIGHTING MATRIX                                           #
# =========================================================================== #

function build_weighting_matrix_v2(dm::Vector{<:Real})
    w = ones(46)
    std_idx = vcat(1:36, [37, 38, 41])
    for k in std_idx
        d = abs(dm[k])
        w[k] = d > 1e-4 ? 1.0/d^2 : 1.0/0.01^2
    end
    w[39] = 1.0 * 0.20     # corr(GDP,π): structurally hard
    w[43] = 1.0             # corr(GDP,Q)
    w[42] = 2.0             # autocorr(Q): identifies rho_pvstar
    w[40] = 0.1             # omG: SS target
    w[44] = 5.0; w[45] = 5.0; w[46] = 5.0   # rank correlations: HIGH
    return Diagonal(w) |> Matrix
end


# =========================================================================== #
#  4.  SS CACHE                                                                #
# =========================================================================== #

mutable struct SSCache
    epsY::Float64
    epsM::Float64
    ss_vec::Vector{Float64}
    valid::Bool
end
SSCache() = SSCache(NaN, NaN, Float64[], false)

const _SS_CACHES_V2 = [SSCache() for _ in 1:max(1, Threads.nthreads())]
const _SS_TOL_V2 = 1e-4

function recompute_ss_cached!(context, epsY, epsM, baseline, endo_names)
    tid = min(Threads.threadid(), length(_SS_CACHES_V2))
    cache = _SS_CACHES_V2[tid]

    if cache.valid && abs(epsY - cache.epsY) < _SS_TOL_V2 && abs(epsM - cache.epsM) < _SS_TOL_V2
        for i in 1:baseline.nsec
            set_param!(context, "epsY_$(i)", epsY)
            set_param!(context, "epsM_$(i)", epsM)
        end
        ss_mut = context.results.model_results[1].trends.endogenous_steady_state
        if length(cache.ss_vec) == length(ss_mut)
            copyto!(ss_mut, cache.ss_vec)
        end
        return true
    end

    ok = recompute_ss!(context, epsY, epsM, baseline, endo_names)
    if ok
        ss_vec = context.results.model_results[1].trends.endogenous_steady_state
        cache.epsY = epsY
        cache.epsM = epsM
        cache.ss_vec = copy(ss_vec)
        cache.valid = true
    end
    return ok
end


# =========================================================================== #
#  5.  IN-PLACE LYAPUNOV SOLVER (avoids temporaries in hot loop)              #
# =========================================================================== #

"""
    local_dlyap_inplace!(P, A, Q; tmp1, tmp2, Ak)

Solve P = A P Aᵀ + Q by matrix-doubling, using pre-allocated temporaries.
Overwrites P in-place.
"""
function local_dlyap_inplace!(P::Matrix{Float64}, A::AbstractMatrix{Float64},
                               Q::Matrix{Float64};
                               tmp1::Matrix{Float64}, tmp2::Matrix{Float64},
                               Ak::Matrix{Float64},
                               maxiter::Int=100, tol::Real=1e-12)
    copyto!(P, Q)
    copyto!(Ak, A)
    n = size(P, 1)
    for _ in 1:maxiter
        # tmp1 = Ak * P
        mul!(tmp1, Ak, P)
        # tmp2 = tmp1 * Ak' = Ak * P * Ak'
        mul!(tmp2, tmp1, Ak')
        # P_new = tmp2 + P
        err = 0.0
        @inbounds for j in 1:n, i in 1:n
            new_val = tmp2[i,j] + P[i,j]
            err = max(err, abs(new_val - P[i,j]))
            P[i,j] = new_val
        end
        # Symmetrize
        @inbounds for j in 1:n, i in 1:j-1
            avg = 0.5 * (P[i,j] + P[j,i])
            P[i,j] = avg; P[j,i] = avg
        end
        # Ak = Ak * Ak
        copyto!(tmp1, Ak)
        mul!(Ak, tmp1, tmp1)
        # Check convergence
        pnorm = norm(P, Inf) + 1e-14
        err / pnorm < tol && break
    end
    # Final symmetrize
    @inbounds for j in 1:n, i in 1:j-1
        avg = 0.5 * (P[i,j] + P[j,i])
        P[i,j] = avg; P[j,i] = avg
    end
    return P
end


# =========================================================================== #
#  6.  PER-THREAD SCRATCH ARRAYS                                               #
# =========================================================================== #

mutable struct SMMScratch
    Σe::Matrix{Float64}       # n_exo × n_exo
    Q_lyap::Matrix{Float64}   # n_state × n_state
    P::Matrix{Float64}        # n_state × n_state
    tmp1::Matrix{Float64}     # n_state × n_state
    tmp2::Matrix{Float64}     # n_state × n_state
    Ak::Matrix{Float64}       # n_state × n_state
    moments::Vector{Float64}  # 46
    endo_idx::Dict{String,Int}
    initialized::Bool
end
SMMScratch() = SMMScratch(
    zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0),
    zeros(46), Dict{String,Int}(), false)

const _SCRATCH_V2 = [SMMScratch() for _ in 1:max(1, Threads.nthreads())]

function _get_scratch(n_exo::Int, n_state::Int, endo_names::Vector{String})
    tid = min(Threads.threadid(), length(_SCRATCH_V2))
    s = _SCRATCH_V2[tid]
    if !s.initialized || size(s.Σe, 1) != n_exo || size(s.P, 1) != n_state
        s.Σe     = zeros(n_exo, n_exo)
        s.Q_lyap = zeros(n_state, n_state)
        s.P      = zeros(n_state, n_state)
        s.tmp1   = zeros(n_state, n_state)
        s.tmp2   = zeros(n_state, n_state)
        s.Ak     = zeros(n_state, n_state)
        s.moments = zeros(46)
        s.endo_idx = Dict(nm => i for (i,nm) in enumerate(endo_names))
        s.initialized = true
    end
    return s
end


# =========================================================================== #
#  7.  OPTIMIZED MOMENT FUNCTION                                               #
# =========================================================================== #

function smm_model_moments_v2(θ::AbstractVector{<:Real}, context, baseline, endo_names)
    nsec = baseline.nsec
    NAN46 = fill(NaN, 46)

    # Use @view to avoid allocating a new vector for isigma_tfp
    ilabcosts = θ[1]; epsY = θ[2]; epsM = θ[3]; kappaV = exp(θ[4])
    rho_om = θ[5]; sigma_om = θ[6]; rho_A = θ[7]
    isigma_tfp = @view θ[8:19]   # NO ALLOCATION
    rho_pvstar = θ[20]; sigma_pvstar = θ[21]; rho_xi = θ[22]; sigma_xi = θ[23]

    # Bounds check
    (!(0<epsY<5)||!(0<epsM<2)||ilabcosts<=0||kappaV<=0||abs(rho_om)>=1||
     sigma_om<0||abs(rho_A)>=1||any(<(0), isigma_tfp)||abs(rho_pvstar)>=1||
     sigma_pvstar<0||abs(rho_xi)>=1||sigma_xi<0) && return NAN46, false

    # Update parameters
    set_param!(context,"ilabcosts",ilabcosts); set_param!(context,"kappaV",kappaV)
    set_param!(context,"rho_om1",rho_om);      set_param!(context,"sigma_om",sigma_om)
    set_param!(context,"rho_tfp1",rho_A);      set_param!(context,"rho_pvstar",rho_pvstar)
    set_param!(context,"sigma_pvstar",sigma_pvstar); set_param!(context,"rho_xi",rho_xi)
    set_param!(context,"sigma_xi",sigma_xi)
    for i in 1:nsec; set_param!(context,"isigma_tfp_$(i)",isigma_tfp[i]); end

    # SS: cached recomputation
    epsY_prev = get_param_val(context,"epsY_1"); epsM_prev = get_param_val(context,"epsM_1")
    need_ss = abs(epsY-epsY_prev)>1e-8 || abs(epsM-epsM_prev)>1e-8
    for i in 1:nsec; set_param!(context,"epsY_$(i)",epsY); set_param!(context,"epsM_$(i)",epsM); end
    if need_ss
        ok = recompute_ss_cached!(context, epsY, epsM, baseline, endo_names)
        !ok && return NAN46, false
    end

    # Klein solve (per-thread cache)
    success, T, R = _resolve_cached!(context, θ)
    !success && return NAN46, false

    n_exo = size(R, 2)
    sr = context.models[1].i_bkwrd_b
    n_state = length(sr)

    # Get pre-allocated scratch arrays
    sc = _get_scratch(n_exo, n_state, endo_names)

    # Build Σe in-place (zero and fill, no allocation)
    fill!(sc.Σe, 0.0)
    sc.Σe[1, 1] = 1.0     # eps_om
    sc.Σe[4, 4] = 1.0     # eps_pvstar
    for i in 5:min(16, n_exo); sc.Σe[i, i] = 1.0; end
    if n_exo >= 17; sc.Σe[17, 17] = 1.0; end

    # Lyapunov: Q_lyap = R[sr,:] * Σe * R[sr,:]'
    Tsr = T[sr, :]
    Rsr = R[sr, :]
    # Rsr is n_state × n_exo (78×17), Σe is n_exo × n_exo (17×17)
    # Product Rsr * Σe is n_state × n_exo → use a temporary of that size
    RΣ = Rsr * sc.Σe   # 78×17 (small, acceptable allocation)
    mul!(sc.Q_lyap, RΣ, Rsr')   # 78×78 = (78×17) × (17×78)
    # Symmetrize
    n_s = n_state
    @inbounds for j in 1:n_s, i in 1:j-1
        avg = 0.5 * (sc.Q_lyap[i,j] + sc.Q_lyap[j,i])
        sc.Q_lyap[i,j] = avg; sc.Q_lyap[j,i] = avg
    end

    # Solve P = A*P*A' + Q in-place
    local_dlyap_inplace!(sc.P, Tsr, sc.Q_lyap;
                          tmp1=sc.tmp1, tmp2=sc.tmp2, Ak=sc.Ak)
    (any(diag(sc.P) .< -1e-10) || any(isnan.(sc.P))) && return NAN46, false

    # Full covariance
    Γ  = T * sc.P * T' + R * sc.Σe * R'
    Γ  = (Γ + Γ') / 2
    Γ1 = T * Tsr * (sc.P * T' + Rsr * sc.Σe * R')

    # Extract moments using pre-computed index dict
    ys = context.results.model_results[1].trends.endogenous_steady_state
    ei = sc.endo_idx

    @inline function pstd(vn::String)
        idx = get(ei, vn, 0)
        idx == 0 && return 0.0
        sqrt(max(Γ[idx,idx], 0.0)) / max(abs(ys[idx]), 1e-12)
    end

    @inline function xcorr(v1::String, v2::String)
        i1 = get(ei, v1, 0); i2 = get(ei, v2, 0)
        (i1 == 0 || i2 == 0) && return NaN
        d = sqrt(max(Γ[i1,i1], 0.0) * max(Γ[i2,i2], 0.0))
        d < 1e-15 ? 0.0 : clamp(Γ[i1,i2]/d, -1.0, 1.0)
    end

    i_Q = get(ei, "Q", 0)
    acQ = (i_Q > 0 && Γ[i_Q,i_Q] > 1e-15) ? Γ1[i_Q,i_Q]/Γ[i_Q,i_Q] : NaN

    # Fill moment vector (pre-allocated in scratch)
    m = sc.moments
    for i in 1:nsec; m[i]      = pstd("Y_$(i)");  end
    for i in 1:nsec; m[12+i]   = pstd("PH_$(i)"); end
    for i in 1:nsec; m[24+i]   = pstd("L_$(i)");  end
    m[37] = pstd("GDP"); m[38] = pstd("pi")
    m[39] = xcorr("GDP","pi"); m[40] = baseline.ombar_val
    m[41] = pstd("Q"); m[42] = acQ; m[43] = xcorr("GDP","Q")

    # Rank correlations
    dY = baseline.data_std_Y; dPH = baseline.data_std_PH; dL = baseline.data_std_L
    std_Y  = @view m[1:12]; std_PH = @view m[13:24]; std_L = @view m[25:36]
    vy = [isfinite(std_Y[i]) && isfinite(dY[i]) for i in 1:nsec]
    vp = [isfinite(std_PH[i]) && isfinite(dPH[i]) for i in 1:nsec]
    vl = [isfinite(std_L[i]) && isfinite(dL[i]) for i in 1:nsec]
    m[44] = sum(vy) >= 3 ? safe_spearman(collect(std_Y[vy]), dY[vy]) : 0.0
    m[45] = sum(vp) >= 3 ? safe_spearman(collect(std_PH[vp]), dPH[vp]) : 0.0
    m[46] = sum(vl) >= 3 ? safe_spearman(collect(std_L[vl]), dL[vl]) : 0.0

    return copy(m), true   # copy so caller gets independent vector
end


# =========================================================================== #
#  8.  MAIN ESTIMATION FUNCTION                                                #
# =========================================================================== #

function smm_run_v2(context::Dynare.Context; endo_names_override=nothing)

    endo_names = if endo_names_override !== nothing && length(endo_names_override) >= 400
        endo_names_override
    else
        Dynare.get_endogenous(context.symboltable)
    end
    @printf "  endo_names: %d variables\n" length(endo_names)

    baseline = build_baseline(context, endo_names, y_d, p_d, l_d, d_TBGDP, d_omG)

    n_threads_active = Threads.nthreads()
    contexts_th = n_threads_active > 1 ?
        [deepcopy(context) for _ in 1:n_threads_active] : [context]
    baselines_th = [deepcopy(baseline) for _ in 1:length(contexts_th)]

    W = build_weighting_matrix_v2(data_moments)
    @printf "  Weighting: proportional + 5x rank correlations.\n\n"

    # Initial θ
    θ0 = default_theta0(context)
    θ0[6]  = max(θ0[6],  0.05)
    for i in 1:12
        θ0[7+i] = max(θ0[7+i], clamp(y_d[i] * 0.4, LB_V2[7+i], UB_V2[7+i]))
    end
    θ0[21] = max(θ0[21], 0.05)
    θ0[23] = max(θ0[23], 0.02)

    θ_warm = load_warm_start(N_THETA_V2)
    θ0 = something(θ_warm, θ0)
    θ0 = clamp.(θ0, LB_V2, UB_V2)

    # Pre-flight
    @printf "=== PRE-FLIGHT CHECK ===\n"
    m_test, ok_test = smm_model_moments_v2(θ0, context, baseline, endo_names)
    if !ok_test || any(isnan, m_test)
        @printf "  FAILED: ok=%s  NaN indices: %s\n" ok_test string(findall(isnan, m_test))
        error("Pre-flight failed.")
    end
    obj_test = dot(data_moments .- m_test, W * (data_moments .- m_test))
    @printf "  obj(θ₀) = %.6f\n" obj_test

    ψ0 = data_moments .- m_test
    obj_Y  = dot(ψ0[1:12],  W[1:12,1:12]  * ψ0[1:12])
    obj_PH = dot(ψ0[13:24], W[13:24,13:24] * ψ0[13:24])
    obj_L  = dot(ψ0[25:36], W[25:36,25:36] * ψ0[25:36])
    obj_agg= dot(ψ0[37:43], W[37:43,37:43] * ψ0[37:43])
    obj_rk = dot(ψ0[44:46], W[44:46,44:46] * ψ0[44:46])
    @printf "  Decomp: Y=%.2f PH=%.2f L=%.2f Agg=%.2f Rank=%.2f\n" obj_Y obj_PH obj_L obj_agg obj_rk
    @printf "=== PRE-FLIGHT PASSED ===\n\n"

    # CMA-ES
    max_evals = 30_000
    @printf "--- CMA-ES v2 ---\n"
    @printf "  %d params | %d moments | max %d evals | %d threads\n" N_THETA_V2 46 max_evals n_threads_active

    best_θ       = Ref(clamp.(θ0, LB_V2, UB_V2))
    best_obj     = Ref(Inf)
    best_lock    = ReentrantLock()
    fail_count   = Threads.Atomic{Int}(0)
    eval_count   = Threads.Atomic{Int}(0)
    t_start      = Ref(time())
    t_last_print = Ref(time())

    _tid_cma() = min(Threads.threadid(), length(contexts_th))

    obj_fn = θ -> begin
        tid = _tid_cma()
        ctx_th = contexts_th[tid]
        base_th = baselines_th[tid]

        moments, ok = smm_model_moments_v2(θ, ctx_th, base_th, endo_names)
        obj = if !ok || any(isnan, moments)
            1e8
        else
            ψ = data_moments .- moments
            dot(ψ, W * ψ)
        end

        Threads.atomic_add!(eval_count, 1)

        if isfinite(obj) && obj < best_obj[]
            lock(best_lock) do
                if isfinite(obj) && obj < best_obj[]
                    best_obj[] = obj
                    best_θ[]   = copy(θ)
                end
            end
        end
        obj >= 1e7 && Threads.atomic_add!(fail_count, 1)

        if eval_count[] % 200 == 0 && tid == 1
            t_now = time()
            ms_per = 1000.0 * (t_now - t_last_print[]) / 200
            total  = _KLEIN_HITS[] + _KLEIN_MISSES[]
            k_pct  = total > 0 ? round(Int, 100*_KLEIN_HITS[]/total) : 0
            @printf "  %-6d  obj=%-12.4f  fail=%-6d  %.0fms/eval  Klein=%d%%\n" eval_count[] best_obj[] fail_count[] ms_per k_pct
            flush(stdout)
            t_last_print[] = t_now
        end
        obj
    end

    # Per-parameter insigma: (UB-LB)/10 — tighter initial search
    insigma_scalar = mean((UB_V2 .- LB_V2) ./ 10)

    CMAEvolutionStrategy.minimize(
        obj_fn,
        clamp.(θ0, LB_V2, UB_V2),
        insigma_scalar;
        lower     = LB_V2,
        upper     = UB_V2,
        maxiter   = max_evals,
        ftol      = 1e-6,
        xtol      = 1e-6,
        seed      = 42,
        verbosity = 0,
        multi_threading = n_threads_active > 1,
    )

    θ_hat      = clamp.(best_θ[], LB_V2, UB_V2)
    total_time = time() - t_start[]
    total_kl   = _KLEIN_HITS[] + _KLEIN_MISSES[]
    @printf "\nCMA-ES done: %d evals | %.1fs | %.0fms/eval | Klein cache=%d%%\n\n" eval_count[] total_time (1000*total_time/max(1,eval_count[])) round(Int,100*_KLEIN_HITS[]/max(1,total_kl))

    # Final evaluation
    moments_hat, ok_final = smm_model_moments_v2(θ_hat, context, baseline, endo_names)
    obj_hat = ok_final ? dot(data_moments .- moments_hat, W * (data_moments .- moments_hat)) : NaN
    ψ_hat = data_moments .- moments_hat

    # Results
    @printf "\n%s\n  SMM v2 RESULTS\n%s\n\n" repeat("=",60) repeat("=",60)
    @printf "  Objective: %.6f\n\n" obj_hat

    @printf "  %-6s  %-18s  %10s  %10s\n" "Idx" "Parameter" "Initial" "Estimate"
    @printf "  %s\n" repeat("-", 50)
    for k in 1:N_THETA_V2
        if k == 4
            @printf "  %3d  %-18s  %10.4f  %10.4f  [log]\n" k PARAM_LABELS_V2[k] θ0[k] θ_hat[k]
            @printf "  %3s  %-18s  %10.2e  %10.2e  [level]\n" "--" "kappaV" exp(θ0[k]) exp(θ_hat[k])
        else
            @printf "  %3d  %-18s  %10.4f  %10.4f\n" k PARAM_LABELS_V2[k] θ0[k] θ_hat[k]
        end
    end

    @printf "\n  %-36s  %9s  %9s  %9s  %6s\n" "Moment" "Data" "Model" "Diff" "Rel%"
    @printf "  %s\n" repeat("-", 75)
    for (i, nm) in enumerate(MOMENT_NAMES)
        re = abs(data_moments[i]) > 1e-6 ? abs(ψ_hat[i])/abs(data_moments[i])*100 : 0.0
        @printf "  %-36s  %9.5f  %9.5f  %+9.5f  %5.1f%%\n" nm data_moments[i] moments_hat[i] ψ_hat[i] re
    end

    # Decomposition
    obj_Y  = dot(ψ_hat[1:12],  W[1:12,1:12]  * ψ_hat[1:12])
    obj_PH = dot(ψ_hat[13:24], W[13:24,13:24] * ψ_hat[13:24])
    obj_L  = dot(ψ_hat[25:36], W[25:36,25:36] * ψ_hat[25:36])
    obj_agg= dot(ψ_hat[37:43], W[37:43,37:43] * ψ_hat[37:43])
    obj_rk = dot(ψ_hat[44:46], W[44:46,44:46] * ψ_hat[44:46])
    @printf "\n  Decomp: Y=%.4f(%.1f%%) PH=%.4f(%.1f%%) L=%.4f(%.1f%%) Agg=%.4f(%.1f%%) Rank=%.4f(%.1f%%)\n" obj_Y 100obj_Y/obj_hat obj_PH 100obj_PH/obj_hat obj_L 100obj_L/obj_hat obj_agg 100obj_agg/obj_hat obj_rk 100obj_rk/obj_hat

    # Save
    df_results = DataFrame(
        param  = vcat(PARAM_LABELS_V2, fill("", 46 - N_THETA_V2)),
        theta  = vcat(θ_hat, fill(NaN, 46 - N_THETA_V2)),
        moment = MOMENT_NAMES,
        data   = data_moments,
        model  = moments_hat,
        diff   = ψ_hat,
    )
    CSV.write(joinpath(DATA_DIR, "smm_results_v2.csv"), df_results)

    df_est = DataFrame(
        param = vcat(["ilabcosts","epsY","epsM","log_kappaV",
                      "rho_om","sigma_om","rho_A"],
                     ["isigma_tfp_$(i)" for i in 1:NSEC],
                     ["rho_pvstar","sigma_pvstar","rho_xi","sigma_xi"]),
        value = θ_hat,
        obj_hat = vcat([obj_hat], fill(NaN, N_THETA_V2-1)),
    )
    CSV.write(joinpath(DATA_DIR, "smm_estimates.csv"), df_est)
    CSV.write(joinpath(DATA_DIR, "smm_checkpoint.csv"), df_est)

    @printf "\nResults saved. Re-run main_SOE_gap.jl (EXERCISE=0) to apply.\n\n"
    return θ_hat, obj_hat, moments_hat
end
