"""
smm_model_moments.jl
====================
Compute theoretical model moments for a given parameter vector θ.

Uses Dynare.jl's linearised state-space to obtain exact unconditional second
moments via the discrete Lyapunov equation, then maps them to the 46-element
target vector used in the paper's SMM estimation.

Translated from smm_model_moments.m (MATLAB / Dynare).

THETA ORDERING (23 elements)
  θ[1]      ilabcosts       inverse aggregate labour adjustment cost
  θ[2]      epsY            elast. of subst. in production (common)
  θ[3]      epsM            elast. of subst. between materials (common)
  θ[4]      log(kappaV)     log of import price adj. cost
  θ[5]      rho_om          AR persistence, goods-services shock
  θ[6]      sigma_om        std dev, goods-services shock
  θ[7]      rho_A           AR persistence, sectoral TFP shocks (common)
  θ[8:19]   isigma_tfp_i    std dev of each sector's TFP shock (i=1,...,12)
  θ[20]     rho_pvstar      AR persistence, import price shock
  θ[21]     sigma_pvstar    std dev of import price shock
  θ[22]     rho_xi          AR persistence, preference/demand shock
  θ[23]     sigma_xi        std dev of preference/demand shock
"""

using LinearAlgebra, Statistics, StatsBase, NLsolve, Dynare, Printf

include("steady_ntwsoe_system.jl")
include("steady_ntwsoe.jl")
include("utils.jl")

# =========================================================================== #
#  MOMENT NAMES (46 elements)                                                 #
# =========================================================================== #

const MOMENT_NAMES = vcat(
    ["std(Y_$(i))"  for i in 1:12],
    ["std(PH_$(i))" for i in 1:12],
    ["std(L_$(i))"  for i in 1:12],
    ["std(GDP)", "std(pi)", "corr(GDP,pi)",
     "mean goods expenditure share",
     "std(Q)", "autocorr(Q)", "corr(GDP,Q)",
     "rank corr: output (model vs data)",
     "rank corr: prices (model vs data)",
     "rank corr: labor  (model vs data)"],
)

# =========================================================================== #
#  PARAMETER HELPER                                                            #
# =========================================================================== #

"""
    param_idx(context, name) -> Int or nothing

Return the 1-based index of parameter `name` in the Dynare.jl context,
or `nothing` if not found.
"""
function param_idx(context::Dynare.Context, name::String)
    names = Dynare.get_parameters(context.symboltable)
    return findfirst(==(name), names)
end

"""
    set_param!(context, name, val)

Update a single named parameter in the Dynare.jl context.
"""
function set_param!(context::Dynare.Context, name::String, val::Real)
    idx = param_idx(context, name)
    idx === nothing && return   # silently skip unknown parameters
    context.models[1].params[idx] = Float64(val)
end

"""
    get_param_val(context, name) -> Float64

Read the current value of a named parameter from the Dynare.jl context.
"""
function get_param_val(context::Dynare.Context, name::String)
    idx = param_idx(context, name)
    idx === nothing && return NaN
    return context.models[1].params[idx]
end


# =========================================================================== #
#  RESOLVE THE LINEARISED MODEL                                                #
#                                                                              #
#  After updating context.models[1].params, call this to recompute the        #
#  first-order decision rule (g1_1 ≈ ghx, g1_2 ≈ ghu) without recompiling.  #
# =========================================================================== #

"""
    resolve_first_order!(context) -> (success::Bool, g1_1, g1_2, Sigma_e)

Re-solve the first-order perturbation solution with the current parameter
values in context.models[1].params.  Equivalent to Dynare/MATLAB's resol().

Returns (true, g1_1, g1_2, Sigma_e) on success, (false, ...) on failure.
"""
function resolve_first_order!(context::Dynare.Context)
    try
        # Dynare.jl API: recompute the first-order solution in-place.
        # This is the Julia equivalent of MATLAB's resol(0, M_, options_, ...).
        Dynare.compute_first_order_solution!(context)

        mr   = context.results.model_results[1]
        lre  = mr.linearrationalexpectations
        g1_1 = lre.g1_1   # n_endo × n_states  (≈ ghx)
        g1_2 = lre.g1_2   # n_endo × n_shocks  (≈ ghu)
        Σe   = context.models[1].Sigma_e
        return true, g1_1, g1_2, Σe
    catch e
        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    end
end


# =========================================================================== #
#  VARIANCE EXTRACTION HELPERS                                                 #
# =========================================================================== #

"""
    endo_dr_idx(context, endo_names, name) -> Int or nothing

Return the decision-rule row index (in Γ) corresponding to endogenous
variable `name`.  Dynare.jl uses `i_bkwrd_b` for state ordering; here we
simply find the column position of `name` in the endo_names vector.
"""
function endo_dr_idx(endo_names::Vector{String}, name::String)
    return findfirst(==(name), endo_names)
end

"""
    pct_std(Γ, idx, ss_val) -> Float64

Percentage std dev of variable at decision-rule row `idx`:
  sqrt(Γ[idx,idx]) / |ss_val|
"""
function pct_std(Γ::AbstractMatrix{<:Real}, idx::Union{Int,Nothing}, ss_val::Real)
    (idx === nothing || ss_val == 0) && return 0.0
    return sqrt(max(Γ[idx, idx], 0.0)) / max(abs(ss_val), 1e-12)
end

"""
    contemporaneous_corr(Γ, i, j) -> Float64

Pearson correlation between variables at DR rows i and j.
"""
function contemporaneous_corr(Γ::AbstractMatrix{<:Real},
                               i::Union{Int,Nothing}, j::Union{Int,Nothing})
    (i === nothing || j === nothing) && return 0.0
    denom = sqrt(max(Γ[i,i], 0.0) * max(Γ[j,j], 0.0))
    denom < 1e-15 && return 0.0
    return clamp(Γ[i,j] / denom, -1.0, 1.0)
end


# =========================================================================== #
#  RECOMPUTE STEADY STATE  (needed when epsY or epsM change)                  #
# =========================================================================== #

"""
    recompute_ss!(context, epsY, epsM, baseline, endo_names) -> Bool

Recompute the full non-linear steady state for new production elasticities
(epsY, epsM) and update both context.models[1].params and the model's
steady-state vector.  Returns true on success.

Mirrors the MATLAB recompute_ss() local function in smm_model_moments.m.
"""
function recompute_ss!(context::Dynare.Context,
                        epsY::Real, epsM::Real,
                        baseline::NamedTuple,
                        endo_names::Vector{String})
    nsec       = baseline.nsec
    modepsY    = fill(epsY, nsec)
    modepsM    = fill(epsM, nsec)
    alpha_vec  = baseline.modalpha
    alphaV_vec = baseline.modalphaV
    beta_mat   = baseline.modbeta
    varrho_val = baseline.modvarrho
    gammag_vec = baseline.modgammag
    gammas_vec = baseline.modgammas
    chiX_vec   = baseline.modchiX
    A_vec      = ones(nsec)
    om_g       = baseline.ombar_val
    om_s       = 1 - om_g
    PVstar_ss  = baseline.PVstar_ss
    epsilon    = baseline.epsilon_val
    sigmaH     = baseline.sigmaH_val
    etastar    = baseline.etastar_val
    omegaX     = baseline.omegaX_val
    Ystar      = baseline.ystar_ss_val
    gamma_val  = baseline.gamma_val
    chi_val    = baseline.chi_val
    psi_val    = baseline.psi_val
    beta_val   = baseline.beta_val
    tb_target  = baseline.tb_target

    # Current SS as warm-start for outer solve
    ss_vec     = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let idx = endo_dr_idx(endo_names, nm)
                     idx === nothing ? 1.0 : ss_vec[idx]
                 end

    pH_guess = [get_ss("PH_$(i)") for i in 1:nsec]
    x0 = [pH_guess; get_ss("w"); get_ss("Q"); get_ss("C")]

    # Outer NLsolve for (pH, w, Q, C)
    res = nlsolve(
        (F, x) -> F .= steady_ntwsoe(x, PVstar_ss, epsilon, varrho_val, sigmaH,
                                      gammag_vec, gammas_vec, om_g, om_s,
                                      chiX_vec, omegaX, etastar, Ystar,
                                      alpha_vec, alphaV_vec, beta_mat,
                                      modepsY, modepsM, gamma_val, chi_val,
                                      psi_val, A_vec, tb_target),
        x0; ftol=1e-12, method=:trust_region, show_trace=false
    )
    !converged(res) && return false

    pH_ss = res.zero[1:nsec]
    w_ss  = res.zero[nsec+1]
    Q_ss  = res.zero[nsec+2]
    C_ss  = res.zero[nsec+3]

    # Derive remaining SS quantities
    PL_ss  = fill(w_ss, nsec)
    PV_ss  = Q_ss * PVstar_ss
    MCi_ss = (epsilon-1)/epsilon .* pH_ss
    PMi_ss = (beta_mat * (pH_ss .^ (1 .- modepsM))) .^ (1 ./ (1 .- modepsM))

    P_ss   = (varrho_val .^ sigmaH .* pH_ss .^ (1-sigmaH)
             .+ (1 .- varrho_val) .^ sigmaH .* PV_ss .^ (1-sigmaH)) .^ (1/(1-sigmaH))
    p_g_ss = prod(P_ss .^ gammag_vec)
    p_s_ss = prod(P_ss .^ gammas_vec)
    C_g_ss = om_g * C_ss / p_g_ss
    C_s_ss = om_s * C_ss / p_s_ss
    C_gi_ss = gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss
    C_si_ss = gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss
    CHg_ss  = varrho_val .^ sigmaH .* (pH_ss ./ P_ss) .^ (-sigmaH) .* C_gi_ss
    CHs_ss  = varrho_val .^ sigmaH .* (pH_ss ./ P_ss) .^ (-sigmaH) .* C_si_ss
    CFg_ss  = (1 .- varrho_val) .^ sigmaH .* (PV_ss ./ P_ss) .^ (-sigmaH) .* C_gi_ss
    CFs_ss  = (1 .- varrho_val) .^ sigmaH .* (PV_ss ./ P_ss) .^ (-sigmaH) .* C_si_ss
    CHi_ss  = CHg_ss .+ CHs_ss
    CFi_ss  = CFg_ss .+ CFs_ss
    PX_ss   = prod(pH_ss .^ chiX_vec)
    X_ss    = omegaX * (PX_ss / Q_ss)^(-etastar) * Ystar
    Xi_ss   = chiX_vec .* X_ss .* PX_ss ./ pH_ss

    # Inner solve for (M, L, Vi, Yi)
    M_g  = get_ss.([("M_$(i)" for i in 1:nsec)...]) |> x -> max.(x, 0.1)
    L_g  = get_ss.([("L_$(i)" for i in 1:nsec)...]) |> x -> max.(x, 0.1)
    Vi_g = get_ss.([("V_$(i)" for i in 1:nsec)...]) |> x -> max.(x, 0.01)
    Yi_g = get_ss.([("Y_$(i)" for i in 1:nsec)...]) |> x -> max.(x, 0.1)

    inner = nlsolve(
        (F, x) -> steady_ntwsoe_system!(F, x, alpha_vec, alphaV_vec, beta_mat,
                                         MCi_ss, PMi_ss, PL_ss, PV_ss,
                                         CHi_ss, Xi_ss, modepsY, modepsM,
                                         A_vec, pH_ss),
        [M_g; L_g; Vi_g; Yi_g];
        ftol=1e-10, method=:trust_region, show_trace=false
    )

    if converged(inner)
        M_ss  = inner.zero[1:nsec]
        L_ss  = inner.zero[nsec+1:2*nsec]
        Vi_ss = inner.zero[2*nsec+1:3*nsec]
        Yi_ss = inner.zero[3*nsec+1:4*nsec]
    else
        # Analytical fallback (no intermediate goods loop)
        M_ss  = (MCi_ss./PMi_ss).^modepsY .* alpha_vec .* (CHi_ss.+Xi_ss)
        L_ss  = (MCi_ss./PL_ss).^modepsY .* (1 .- alpha_vec .- alphaV_vec) .* (CHi_ss.+Xi_ss)
        Vi_ss = (MCi_ss./PV_ss).^modepsY .* alphaV_vec .* (CHi_ss.+Xi_ss)
        Yi_ss = A_vec .* (alpha_vec.^(1 ./modepsY).*max.(M_ss,1e-20).^((modepsY.-1)./modepsY)
                        .+ alphaV_vec.^(1 ./modepsY).*max.(Vi_ss,1e-20).^((modepsY.-1)./modepsY)
                        .+ (1 .-alphaV_vec.-alpha_vec).^(1 ./modepsY).*max.(L_ss,1e-20).^((modepsY.-1)./modepsY)
                        ).^(modepsY./(modepsY.-1))
    end

    # Aggregate quantities
    IMP_tot_ss = sum(Vi_ss) + sum(CFi_ss)
    TB_ss      = PX_ss*X_ss - PV_ss*IMP_tot_ss
    GDP_ss     = C_ss + TB_ss
    N_ss       = sum(L_ss)
    Y_tot_ss   = sum(Yi_ss)
    M_tot_ss   = sum(M_ss)
    VA_ss      = sum(Yi_ss .- M_ss)
    Ctotg_ss   = sum(gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss)
    Ctots_ss   = sum(gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss)
    Ctot_ss    = Ctotg_ss + Ctots_ss
    Pistar_ss  = 1.0
    r_star_ss  = Pistar_ss / beta_val
    Bstar_ss   = -TB_ss / (Q_ss * (1 - r_star_ss/Pistar_ss))
    bbar_new   = Q_ss * Bstar_ss / GDP_ss

    # Update Dynare context: scalar model parameters
    for (nm, val) in [("Ctot_ss",Ctot_ss), ("Ctotg_ss",Ctotg_ss), ("Ctots_ss",Ctots_ss),
                       ("VA_ss",VA_ss), ("M_tot_ss",M_tot_ss), ("Y_ss",Y_tot_ss),
                       ("IMP_ss",IMP_tot_ss), ("bbar",bbar_new)]
        set_param!(context, nm, val)
    end
    for i in 1:nsec
        set_param!(context, "PL_ss$(i)", w_ss)
        set_param!(context, "epsY_$(i)", epsY)
        set_param!(context, "epsM_$(i)", epsM)
    end

    # Update steady-state vector in context
    ss_mut = context.results.model_results[1].trends.endogenous_steady_state
    function upd_ss(nm, val)
        idx = endo_dr_idx(endo_names, nm)
        idx !== nothing && (ss_mut[idx] = val)
    end
    upd_ss("w", w_ss); upd_ss("Q", Q_ss); upd_ss("C", C_ss)
    upd_ss("N", N_ss); upd_ss("GDP", GDP_ss); upd_ss("TB", TB_ss)
    for i in 1:nsec
        upd_ss("PH_$(i)", pH_ss[i]);  upd_ss("MC_$(i)", MCi_ss[i])
        upd_ss("PM_$(i)", PMi_ss[i]); upd_ss("PL_$(i)", PL_ss[i])
        upd_ss("P_$(i)",  P_ss[i]);   upd_ss("Y_$(i)",  Yi_ss[i])
        upd_ss("L_$(i)",  L_ss[i]);   upd_ss("M_$(i)",  M_ss[i])
        upd_ss("V_$(i)",  Vi_ss[i])
    end

    return true
end


# =========================================================================== #
#  MAIN MOMENT FUNCTION                                                        #
# =========================================================================== #

"""
    smm_model_moments(θ, context, baseline, endo_names) -> (moments, success)

Compute the 46-element theoretical moment vector for parameter vector θ.

Returns:
  moments :: Vector{Float64}  — 46 elements (NaN on failure)
  success :: Bool             — false if model failed to solve
"""
function smm_model_moments(θ::AbstractVector{<:Real},
                            context::Dynare.Context,
                            baseline::NamedTuple,
                            endo_names::Vector{String})

    nsec = baseline.nsec
    NAN46 = fill(NaN, 46)

    # ---- 1. Unpack θ and check feasibility -------------------------------- #
    ilabcosts  = θ[1]
    epsY       = θ[2]
    epsM       = θ[3]
    kappaV     = exp(θ[4])   # stored in log-space
    rho_om     = θ[5]
    sigma_om   = θ[6]
    rho_A      = θ[7]
    isigma_tfp = θ[8:19]
    rho_pvstar   = θ[20]
    sigma_pvstar = θ[21]
    rho_xi       = θ[22]
    sigma_xi     = θ[23]

    # Hard feasibility: return NaN if outside valid region
    if !(0 < epsY < 5) || !(0 < epsM < 2) || ilabcosts <= 0 || kappaV <= 0 ||
       abs(rho_om) >= 1 || sigma_om < 0 || abs(rho_A) >= 1 ||
       any(isigma_tfp .< 0) || abs(rho_pvstar) >= 1 || sigma_pvstar < 0 ||
       abs(rho_xi) >= 1 || sigma_xi < 0
        return NAN46, false
    end

    # ---- 2. Update Dynare parameters -------------------------------------- #
    set_param!(context, "ilabcosts",    ilabcosts)
    set_param!(context, "kappaV",       kappaV)
    set_param!(context, "rho_om1",      rho_om)
    set_param!(context, "sigma_om",     sigma_om)
    set_param!(context, "rho_tfp1",     rho_A)
    set_param!(context, "rho_pvstar",   rho_pvstar)
    set_param!(context, "sigma_pvstar", sigma_pvstar)
    set_param!(context, "rho_xi",       rho_xi)
    set_param!(context, "sigma_xi",     sigma_xi)
    for i in 1:nsec
        set_param!(context, "isigma_tfp_$(i)", isigma_tfp[i])
    end
    # Activate shock variances
    context.models[1].Sigma_e[4,  4]  = 1.0   # PVstar shock
    context.models[1].Sigma_e[17, 17] = 1.0   # preference (xi) shock

    # ---- 3. Recompute SS if epsY or epsM changed -------------------------- #
    epsY_prev = get_param_val(context, "epsY_1")
    epsM_prev = get_param_val(context, "epsM_1")
    need_ss   = (abs(epsY - epsY_prev) > 1e-8) || (abs(epsM - epsM_prev) > 1e-8)

    for i in 1:nsec
        set_param!(context, "epsY_$(i)", epsY)
        set_param!(context, "epsM_$(i)", epsM)
    end

    if need_ss
        ok = recompute_ss!(context, epsY, epsM, baseline, endo_names)
        !ok && return NAN46, false
    end

    # ---- 4. Re-solve first-order perturbation ----------------------------- #
    success, T, R, Σe = resolve_first_order!(context)
    !success && return NAN46, false

    # ---- 5. Discrete Lyapunov equation ------------------------------------ #
    # State rows: i_bkwrd_b gives 1-based indices of backward variables in T
    state_rows = context.models[1].i_bkwrd_b
    A_lyap = T[state_rows, :]
    B_lyap = R[state_rows, :]

    Q_lyap = B_lyap * Σe * B_lyap'
    Q_lyap = (Q_lyap + Q_lyap') / 2
    P_st   = local_dlyap(A_lyap, Q_lyap)

    # Check numerical health
    if any(diag(P_st) .< -1e-10) || any(isnan.(P_st)) || any(isinf.(P_st))
        return NAN46, false
    end
    P_st = (P_st + P_st') / 2

    # Contemporaneous variance-covariance of ALL endogenous variables
    Γ     = T * P_st * T' + R * Σe * R'
    Γ     = (Γ + Γ') / 2

    # Lag-1 cross-covariance (for autocorrelations)
    Γ_1   = T * A_lyap * (P_st * T' + B_lyap * Σe * R')

    # ---- 6. Extract steady-state vector ----------------------------------- #
    ss_vec = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let idx = endo_dr_idx(endo_names, nm)
                     idx === nothing ? 1.0 : ss_vec[idx]
                 end

    # ---- 7. Compute moments ----------------------------------------------- #
    # Percentage std devs for sectoral variables
    std_Y  = [pct_std(Γ, endo_dr_idx(endo_names, "Y_$(i)"),  get_ss("Y_$(i)"))  for i in 1:nsec]
    std_PH = [pct_std(Γ, endo_dr_idx(endo_names, "PH_$(i)"), get_ss("PH_$(i)")) for i in 1:nsec]
    std_L  = [pct_std(Γ, endo_dr_idx(endo_names, "L_$(i)"),  get_ss("L_$(i)"))  for i in 1:nsec]

    # Aggregate std devs
    m_std_GDP = pct_std(Γ, endo_dr_idx(endo_names, "GDP"), get_ss("GDP"))
    m_std_pi  = pct_std(Γ, endo_dr_idx(endo_names, "pi"),  get_ss("pi"))
    m_std_Q   = pct_std(Γ, endo_dr_idx(endo_names, "Q"),   get_ss("Q"))

    # Correlations
    m_corr_GDPpi = contemporaneous_corr(Γ,
        endo_dr_idx(endo_names, "GDP"),
        endo_dr_idx(endo_names, "pi"))

    # AR(1) autocorrelation of Q: Γ_1(i,i) / Γ(i,i)
    i_Q = endo_dr_idx(endo_names, "Q")
    m_autocorr_Q = if i_Q !== nothing && Γ[i_Q, i_Q] > 1e-15
        Γ_1[i_Q, i_Q] / Γ[i_Q, i_Q]
    else
        NaN
    end

    # corr(GDP, Q) contemporaneous
    m_corr_GDPQ = contemporaneous_corr(Γ,
        endo_dr_idx(endo_names, "GDP"),
        endo_dr_idx(endo_names, "Q"))

    m_omG = baseline.ombar_val   # SS calibration target (passive)

    # ---- 8. Rank correlations (cross-sectional, model vs data) ------------ #
    d_Y  = baseline.data_std_Y
    d_PH = baseline.data_std_PH
    d_L  = baseline.data_std_L

    valid_y = isfinite.(std_Y)  .& isfinite.(d_Y)
    valid_p = isfinite.(std_PH) .& isfinite.(d_PH)
    valid_l = isfinite.(std_L)  .& isfinite.(d_L)

    m_rho_Y  = sum(valid_y) >= 3 ? safe_spearman(std_Y[valid_y],  d_Y[valid_y])  : 0.0
    m_rho_PH = sum(valid_p) >= 3 ? safe_spearman(std_PH[valid_p], d_PH[valid_p]) : 0.0
    m_rho_L  = sum(valid_l) >= 3 ? safe_spearman(std_L[valid_l],  d_L[valid_l])  : 0.0

    # ---- 9. Stack 46-element moment vector -------------------------------- #
    moments = [
        std_Y;     # 1-12
        std_PH;    # 13-24
        std_L;     # 25-36
        m_std_GDP; # 37
        m_std_pi;  # 38
        m_corr_GDPpi; # 39
        m_omG;     # 40
        m_std_Q;   # 41
        m_autocorr_Q; # 42
        m_corr_GDPQ;  # 43
        m_rho_Y;   # 44
        m_rho_PH;  # 45
        m_rho_L;   # 46
    ]

    return moments, true
end
