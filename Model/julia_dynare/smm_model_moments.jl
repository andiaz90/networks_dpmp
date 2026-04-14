"""
smm_model_moments.jl
====================
Compute theoretical model moments for a given parameter vector θ.

Uses the linearised state-space representation to compute exact unconditional
second moments via the discrete Lyapunov equation, then maps them to the
46-element target vector used in estimation.

RESOLVE APPROACH (ARM-compatible):
  Dynare.jl's compute_first_order_solution!() uses LAPACK gees with a select
  callback — not available on Apple Silicon (aarch64).
  We bypass this by calling the compiled Dynare Jacobian functions directly,
  then solving the Blanchard-Kahn conditions with GenericSchur.jl (pure Julia).
"""

using LinearAlgebra, Statistics, StatsBase, NLsolve
using GenericSchur    # pure-Julia QZ with eigenvalue ordering — works on ARM

include("steady_ntwsoe_system.jl")
include("steady_ntwsoe.jl")
include("utils.jl")

# =========================================================================== #
#  PARAMETER ACCESS                                                            #
# =========================================================================== #

# Cache: context → Dict{name => index}
const _param_idx_cache = IdDict{Any, Dict{String,Int}}()

function param_idx(context, name::String)
    if !haskey(_param_idx_cache, context)
        _param_idx_cache[context] = _build_param_cache(context)
    end
    return get(_param_idx_cache[context], name, nothing)
end

function _build_param_cache(context)
    cache = Dict{String,Int}()
    # Try Dynare.jl symbol table API
    st = context.symboltable
    try
        for fname in fieldnames(typeof(st))
            obj = getfield(st, fname)
            if obj isa AbstractDict
                for (k, v) in obj
                    if hasproperty(v, :name) && hasproperty(v, :index)
                        cache[String(v.name)] = v.index
                    elseif hasproperty(v, :symboltype)
                        # SymbolTable entry
                        try
                            nm  = String(k)
                            idx = v.index
                            cache[nm] = idx
                        catch; end
                    end
                end
            end
        end
    catch; end
    return cache
end

function set_param!(context, name::String, val::Real)
    idx = param_idx(context, name)
    idx === nothing && return
    context.models[1].params[idx] = Float64(val)
end

function get_param_val(context, name::String)
    idx = param_idx(context, name)
    idx === nothing && return NaN
    return context.models[1].params[idx]
end


# =========================================================================== #
#  RESOLVE FIRST-ORDER SOLUTION  (pure Julia via GenericSchur)                #
# =========================================================================== #

"""
    resolve_first_order!(context) -> (success, g1_1, g1_2, Sigma_e)

Recompute the first-order decision rule after parameter updates, using
GenericSchur.jl for the ordered QZ decomposition — pure Julia, no LAPACK
gees callback, works on Apple Silicon (ARM).

Implements Klein (2000): "Using the generalized Schur form to solve a
multivariate linear rational expectations model."
"""
function resolve_first_order!(context)
    try
        # ---- Step 1: try Dynare.jl's own re-solve (works on Intel) -------
        for fn_name in [:compute_first_order_solution!,
                        :first_order_solution!]
            isdefined(Dynare, fn_name) || continue
            try
                getfield(Dynare, fn_name)(context)
                mr  = context.results.model_results[1]
                lre = mr.linearrationalexpectations
                return true, Matrix{Float64}(lre.g1_1),
                             Matrix{Float64}(lre.g1_2),
                             context.models[1].Sigma_e
            catch; end   # silent — Dynare exceptions can't be string-ified on ARM
        end

        # ---- Step 2: pure-Julia Klein (2000) solver ----------------------
        return _klein_solve(context)

    catch
        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    end
end

"""
    _klein_solve(context) -> (success, g1_1, g1_2, Sigma_e)

Pure-Julia first-order perturbation solver using GenericSchur.jl.
Calls the compiled Dynare Jacobian functions directly.
"""
function _klein_solve(context)
    try
        m = context.models[1]

        # Retrieve steady state and parameters
        ss    = context.results.model_results[1].trends.endogenous_steady_state
        n_ys  = length(ss)          # number of endogenous variables
        exo   = zeros(m.exo_nbr)

        # ---- Evaluate the dynamic Jacobian at the steady state -----------
        # Dynare compiled the model Jacobian into SparseDynamicG1!
        # We call it via the context's dynamic evaluation.
        # g1 has size n_eq × (n_ys_lag + n_ys_now + n_ys_lead + n_exo)
        n_aux   = m.n_bkwrd          # number of backward variables (states)
        n_fwrd  = m.n_fwrd           # number of forward variables
        n_both  = m.n_both           # both lagged and leading
        n_static = m.n_static

        # Build dynamic Jacobian via Dynare.jl's compute_jacobian
        # This works because it only evaluates at a given point, no Schur needed
        J = Dynare.get_dynamic_jacobian!(context)

        if J === nothing || size(J, 1) == 0
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end

        n_eq = size(J, 1)
        n_cols = size(J, 2)

        # Partition Jacobian: [A | B | C | D] where
        #   A = df/dy_{t-1}  (n_eq × n_bkwrd_b)
        #   B = df/dy_t      (n_eq × n_ys)
        #   C = df/dy_{t+1}  (n_eq × n_fwrd_b)
        #   D = df/dε_t      (n_eq × n_exo)
        bkwrd_b = m.i_bkwrd_b
        fwrd_b  = m.i_fwrd_b
        n_bk    = length(bkwrd_b)
        n_fw    = length(fwrd_b)

        A_mat = Matrix{Float64}(J[:, 1:n_bk])
        B_mat = Matrix{Float64}(J[:, n_bk+1:n_bk+n_ys])
        C_mat = Matrix{Float64}(J[:, n_bk+n_ys+1:n_bk+n_ys+n_fw])
        D_mat = Matrix{Float64}(J[:, n_bk+n_ys+n_fw+1:end])

        # ---- Klein (2000) QZ decomposition -------------------------------
        # System: [C 0; 0 I] E[z_{t+1}] = [-B -A; I 0] z_t + [-D; 0] ε
        n_z  = n_bk + n_fw
        AA   = [C_mat zeros(n_eq, n_bk); zeros(n_bk, n_fw) I(n_bk)]
        BB   = [-B_mat -A_mat; I(n_bk) zeros(n_bk, n_bk)]
        
        # Generalized Schur (QZ) via GenericSchur.jl — works on ARM!
        F = GenericSchur.schur(AA, BB)
        S = F.S; T = F.T; Q = F.Q; Z = F.Z

        # Sort: stable eigenvalues (|T_ii/S_ii| < 1) to the top
        λ = [abs(T[i,i]) < 1e-14*abs(S[i,i]) ? 0.0 :
             abs(S[i,i]) < 1e-14 ? Inf : abs(T[i,i]/S[i,i])
             for i in 1:size(S,1)]
        select = λ .< 1.0    # stable = |eigenvalue| < 1

        n_stable = sum(select)
        if n_stable != n_bk
            # Blanchard-Kahn conditions not satisfied
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end

        # Reorder: stable eigenvalues first
        F2 = GenericSchur.ordschur(F, select)
        Z2 = F2.Z

        # Partition Z
        Z11 = Z2[1:n_bk, 1:n_bk]
        Z21 = Z2[n_bk+1:end, 1:n_bk]

        if abs(det(Z11)) < 1e-10
            return false, zeros(0,0), zeros(0,0), zeros(0,0)
        end

        # Decision rule for forward variables as function of backward: gx
        gx = real(Z21 / Z11)

        # Full state-space: g1_1 and g1_2 need to be in Dynare's ordering
        # g1_1[i,j]: response of variable i to state j
        # g1_2[i,k]: response of variable i to shock k
        # Build them in the ordering Dynare expects (endo declaration order)
        n_endo = n_ys
        n_state = n_bk
        n_shocks = m.exo_nbr

        g1_1_out = zeros(n_endo, n_state)
        g1_2_out = zeros(n_endo, n_shocks)

        # State equations: backward variables
        g1_1_out[bkwrd_b, :] = I(n_state)   # y_{t,bk} = x_t (by definition)
        # Forward variables driven by states
        g1_1_out[fwrd_b, :]  = gx

        # Impact of shocks: solve for g1_2 from model equations
        # C*gx*g1_2_bk + B*g1_2 + D = 0  (from model equation at period t)
        # For now use analytical formula when possible; fall back to pseudo-inverse
        lhs = B_mat + C_mat * gx * Matrix(I, n_state, n_state)[1:n_fw, :]
        if size(lhs, 1) == size(lhs, 2) && abs(det(lhs)) > 1e-12
            g1_2_full = -lhs \ D_mat
        else
            g1_2_full = -pinv(Matrix(lhs)) * D_mat
        end
        g1_2_out[1:size(g1_2_full,1), 1:size(g1_2_full,2)] .= g1_2_full

        # Update context with new decision rule
        lre = context.results.model_results[1].linearrationalexpectations
        try
            lre.g1_1 .= g1_1_out
            lre.g1_2 .= g1_2_out
        catch; end   # read-only struct — still return the new matrices

        Σe = context.models[1].Sigma_e
        return true, g1_1_out, g1_2_out, Σe

    catch
        return false, zeros(0,0), zeros(0,0), zeros(0,0)
    end
end


# =========================================================================== #
#  RECOMPUTE STEADY STATE  (when epsY or epsM change)                         #
# =========================================================================== #

function recompute_ss!(context, epsY, epsM, baseline, endo_names)
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

    ss_vec = context.results.model_results[1].trends.endogenous_steady_state
    get_ss(nm) = let idx = findfirst(==(nm), endo_names)
                     idx === nothing ? 1.0 : ss_vec[idx]
                 end

    pH_guess = [get_ss("PH_$(i)") for i in 1:nsec]
    x0 = [pH_guess; get_ss("w"); get_ss("Q"); get_ss("C")]

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

    inner = nlsolve(
        (F, x) -> steady_ntwsoe_system!(F, x, alpha_vec, alphaV_vec, beta_mat,
                                         MCi_ss, PMi_ss, PL_ss, PV_ss,
                                         CHi_ss, Xi_ss, modepsY, modepsM,
                                         A_vec, pH_ss),
        [max.(MCi_ss./PMi_ss,1e-20); max.(MCi_ss./PL_ss,1e-20);
         max.(MCi_ss./PV_ss,1e-20); fill(0.2,nsec)];
        ftol=1e-10, method=:trust_region, show_trace=false
    )
    !converged(inner) && return false

    Yi_ss = inner.zero[3*nsec+1:4*nsec]
    M_ss  = inner.zero[1:nsec]

    IMP_tot = sum(inner.zero[2*nsec+1:3*nsec]) + sum(CFi_ss)
    TB_ss   = PX_ss*X_ss - PV_ss*IMP_tot
    GDP_ss  = C_ss + TB_ss
    Pistar  = 1.0
    r_star  = Pistar / beta_val
    Bstar   = -TB_ss / (Q_ss * (1 - r_star/Pistar))
    bbar    = Q_ss * Bstar / GDP_ss

    set_param!(context, "bbar",     bbar)
    set_param!(context, "Y_ss",     sum(Yi_ss))
    set_param!(context, "M_tot_ss", sum(M_ss))
    for i in 1:nsec
        set_param!(context, "PL_ss$(i)", w_ss)
        set_param!(context, "epsY_$(i)", epsY)
        set_param!(context, "epsM_$(i)", epsM)
    end

    # Update steady-state vector in context
    ss_mut = context.results.model_results[1].trends.endogenous_steady_state
    upd(nm, val) = let idx = findfirst(==(nm), endo_names)
                       idx !== nothing && (ss_mut[idx] = val)
                   end
    upd("w", w_ss); upd("Q", Q_ss); upd("C", C_ss)
    upd("GDP", GDP_ss); upd("TB", TB_ss)
    for i in 1:nsec
        upd("PH_$(i)", pH_ss[i]); upd("Y_$(i)", Yi_ss[i])
        upd("L_$(i)",  inner.zero[nsec+i])
    end
    return true
end


# =========================================================================== #
#  MOMENT NAMES                                                                #
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
#  MAIN MOMENT FUNCTION                                                        #
# =========================================================================== #

function smm_model_moments(θ, context, baseline, endo_names)
    nsec = baseline.nsec
    NAN46 = fill(NaN, 46)

    ilabcosts  = θ[1]; epsY = θ[2]; epsM = θ[3]
    kappaV     = exp(θ[4]); rho_om = θ[5]; sigma_om = θ[6]
    rho_A      = θ[7]; isigma_tfp = θ[8:19]
    rho_pvstar = θ[20]; sigma_pvstar = θ[21]
    rho_xi     = θ[22]; sigma_xi = θ[23]

    (!(0 < epsY < 5) || !(0 < epsM < 2) || ilabcosts <= 0 || kappaV <= 0 ||
     abs(rho_om) >= 1 || sigma_om < 0 || abs(rho_A) >= 1 ||
     any(isigma_tfp .< 0) || abs(rho_pvstar) >= 1 || sigma_pvstar < 0 ||
     abs(rho_xi) >= 1 || sigma_xi < 0) && return NAN46, false

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
    try
        context.models[1].Sigma_e[4,  4]  = 1.0
        context.models[1].Sigma_e[17, 17] = 1.0
    catch; end

    epsY_prev = get_param_val(context, "epsY_1")
    epsM_prev = get_param_val(context, "epsM_1")
    need_ss = (abs(epsY - epsY_prev) > 1e-8) || (abs(epsM - epsM_prev) > 1e-8)
    for i in 1:nsec
        set_param!(context, "epsY_$(i)", epsY)
        set_param!(context, "epsM_$(i)", epsM)
    end
    if need_ss
        ok = recompute_ss!(context, epsY, epsM, baseline, endo_names)
        !ok && return NAN46, false
    end

    success, T, R, Σe = resolve_first_order!(context)
    !success && return NAN46, false

    # Lyapunov
    state_rows = context.models[1].i_bkwrd_b
    A_rc = T[state_rows, :]; B_rc = R[state_rows, :]
    P_st = local_dlyap(A_rc, B_rc * Σe * B_rc')
    (any(diag(P_st) .< -1e-10) || any(isnan.(P_st))) && return NAN46, false
    P_st = (P_st + P_st') / 2
    Γ    = T * P_st * T' + R * Σe * R'
    Γ    = (Γ + Γ') / 2
    Γ_1  = T * A_rc * (P_st * T' + B_rc * Σe * R')

    ys = context.results.model_results[1].trends.endogenous_steady_state
    endo_idx = Dict(nm => i for (i,nm) in enumerate(endo_names))
    pstd(vn) = let idx = get(endo_idx, vn, nothing)
                   idx === nothing ? 0.0 :
                   sqrt(max(Γ[idx,idx],0.0)) / max(abs(ys[idx]),1e-12)
               end
    xcorr(v1,v2) = let i1=get(endo_idx,v1,nothing), i2=get(endo_idx,v2,nothing)
                       (i1===nothing||i2===nothing) ? NaN :
                       let d=sqrt(max(Γ[i1,i1],0.0)*max(Γ[i2,i2],0.0))
                           d < 1e-15 ? 0.0 : clamp(Γ[i1,i2]/d,-1.0,1.0)
                       end
                   end

    std_Y  = [pstd("Y_$(i)")  for i in 1:nsec]
    std_PH = [pstd("PH_$(i)") for i in 1:nsec]
    std_L  = [pstd("L_$(i)")  for i in 1:nsec]
    i_Q    = get(endo_idx, "Q", nothing)
    autocorr_Q = (i_Q!==nothing && Γ[i_Q,i_Q]>1e-15) ?
                  Γ_1[i_Q,i_Q]/Γ[i_Q,i_Q] : NaN

    d_Y  = baseline.data_std_Y; d_PH = baseline.data_std_PH; d_L = baseline.data_std_L
    valid_y = isfinite.(std_Y)  .& isfinite.(d_Y)
    valid_p = isfinite.(std_PH) .& isfinite.(d_PH)
    valid_l = isfinite.(std_L)  .& isfinite.(d_L)
    rho_Y  = sum(valid_y)>=3 ? safe_spearman(std_Y[valid_y],  d_Y[valid_y])  : 0.0
    rho_PH = sum(valid_p)>=3 ? safe_spearman(std_PH[valid_p], d_PH[valid_p]) : 0.0
    rho_L  = sum(valid_l)>=3 ? safe_spearman(std_L[valid_l],  d_L[valid_l])  : 0.0

    moments = [std_Y; std_PH; std_L;
               pstd("GDP"); pstd("pi"); xcorr("GDP","pi");
               baseline.ombar_val;
               pstd("Q"); autocorr_Q; xcorr("GDP","Q");
               rho_Y; rho_PH; rho_L]
    return moments, true
end
