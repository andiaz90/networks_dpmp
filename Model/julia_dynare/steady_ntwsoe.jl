"""
steady_ntwsoe.jl
================
Outer steady-state residuals for the NK-IOSOE model.
Translated from steady_ntwsoe.m (MATLAB).

Unknowns: x_vec = [pH_1,...,pH_nsec, w, Q, C]  (length nsec+3)

The function:
  1. Derives all prices and demands from (pH, w, Q, C)
  2. Calls steady_ntwsoe_system! (inner NLsolve) for (M, L, Vi, Yi)
  3. Returns nsec+3 market-clearing residuals

Returns residual vector of length nsec+3:
  residual[1:nsec]   = output market clearing (Yi = CHi + Xi + intermediate_use)
  residual[nsec+1]   = trade balance target   (TB/GDP = tb_target)
  residual[nsec+2]   = labor market clearing  (N = Σ L_i)
  residual[nsec+3]   = Cobb-Douglas consumption identity
"""

using NLsolve, Statistics

function steady_ntwsoe(
    x_vec    :: AbstractVector{<:Real},
    PVstar   :: Real,
    epsilon  :: Real,
    varrho_vec  :: Vector{<:Real},
    sigmaH   :: Real,
    gammag_vec  :: Vector{<:Real},
    gammas_vec  :: Vector{<:Real},
    om_g     :: Real,
    om_s     :: Real,
    chiX_vec :: Vector{<:Real},
    omegaX   :: Real,
    etastar  :: Real,
    Ystar    :: Real,
    alpha_vec   :: Vector{<:Real},
    alphaV_vec  :: Vector{<:Real},
    beta_mat :: Matrix{<:Real},
    epsY_vec :: Vector{<:Real},
    epsM_vec :: Vector{<:Real},
    GAMMA    :: Real,
    CHI      :: Real,
    PSI      :: Real,
    A_vec    :: Vector{<:Real},
    tb_target :: Real,
)
    nsec  = length(varrho_vec)
    # Clamp price-like unknowns to positive values so that fractional
    # exponentiation never receives a negative base (the NLsolve trust-region
    # solver can explore negative regions).
    pHvec = max.(x_vec[1:nsec], 1e-20)
    w     = max(x_vec[nsec+1], 1e-20)
    q     = max(x_vec[nsec+2], 1e-20)
    C     = max(x_vec[nsec+3], 1e-20)

    PL  = fill(w, nsec)
    PV  = q * PVstar

    # Marginal costs (Bertrand pricing: MC = (ε-1)/ε × PH)
    MCi = (epsilon - 1) / epsilon .* pHvec

    # Intermediate price indices: PMi[i] = (Σ_j beta[i,j] * pH_j^(1-ε_m_j))^(1/(1-ε_m_i))
    PMi = max.((beta_mat * (pHvec .^ (1 .- epsM_vec))) .^ (1 ./ (1 .- epsM_vec)), 1e-20)

    # Consumer price of each sector (Armington aggregator)
    pvec = max.((varrho_vec .^ sigmaH .* pHvec .^ (1 - sigmaH)
          .+ (1 .- varrho_vec) .^ sigmaH .* PV .^ (1 - sigmaH)) .^ (1 / (1 - sigmaH)), 1e-20)

    # Goods / services price indices (Cobb-Douglas)
    p_g = prod(pvec .^ gammag_vec)
    p_s = prod(pvec .^ gammas_vec)

    # Goods / services demands
    C_g = om_g * C / p_g
    C_s = om_s * C / p_s

    # Sectoral consumer demands within each basket
    C_gi = gammag_vec .* (p_g ./ pvec) .* C_g
    C_si = gammas_vec .* (p_s ./ pvec) .* C_s

    # Home vs. foreign split in consumption
    CH_gi = varrho_vec .^ sigmaH .* (pHvec ./ pvec) .^ (-sigmaH) .* C_gi
    CH_si = varrho_vec .^ sigmaH .* (pHvec ./ pvec) .^ (-sigmaH) .* C_si
    CF_gi = (1 .- varrho_vec) .^ sigmaH .* (PV ./ pvec) .^ (-sigmaH) .* C_gi
    CF_si = (1 .- varrho_vec) .^ sigmaH .* (PV ./ pvec) .^ (-sigmaH) .* C_si

    CHi = CH_gi .+ CH_si  # total home consumption by sector
    CFi = CF_gi .+ CF_si  # total import consumption by sector

    # Exports
    PX = prod(pHvec .^ chiX_vec)
    X  = omegaX * (PX / q)^(-etastar) * Ystar
    Xi = chiX_vec .* X .* PX ./ pHvec

    # ---- Inner solve for (M, L, Vi, Yi) ----
    # Initial guess: ignore IO intermediate use first
    ig = zeros(nsec)
    for i in 1:nsec, j in 1:nsec
        ig[i] += beta_mat[j, i]
    end
    ig .*= mean(CHi .+ Xi)

    M_init  = max.(MCi ./ PMi, 1e-20) .^ epsY_vec .* alpha_vec .* (CHi .+ Xi .+ ig)
    L_init  = max.(MCi ./ PL,  1e-20) .^ epsY_vec .* (1 .- alpha_vec .- alphaV_vec) .* (CHi .+ Xi .+ ig)
    Vi_init = max.(MCi ./ PV,  1e-20) .^ epsY_vec .* alphaV_vec .* (CHi .+ Xi .+ ig)
    Yi_init = A_vec .* (
        alpha_vec        .^ (1 ./ epsY_vec) .* max.(M_init,  1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ alphaV_vec      .^ (1 ./ epsY_vec) .* max.(Vi_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ (1 .- alphaV_vec .- alpha_vec) .^ (1 ./ epsY_vec) .* max.(L_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
    ) .^ (epsY_vec ./ (epsY_vec .- 1))

    x0_inner = [M_init; L_init; Vi_init; Yi_init]

    sol = nlsolve(
        (F, x) -> steady_ntwsoe_system!(
            F, x, alpha_vec, alphaV_vec, beta_mat,
            MCi, PMi, PL, PV, CHi, Xi,
            epsY_vec, epsM_vec, A_vec, pHvec
        ),
        x0_inner;
        ftol      = 1e-10,
        show_trace = false,
        method    = :trust_region,
    )

    M_sol  = sol.zero[1:nsec]
    L_sol  = sol.zero[nsec+1:2*nsec]
    Vi_sol = sol.zero[2*nsec+1:3*nsec]
    Yi_sol = sol.zero[3*nsec+1:4*nsec]

    # Intermediate use with converged solution (for market clearing)
    intermediate_use_final = zeros(nsec)
    for i in 1:nsec, j in 1:nsec
        intermediate_use_final[i] += beta_mat[j, i] * (PMi[j] / pHvec[i])^epsM_vec[j] * M_sol[j]
    end

    # Macro aggregates
    TB  = PX * X - PV * (sum(CFi) + sum(Vi_sol))
    GDP = C + TB
    N   = (C^(-GAMMA) * w / CHI)^(1 / PSI)

    # Residuals
    residual = zeros(nsec + 3)
    residual[1:nsec]   .= Yi_sol .- CHi .- Xi .- intermediate_use_final   # output market clearing
    residual[nsec+1]    = TB / GDP - tb_target                             # trade balance target
    residual[nsec+2]    = N - sum(L_sol)                                   # labor market clearing
    residual[nsec+3]    = C - (C_g / om_g)^om_g * (C_s / om_s)^om_s     # consumption identity

    return residual
end
