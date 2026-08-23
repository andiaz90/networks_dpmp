"""
steady_ntwsoe.jl
================
Outer steady-state residuals for the NK-IOSOE model.
Translated from steady_ntwsoe.m (MATLAB).

Unknowns: x_vec = [pH_1,...,pH_nsec, w, Q, C, Gscale]  (length nsec+4)

MARKUP SUBSIDY (LPR 2024 eq. 13, added 2026-08-20). Their subsidy
1 - tau*_i = (eps_i - 1)/eps_i makes price equal marginal cost in steady state,
so revenue equals total cost and each sector's factor COST shares become its
factor REVENUE shares. Without it the model's labour share of value added was
0.336 and capital 0.490, summing to 0.826 — the missing 17.4% being exactly the
eps = 10 monopoly rent, 1/eps of gross output. Callers pass mc_over_ph = 1 under
the subsidy and (eps-1)/eps without it.

The function:
  1. Derives all prices and demands from (pH, w, Q, C)
  2. Calls steady_ntwsoe_system! (inner NLsolve) for (M, L, Vi, Yi, K)
  3. Returns nsec+3 market-clearing residuals

GOVERNMENT (2026-08-20): G_vec is real government consumption by sector, held
EXOGENOUS. It enters goods-market clearing and GDP but not the household CES —
27% of Chilean final consumption is government, concentrated in administración
pública (98.7% of its final demand) and servicios personales (53.4%), and making
it price-elastic would distort exactly the object this model is for. Leaving it
out made administración pública 16x too small in steady state: 0.23% of gross
output against 3.63% in the IO table. G is added to CHi before the inner solve —
the inner system uses CHi only inside `demand`, so it needs no signature change.

Its LEVEL is an unknown of this system, pinned by the extra residual
Gnom/GDP = gshare_target, because GDP is only available here. Solving it outside
would require targeting G/C instead, which is a different (and wrong) ratio —
C is about two thirds of GDP once investment and the trade balance are counted.

Capital (LPR 2024, see steady_ntwsoe_system.jl): the inner solve runs in
CALIBRATION mode, so K is solved jointly under K_i = alpha_Ki * Y_i. The
investment expenditure it generates, EInv = nu/(1+nu) * sum(R_i K_i), is real
final demand and enters BOTH the goods-market residual and GDP.

Returns residual vector of length nsec+4:
  residual[1:nsec]   = output market clearing (Yi = CHi + Xi + intermediate_use)
  residual[nsec+1]   = trade balance target   (TB/GDP = tb_target)
  residual[nsec+2]   = labor market clearing  (N = Σ L_i)
  residual[nsec+3]   = Cobb-Douglas consumption identity
  residual[nsec+4]   = government consumption share  (Gnom/GDP = gshare_target)
"""

using NLsolve, Statistics

function steady_ntwsoe(
    x_vec    :: AbstractVector{<:Real},
    PVstar   :: Real,
    mc_over_ph :: Real,             # steady-state MC/PH = (eps-1)/(eps*subsMC)
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
    gammaG_vec :: Vector{<:Real},   # government demand composition, sums to 1
    gshare_target :: Real,          # target G / GDP
    alphaK_vec :: Vector{<:Real},
    chiI_vec   :: Vector{<:Real},
    nuK        :: Real,
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
    Gscale = max(x_vec[nsec+4], 1e-20)
    G_vec  = gammaG_vec .* Gscale

    PL  = fill(w, nsec)
    PV  = q * PVstar

    # Marginal costs (Bertrand pricing: MC = (ε-1)/ε × PH)
    # MC/PH is passed in rather than derived from epsilon, because the LPR
    # markup subsidy changes it: 1 - tau*_i = (eps_i-1)/eps_i sets MC/PH = 1.
    MCi = mc_over_ph .* pHvec

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
    ig .*= mean(CHi .+ G_vec .+ Xi)

    M_init  = max.(MCi ./ PMi, 1e-20) .^ epsY_vec .* alpha_vec .* (CHi .+ G_vec .+ Xi .+ ig)
    L_init  = max.(MCi ./ PL,  1e-20) .^ epsY_vec .* (1 .- alpha_vec .- alphaV_vec .- alphaK_vec) .* (CHi .+ G_vec .+ Xi .+ ig)
    Vi_init = max.(MCi ./ PV,  1e-20) .^ epsY_vec .* alphaV_vec .* (CHi .+ G_vec .+ Xi .+ ig)
    Yi_init = A_vec .* (
        alpha_vec        .^ (1 ./ epsY_vec) .* max.(M_init,  1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ alphaV_vec      .^ (1 ./ epsY_vec) .* max.(Vi_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ (1 .- alphaV_vec .- alpha_vec .- alphaK_vec) .^ (1 ./ epsY_vec) .* max.(L_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
    ) .^ (epsY_vec ./ (epsY_vec .- 1))

    x0_inner = [M_init; L_init; Vi_init; Yi_init; alphaK_vec .* Yi_init]

    # Government demand rides along with household demand into the inner system.
    CHi_G = CHi .+ G_vec

    sol = nlsolve(
        (F, x) -> steady_ntwsoe_system!(
            F, x, alpha_vec, alphaV_vec, beta_mat,
            MCi, PMi, PL, PV, CHi_G, Xi,
            epsY_vec, epsM_vec, A_vec, pHvec,
            alphaK_vec, chiI_vec, nuK
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
    K_sol  = sol.zero[4*nsec+1:5*nsec]

    # Investment demand. R_i from the capital FOC (at the inner solution
    # K = alpha_K Y this equals MC_i, but computing it from the FOC keeps the
    # two in step if the normalisation ever changes).
    R_sol = A_vec .^ ((epsY_vec .- 1) ./ epsY_vec) .* MCi .*
            (alphaK_vec .* max.(Yi_sol, 1e-20) ./ max.(K_sol, 1e-20)) .^ (1 ./ epsY_vec)
    EInv  = (nuK / (1 + nuK)) * sum(R_sol .* max.(K_sol, 1e-20))
    inv_use = chiI_vec .* EInv ./ pHvec

    # Intermediate use with converged solution (for market clearing)
    intermediate_use_final = zeros(nsec)
    for i in 1:nsec, j in 1:nsec
        intermediate_use_final[i] += beta_mat[j, i] * (PMi[j] / pHvec[i])^epsM_vec[j] * M_sol[j]
    end

    # Macro aggregates
    TB  = PX * X - PV * (sum(CFi) + sum(Vi_sol))
    Gnom = sum(pHvec .* G_vec)          # nominal government consumption
    GDP = C + Gnom + EInv + TB
    N   = (C^(-GAMMA) * w / CHI)^(1 / PSI)

    # Residuals
    residual = zeros(nsec + 4)
    residual[1:nsec]   .= Yi_sol .- CHi_G .- Xi .- intermediate_use_final .- inv_use  # output market clearing
    residual[nsec+1]    = TB / GDP - tb_target                             # trade balance target
    residual[nsec+2]    = N - sum(L_sol)                                   # labor market clearing
    residual[nsec+3]    = C - (C_g / om_g)^om_g * (C_s / om_s)^om_s     # consumption identity
    residual[nsec+4]    = Gnom / GDP - gshare_target                       # government share

    return residual
end
