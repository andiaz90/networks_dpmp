"""
steady_ntwsoe_system.jl
=======================
In-place residual function for the sectoral production block in steady state.
Translated from steady_ntwsoe_system.m (MATLAB).

PRODUCTION
----------
Four CES limbs: domestic materials M, imported inputs V, capital K, and labour L
with the residual weight 1 - alpha_i - alpha_Vi - alpha_Ki.

Capital follows Luttini, Pastén & Rubbo (2024). Their asset f combines a fixed
endowment Kbar_f with an investment good I_f,

    K_f = [(1+phi_f) I_f]^{1/(1+phi_f)} Kbar_f                            (4)

and retailer optimisation gives the capital supply curve

    U_f^{phi_f} = R_f Kbar_f / P^I_f ,     U_f = K_f / Kbar_f             (6)

so nu = 1/phi is the elasticity of capital services to the real rental. Two
accounting identities follow: investment expenditure is R_f K_f/(1+phi) (eq. 9)
and the complement R_f K_f phi/(1+phi) is retailer profit rebated to households
(eq. 8). Investment is therefore a genuine block of FINAL DEMAND and enters
market clearing here.

WHY CAPITAL IS IN THE MODEL AT ALL. With constant returns and labour as the only
primary factor, labour must absorb the entirety of value added — Rubbo (2023,
Econometrica) Remark 3. That implies a labour share of VA of 1.000 against 0.416
in the Chilean accounts, and excedente bruto de explotación is 56.8% of Chilean
VA, concentrated in mining (EBE/GO = 0.58) and housing (0.70) where calling it
pure profit is least defensible. With the nest the labour share of VA is 0.432.

CALIBRATION. alpha_K_i = EBE_i/GO_i and the investment bundle chi_I come from
Data/sector_calibration.csv; nu comes from Data/capital_calibration.csv. LPR
never report a phi, so it is recovered from their own eq. 9 identity applied to
Chilean data: capital income sum(alpha_K*GO) = 121,483 against FBCF = 36,458
(Cuadro 20), giving 1/(1+phi) = 0.3001, phi = 2.332, nu = 0.4288.

TWO SOLVER MODES
----------------
Both solve for [M; L; Vi; Yi], and the difference is only what happens to the
endowment:

  * CALIBRATION (Kfix_vec === nothing) — 5*nsec system, K solved jointly under
    the unit normalisation K_i = alpha_Ki * Y_i. That makes the steady-state
    capital cost share equal alpha_Ki exactly: with A = 1 the capital FOC then
    gives R_i = MC_i, so alpha_Ki (MC_i/R_i)^{epsY-1} = alpha_Ki. Run once by
    main_SOE_gap.jl; the solved K becomes the endowment Kbar written to the .mod.

  * ESTIMATION (Kfix_vec given) — 4*nsec system with the endowment held FIXED at
    its calibrated value. The endowment is a datum, not a function of the
    estimated parameters, so the capital cost share is free to drift as the
    steady state moves — as in Baqaee-Farhi and in LPR, where Kbar_f is fixed.

Economics of the residuals:
  - F1: Material input demand (CES)
  - F2: Labor demand (CES)
  - F3: Imported input demand (CES)
  - F4: CES production function
  - F5: capital endowment normalisation (calibration mode only)
"""

function steady_ntwsoe_system!(
    F        :: AbstractVector{<:Real},
    x        :: AbstractVector{<:Real},
    alpha_vec  :: Vector{<:Real},   # material input shares
    alphaV_vec :: Vector{<:Real},   # import input shares
    beta_mat   :: Matrix{<:Real},   # IO matrix: beta_mat[i,j] = share of sector j inputs from sector i
    MCi        :: Vector{<:Real},   # nominal marginal costs
    PMi        :: Vector{<:Real},   # intermediate input price indices
    PL         :: Vector{<:Real},   # sectoral labor prices (= w in SS)
    PV         :: Real,             # import price (= Q * PVstar)
    CHi        :: Vector{<:Real},   # home consumption demand by sector
    Xi         :: Vector{<:Real},   # export demand by sector
    epsY_vec   :: Vector{<:Real},   # production input substitution elasticities
    epsM_vec   :: Vector{<:Real},   # materials substitution elasticities
    A_vec      :: Vector{<:Real},   # TFP by sector
    pHvec      :: Vector{<:Real},   # home prices by sector
    alphaK_vec :: Vector{<:Real},   # capital shares, EBE_i/GO_i
    chiI_vec   :: Vector{<:Real},   # investment bundle weights (FBCF, sum to 1)
    nuK        :: Real,             # capital supply elasticity = 1/phi
    Kfix_vec   :: Union{Nothing,Vector{<:Real}} = nothing,  # endowment (estimation mode)
)
    nsec = length(alpha_vec)
    M   = x[1:nsec]
    L   = x[nsec+1:2*nsec]
    Vi  = x[2*nsec+1:3*nsec]
    Yi  = x[3*nsec+1:4*nsec]

    solve_K = Kfix_vec === nothing
    if solve_K && length(x) < 5 * nsec
        error("steady_ntwsoe_system!: calibration mode needs x of length " *
              "$(5*nsec), got $(length(x)). Extend the initial guess.")
    end
    K = solve_K ? x[4*nsec+1:5*nsec] : convert(Vector{eltype(x)}, Kfix_vec)

    # Intermediate use: sector i's output used as input in sector j
    # intermediate_use[i] = Σ_j beta_mat[j,i] * (PMi[j]/pHvec[i])^epsM_vec[j] * M[j]
    intermediate_use = zeros(eltype(x), nsec)
    for i in 1:nsec
        for j in 1:nsec
            intermediate_use[i] += beta_mat[j, i] * (PMi[j] / pHvec[i])^epsM_vec[j] * M[j]
        end
    end

    # Investment demand. R_i comes from the capital FOC rather than from R = MC,
    # because R = MC only holds AT the solution (where K = alpha_K Y), not during
    # the iteration.
    Kpos = max.(K, 1e-20)
    R_i = A_vec .^ ((epsY_vec .- 1) ./ epsY_vec) .* MCi .*
          (alphaK_vec .* max.(Yi, 1e-20) ./ Kpos) .^ (1 ./ epsY_vec)
    EInv = (nuK / (1 + nuK)) * sum(R_i .* Kpos)
    inv_demand = chiI_vec .* EInv ./ pHvec

    demand = CHi .+ Xi .+ intermediate_use .+ inv_demand

    # Clamp price ratios to avoid complex exponentiation with fractional epsY
    rM  = max.(MCi ./ PMi, 1e-20)
    rL  = max.(MCi ./ PL,  1e-20)
    rV  = max.(MCi ./ PV,  1e-20)

    # F1: material input demand
    F[1:nsec] .= M .- rM .^ epsY_vec .* alpha_vec .* demand

    # F2: labor demand — weight net of the capital share
    F[nsec+1:2*nsec] .= L .- rL .^ epsY_vec .* (1 .- alpha_vec .- alphaV_vec .- alphaK_vec) .* demand

    # F3: imported input demand
    F[2*nsec+1:3*nsec] .= Vi .- rV .^ epsY_vec .* alphaV_vec .* demand

    # F4: CES production function (clamp negatives for robustness)
    Mpos  = max.(M,  1e-20)
    Vipos = max.(Vi, 1e-20)
    Lpos  = max.(L,  1e-20)

    F[3*nsec+1:4*nsec] .= Yi .- A_vec .* (
         alpha_vec        .^ (1 ./ epsY_vec) .* Mpos  .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ alphaV_vec       .^ (1 ./ epsY_vec) .* Vipos .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ alphaK_vec       .^ (1 ./ epsY_vec) .* Kpos  .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ (1 .- alphaV_vec .- alpha_vec .- alphaK_vec) .^ (1 ./ epsY_vec) .* Lpos .^ ((epsY_vec .- 1) ./ epsY_vec)
    ) .^ (epsY_vec ./ (epsY_vec .- 1))

    # F5: capital endowment normalisation, K_i = alpha_Ki * Y_i.
    # A CHOICE OF UNITS for the fixed factor, not a behavioural equation.
    if solve_K
        F[4*nsec+1:5*nsec] .= K .- alphaK_vec .* Yi
    end

    return F
end
