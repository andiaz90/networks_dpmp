"""
steady_ntwsoe_system.jl
=======================
In-place residual function for the sectoral production block in steady state.
Translated from steady_ntwsoe_system.m (MATLAB).

Solves for [M; L; Vi; Yi] (4×nsec unknowns) given prices and demands.

Economics:
  - F1: Material input demand (CES)
  - F2: Labor demand (CES)
  - F3: Imported input demand (CES)
  - F4: CES production function
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
)
    nsec = length(alpha_vec)
    M   = x[1:nsec]
    L   = x[nsec+1:2*nsec]
    Vi  = x[2*nsec+1:3*nsec]
    Yi  = x[3*nsec+1:4*nsec]

    # Intermediate use: sector i's output used as input in sector j
    # intermediate_use[i] = Σ_j beta_mat[j,i] * (PMi[j]/pHvec[i])^epsM_vec[j] * M[j]
    intermediate_use = zeros(nsec)
    for i in 1:nsec
        for j in 1:nsec
            intermediate_use[i] += beta_mat[j, i] * (PMi[j] / pHvec[i])^epsM_vec[j] * M[j]
        end
    end

    demand = CHi .+ Xi .+ intermediate_use

    # Clamp price ratios to avoid complex exponentiation with fractional epsY
    rM  = max.(MCi ./ PMi, 1e-20)
    rL  = max.(MCi ./ PL,  1e-20)
    rV  = max.(MCi ./ PV,  1e-20)

    # F1: material input demand
    F[1:nsec] .= M .- rM .^ epsY_vec .* alpha_vec .* demand

    # F2: labor demand
    F[nsec+1:2*nsec] .= L .- rL .^ epsY_vec .* (1 .- alpha_vec .- alphaV_vec) .* demand

    # F3: imported input demand
    F[2*nsec+1:3*nsec] .= Vi .- rV .^ epsY_vec .* alphaV_vec .* demand

    # F4: CES production function (clamp negatives for robustness)
    Mpos  = max.(M,  1e-20)
    Vipos = max.(Vi, 1e-20)
    Lpos  = max.(L,  1e-20)

    F[3*nsec+1:4*nsec] .= Yi .- A_vec .* (
         alpha_vec        .^ (1 ./ epsY_vec) .* Mpos  .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ alphaV_vec       .^ (1 ./ epsY_vec) .* Vipos .^ ((epsY_vec .- 1) ./ epsY_vec)
      .+ (1 .- alphaV_vec .- alpha_vec) .^ (1 ./ epsY_vec) .* Lpos .^ ((epsY_vec .- 1) ./ epsY_vec)
    ) .^ (epsY_vec ./ (epsY_vec .- 1))

    return F
end
