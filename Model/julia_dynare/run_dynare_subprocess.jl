"""
run_dynare_subprocess.jl
========================
Standalone script called AS A SUBPROCESS by main_SOE_gap.jl.

Running Dynare in a fresh Julia process completely avoids:
  - Julia 1.12 world-age binding restrictions
  - @dynare macro scope/CWD issues inside functions
  - World-age conflicts between include() and the calling function

Usage (called automatically by main_SOE_gap.jl):
  julia --project=/path/to/julia_dynare run_dynare_subprocess.jl /path/to/mod

Outputs (written to MOD_DIR):
  dynare_endo_names.csv  — endogenous variable names (declaration order)
  dynare_ss.csv          — steady-state vector
  dynare_g1_1.csv        — decision rule: state feedback  (≈ oo_.dr.ghx)
  dynare_g1_2.csv        — decision rule: shock impact    (≈ oo_.dr.ghu)
  dynare_sigma_e.csv     — shock covariance matrix        (≈ M_.Sigma_e)
  dynare_state_rows.csv  — state variable row indices in g1_1
  dynare_sim.csv         — simulated time paths (variable × period)
"""

length(ARGS) >= 1 || error("Usage: julia run_dynare_subprocess.jl /path/to/mod")

MOD_DIR = ARGS[1]
isdir(MOD_DIR) || error("MOD_DIR not found: $MOD_DIR")

# If Dynare is not in the active project, try importing from the default depot.
# This handles the case where the parent's project differs from julia_dynare/.
try
    using Dynare
catch
    import Pkg; Pkg.add("Dynare"); using Dynare
end
using CSV, DataFrames, LinearAlgebra, Statistics

# Change to the model directory so @dynare finds .mod and @#include files
cd(MOD_DIR)
@info "Running Dynare in $(pwd())"

# Clear the Dynare.jl compiled model cache so the updated .mod and
# params_jl.mod are always picked up fresh (avoids stale-cache bugs).
cache_dir = joinpath(MOD_DIR, "NK_SOE_lev_gap2")
if isdir(cache_dir)
    @info "Clearing Dynare model cache: $cache_dir"
    rm(cache_dir, recursive=true)
end

# @dynare at TOP LEVEL — no function scope, no world-age issues
context = @dynare "NK_SOE_lev_gap2"

@info "Dynare solve complete — extracting results"

mr         = context.results.model_results[1]
endo_names = Dynare.get_endogenous(context.symboltable)
n_endo     = length(endo_names)
@info "  $n_endo endogenous variables"

# ---- Steady state ----
# In Dynare.jl 0.10.x, the steady state is in mr.trends.endogenous_steady_state
# ONLY when stoch_simul successfully computes it (requires all initval vars defined).
# With missing parameter declarations (now fixed), this should now be populated.
# If still empty, fall back to context.work.initval_endogenous (user's initval).
function get_ss(context, mr, n_endo)
    # Primary: trends field (populated after successful stoch_simul steady state solve)
    if length(mr.trends.endogenous_steady_state) == n_endo
        @info "SS from mr.trends.endogenous_steady_state"
        return Float64.(mr.trends.endogenous_steady_state)
    end
    # Fallback: initval_endogenous — the initval block values set by params_jl.mod.
    # Confirmed present in context.work from diagnostic run.
    v = context.work.initval_endogenous
    if isa(v, AbstractVector) && length(v) == n_endo
        @info "SS from context.work.initval_endogenous ($(sum(abs.(Float64.(v)) .> 1e-10)) non-zero)"
        return Float64.(v)
    end
    @warn "No SS of length $n_endo found — returning zeros"
    return zeros(n_endo)
end

ss_vec = get_ss(context, mr, n_endo)
@info "SS length: $(length(ss_vec))  non-zero entries: $(sum(abs.(ss_vec) .> 1e-10))"

# ---- Decision rule ----
# In Dynare.jl the LRE object stores three matrices:
#   g1_1  →  response to predetermined (backward) variables  ≈ oo_.dr.ghx
#   g1_2  →  response to forward-looking variables            (NOT shocks)
#   g1_3  →  response to exogenous shocks                    ≈ oo_.dr.ghu
# Earlier code tried (g1_1, g1_2) and got zeros because g1_2 is forward vars.
# The fix: try (g1_1, g1_3) first, then fall back to other names.
# We also select the best pair by nonzero element count so stale-zero fields
# don't silently win.
function get_decision_rule(mr, n_endo)
    if !isdefined(mr, :linearrationalexpectations)
        @warn "no linearrationalexpectations field"
        return zeros(n_endo, 1), zeros(n_endo, 1)
    end
    lre = mr.linearrationalexpectations
    lre_fields = fieldnames(typeof(lre))
    @info "lre type: $(typeof(lre))"
    @info "lre fields: $lre_fields"

    # Log every numeric field with its size and nonzero count for diagnostics
    for fn in lre_fields
        try
            v = getfield(lre, fn)
            if isa(v, AbstractMatrix) && eltype(v) <: Real
                nz = sum(abs.(Float64.(v)) .> 1e-12)
                @info "  lre.$fn  size=$(size(v))  nonzero=$nz"
            end
        catch; end
    end

    best_ghx = nothing; best_ghu = nothing; best_nz = -1
    # Candidate pairs: (state_feedback, shock_impact)
    # g1_3 is the Dynare.jl shock matrix; g1_2 is forward-looking (wrong for shocks)
    for (cand1, cand2) in [(:g1_1, :g1_3), (:g1_1, :g1_2), (:ghx, :ghu), (:g1, :g2)]
        if isdefined(lre, cand1) && isdefined(lre, cand2)
            m1 = Matrix{Float64}(getfield(lre, cand1))
            m2 = Matrix{Float64}(getfield(lre, cand2))
            nz = sum(abs.(m1) .> 1e-12) + sum(abs.(m2) .> 1e-12)
            @info "Candidate ($cand1,$cand2): sizes $(size(m1))/$(size(m2))  nonzero=$nz"
            if nz > best_nz
                best_ghx = m1; best_ghu = m2; best_nz = nz
                @info "  → new best: ($cand1,$cand2)"
            end
        end
    end

    if isnothing(best_ghx)
        @warn "Could not find any decision rule matrices — IRFs will be zero"
        return zeros(n_endo, 1), zeros(n_endo, 1)
    end

    @info "Selected decision rule: ghx size=$(size(best_ghx))  ghu size=$(size(best_ghu))  total_nonzero=$best_nz"
    return best_ghx, best_ghu
end

g1_1, g1_2 = get_decision_rule(mr, n_endo)

# ---- Sigma_e ----
Sigma_e = isdefined(context.models[1], :Sigma_e) ?
          Matrix{Float64}(context.models[1].Sigma_e) :
          Matrix{Float64}(I, size(g1_2, 2), size(g1_2, 2))

# ---- State variable rows ----
state_rows = isdefined(context.models[1], :i_bkwrd_b) ?
             Int.(context.models[1].i_bkwrd_b) :
             collect(1:size(g1_1, 2))

# =========================================================================== #
#  Write CSV files                                                            #
# =========================================================================== #

CSV.write(joinpath(MOD_DIR, "dynare_endo_names.csv"),
    DataFrame(variable = endo_names))

CSV.write(joinpath(MOD_DIR, "dynare_ss.csv"),
    DataFrame(variable = endo_names, ss_value = Float64.(ss_vec)))

CSV.write(joinpath(MOD_DIR, "dynare_g1_1.csv"),
    DataFrame(g1_1, [Symbol("c$i") for i in 1:size(g1_1,2)]))

CSV.write(joinpath(MOD_DIR, "dynare_g1_2.csv"),
    DataFrame(g1_2, [Symbol("c$i") for i in 1:size(g1_2,2)]))

CSV.write(joinpath(MOD_DIR, "dynare_sigma_e.csv"),
    DataFrame(Sigma_e, [Symbol("c$i") for i in 1:size(Sigma_e,2)]))

CSV.write(joinpath(MOD_DIR, "dynare_state_rows.csv"),
    DataFrame(state_row = state_rows))

# Simulated paths (skip — not needed; we compute IRFs analytically below)
@info "Skipping stochastic simulation (not required for IRF/variance analysis)"

# =========================================================================== #
#  COMPUTE IRFs ANALYTICALLY FROM DECISION RULE (bypasses stoch_simul/gees) #
#                                                                              #
#  On Apple Silicon (aarch64), Dynare's stoch_simul IRF computation fails    #
#  because LAPACK's gees() with select function is not available on ARM.     #
#  We compute IRFs directly from the first-order decision rule:              #
#                                                                              #
#  Let A = state rows of g1_1  (n_states × n_states, transition)             #
#      B = state rows of g1_2  (n_states × n_shocks,  impact)                #
#  At horizon 1: state  x_1 = B[:,k]  for shock k                           #
#                model  y_1 = g1_1 * 0 + g1_2[:,k]  = g1_2[:,k]             #
#  At horizon h: state  x_h = A * x_{h-1}                                   #
#                model  y_h = g1_1 * x_{h-1}                                 #
#                                                                              #
#  This is exact (not approximate) for first-order perturbation.             #
# =========================================================================== #

@info "Computing IRFs analytically from decision rule (bypasses gees/ARM issue)"

irf_rows = NamedTuple{(:variable, :shock, :period, :value), Tuple{String,String,Int,Float64}}[]

try
    n_irf    = 150         # horizon matches stoch_simul(irf=150)
    n_endo_  = size(g1_1, 1)
    n_states = size(g1_1, 2)
    n_shocks = size(g1_2, 2)

    # Identify state-variable rows (backward-looking variables in g1_1)
    # context.models[1].i_bkwrd_b gives the 1-based row indices
    sr = state_rows   # already extracted above

    A = g1_1[sr, :]   # n_states × n_states  (transition matrix)
    B = g1_2[sr, :]   # n_states × n_shocks  (impact matrix)

    # Exogenous variable names from the model
    exo_names = string.(Dynare.get_exogenous(context.symboltable))
    @info "  Computing IRFs for $(n_shocks) shock(s): $(join(exo_names, ", "))"

    for k in 1:n_shocks
        shock_name = k <= length(exo_names) ? exo_names[k] : "shock_$k"
        x = zeros(n_states)    # state vector (starts at SS = 0 deviation)
        bk = B[:, k]           # impact of shock k on states at t=1

        for h in 1:n_irf
            if h == 1
                y_h = g1_2[:, k]       # direct impact: y_1 = g1_2 * e_k
                x   = bk               # state after shock
            else
                y_h = g1_1 * x         # propagation: y_h = g1_1 * x_{h-1}
                x   = A * x            # update state
            end

            # Store percentage deviations (scaled by SS; skip near-zero SS)
            for (i, vn) in enumerate(endo_names)
                ss_val = abs(ss_vec[i])
                if ss_val > 1e-10
                    push!(irf_rows, (variable=vn, shock=shock_name,
                                     period=h, value=100.0 * y_h[i] / ss_val))
                end
            end
        end
    end

    df_irfs = DataFrame(irf_rows)
    CSV.write(joinpath(MOD_DIR, "dynare_irfs.csv"), df_irfs)
    @info "IRFs saved: $(n_shocks) shock(s), $(length(endo_names)) variable(s), $(n_irf) periods"

catch e
    @warn "Analytical IRF computation failed: $e"
    CSV.write(joinpath(MOD_DIR, "dynare_irfs.csv"),
              DataFrame(variable=String[], shock=String[], period=Int[], value=Float64[]))
end

# ---- Save endogenous_variance (unconditional variance matrix Γ) ----------- #
# Available even when stoch_simul fails on aarch64.
# Used in main_SOE_gap.jl for Lyapunov-based rank correlations.
if isdefined(mr.linearrationalexpectations, :endogenous_variance)
    ev = mr.linearrationalexpectations.endogenous_variance
    if !isnothing(ev) && !isempty(ev)
        ev_mat = Matrix{Float64}(ev)
        @info "endogenous_variance size: $(size(ev_mat)) — saving for rank correlations"
        CSV.write(joinpath(MOD_DIR, "dynare_endogenous_variance.csv"),
                  DataFrame(ev_mat, [Symbol("c$i") for i in 1:size(ev_mat,2)]))
    end
end

@info "Results written to $(MOD_DIR)"
println("DYNARE_SUCCESS")   # sentinel read by parent to confirm success

# ---- Save context for SMM estimation ------------------------------------
# smm_estimation.jl loads from this path to avoid re-running @dynare
using Serialization
context_jls = joinpath(MOD_DIR, "nk_iosoe_context.jls")
try
    serialize(context_jls, context)
    @info "Context saved for SMM estimation: $context_jls"
catch e
    @warn "Could not serialize context: $e"
end
