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
function get_decision_rule(mr, n_endo)
    if !isdefined(mr, :linearrationalexpectations)
        @warn "no linearrationalexpectations field"
        return zeros(n_endo, 1), zeros(n_endo, 1)
    end
    lre = mr.linearrationalexpectations
    @info "lre fields: $(fieldnames(typeof(lre)))"

    g1_1 = nothing; g1_2 = nothing
    for (cand1, cand2) in [(:g1_1, :g1_2), (:ghx, :ghu), (:g1, :g2)]
        if isdefined(lre, cand1) && isdefined(lre, cand2)
            g1_1 = Matrix{Float64}(getfield(lre, cand1))
            g1_2 = Matrix{Float64}(getfield(lre, cand2))
            @info "Found decision rule as lre.$cand1 / lre.$cand2 — sizes: $(size(g1_1)) / $(size(g1_2))"
            break
        end
    end
    if isnothing(g1_1)
        @warn "Could not find decision rule matrices"
        g1_1 = zeros(n_endo, 1); g1_2 = zeros(n_endo, 1)
    end
    return g1_1, g1_2
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

# Simulated paths (if available)
if !isempty(mr.simulations)
    sim_data = Matrix{Float64}(mr.simulations[1].data)  # periods × variables
    df_sim   = DataFrame(sim_data, [Symbol(n) for n in endo_names])
    CSV.write(joinpath(MOD_DIR, "dynare_sim.csv"), df_sim)
    @info "Simulation saved ($(size(sim_data,1)) periods × $(size(sim_data,2)) variables)"
else
    @info "No simulation data"
end

# ---- IRFs (impulse response functions) ----
# Diagnose mr.irfs type first so we can handle any Dynare.jl 0.10.x structure
@info "mr.irfs type: $(typeof(mr.irfs))  isempty: $(isempty(mr.irfs))"

irf_rows = NamedTuple{(:variable, :shock, :period, :value), Tuple{String,String,Int,Float64}}[]

if !isempty(mr.irfs)
    # Determine structure: Dict or AxisArray?
    first_val = first(values(mr.irfs))
    @info "  first shock value type: $(typeof(first_val))"

    try
        for (shock_name, var_data) in mr.irfs
            if var_data isa AbstractDict
                # Dict{varname => Vector}
                for (var_name, irf_vec) in var_data
                    for (t, v) in enumerate(irf_vec)
                        push!(irf_rows, (variable=string(var_name), shock=string(shock_name),
                                         period=t, value=Float64(v)))
                    end
                end
            else
                # AxisArrayTable or Matrix: rows=periods, cols=variables
                mat = Matrix{Float64}(var_data)
                vnames = string.(names(var_data))   # column names = variable names
                for (ci, vn) in enumerate(vnames)
                    for t in 1:size(mat, 1)
                        push!(irf_rows, (variable=vn, shock=string(shock_name),
                                         period=t, value=mat[t, ci]))
                    end
                end
            end
        end
    catch e
        @warn "IRF extraction failed: $e"
    end

    if !isempty(irf_rows)
        df_irfs = DataFrame(irf_rows)
        CSV.write(joinpath(MOD_DIR, "dynare_irfs.csv"), df_irfs)
        @info "IRFs saved: $(length(mr.irfs)) shock(s), $(length(unique(df_irfs.variable))) variable(s)"
    else
        @warn "IRF extraction produced no rows despite non-empty mr.irfs"
        CSV.write(joinpath(MOD_DIR, "dynare_irfs.csv"),
                  DataFrame(variable=String[], shock=String[], period=Int[], value=Float64[]))
    end
else
    @info "No IRF data in mr.irfs — stoch_simul may not have completed"
    CSV.write(joinpath(MOD_DIR, "dynare_irfs.csv"),
              DataFrame(variable=String[], shock=String[], period=Int[], value=Float64[]))
end

# Also save endogenous_variance if available (unconditional variance matrix)
if isdefined(mr.linearrationalexpectations, :endogenous_variance)
    ev = mr.linearrationalexpectations.endogenous_variance
    if !isnothing(ev) && !isempty(ev)
        ev_mat = Matrix{Float64}(ev)
        @info "endogenous_variance size: $(size(ev_mat))"
        CSV.write(joinpath(MOD_DIR, "dynare_endogenous_variance.csv"),
                  DataFrame(ev_mat, [Symbol("c$i") for i in 1:size(ev_mat,2)]))
    end
end

@info "Results written to $(MOD_DIR)"
println("DYNARE_SUCCESS")   # sentinel read by parent to confirm success
