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
using CSV, DataFrames, LinearAlgebra

# Change to the model directory so @dynare finds .mod and @#include files
cd(MOD_DIR)
@info "Running Dynare in $(pwd())"

# @dynare at TOP LEVEL — no function scope, no world-age issues
context = @dynare "NK_SOE_lev_gap2"

@info "Dynare solve complete — inspecting context structure"

# =========================================================================== #
#  Inspect Dynare.jl context to find the correct field names in v0.10.x      #
# =========================================================================== #

mr         = context.results.model_results[1]
endo_names = Dynare.get_endogenous(context.symboltable)
n_endo     = length(endo_names)

# Print context structure for debugging
@info "context type: $(typeof(context))"
@info "context fields: $(fieldnames(typeof(context)))"
@info "mr type: $(typeof(mr))"
@info "mr fields: $(fieldnames(typeof(mr)))"
if isdefined(mr, :trends)
    @info "mr.trends type: $(typeof(mr.trends))"
    @info "mr.trends fields: $(fieldnames(typeof(mr.trends)))"
    @info "mr.trends.endogenous_steady_state length: $(length(mr.trends.endogenous_steady_state))"
end
if isdefined(context, :work)
    @info "context.work fields: $(fieldnames(typeof(context.work)))"
    if isdefined(context.work, :steady_state)
        @info "context.work.steady_state length: $(length(context.work.steady_state))"
    end
end

# ---- Steady state ----
# In Dynare.jl 0.10.x, the steady state is in mr.trends.endogenous_steady_state
# ONLY when stoch_simul successfully computes it (requires all initval vars defined).
# With missing parameter declarations (now fixed), this should now be populated.
# If still empty, fall back to context.work.initval_endogenous (user's initval).
function get_ss(context, mr, n_endo)
    # Primary: trends field (populated after successful stoch_simul solve)
    if isdefined(mr, :trends) && length(mr.trends.endogenous_steady_state) == n_endo
        @info "SS from mr.trends.endogenous_steady_state"
        return Float64.(mr.trends.endogenous_steady_state)
    end
    # Fallback 1: initval_endogenous from context.work
    if isdefined(context, :work) && isdefined(context.work, :initval_endogenous)
        v = context.work.initval_endogenous
        if isa(v, AbstractVector) && length(v) == n_endo
            @info "SS from context.work.initval_endogenous"
            return Float64.(v)
        end
    end
    # Fallback 2: any non-boolean vector field of correct length on mr or lre
    for obj in [mr, isdefined(mr, :linearrationalexpectations) ? mr.linearrationalexpectations : nothing]
        obj === nothing && continue
        for fname in fieldnames(typeof(obj))
            v = getfield(obj, fname)
            if isa(v, AbstractVector) && length(v) == n_endo
                vf = Float64.(v)
                # Skip boolean flags (stationary_variables has range [0,1] with integer values)
                if any(x -> !(x ≈ 0.0 || x ≈ 1.0), vf) || std(vf) > 0.01
                    @info "SS from $(typeof(obj)).$fname"
                    return vf
                end
            end
        end
    end
    @warn "Could not find steady state — returning zeros. Check that pi_ss, r_ss, etc. are declared in parameters block."
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
    @info "No simulation data available"
end

@info "Results written to $(MOD_DIR)"
println("DYNARE_SUCCESS")   # sentinel read by parent to confirm success
