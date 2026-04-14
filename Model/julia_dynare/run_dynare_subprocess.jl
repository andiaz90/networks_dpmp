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

@info "Dynare solve complete — extracting results"

# =========================================================================== #
#  Extract results from context                                               #
# =========================================================================== #

mr         = context.results.model_results[1]
endo_names = Dynare.get_endogenous(context.symboltable)
ss_vec     = mr.trends.endogenous_steady_state
lre        = mr.linearrationalexpectations
g1_1       = Matrix{Float64}(lre.g1_1)   # n_endo × n_states
g1_2       = Matrix{Float64}(lre.g1_2)   # n_endo × n_shocks
Sigma_e    = Matrix{Float64}(context.models[1].Sigma_e)
state_rows = Int.(context.models[1].i_bkwrd_b)

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
