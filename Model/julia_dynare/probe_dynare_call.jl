"""
probe_dynare_call.jl
====================
Dynare's compute_first_order_solution! throws `AssertionError: length(residual) == 6`.
591 equations, so 6 comes from somewhere unexpected in the deserialized context.
This reproduces the call with a FULL backtrace and dumps the sizes of everything
being passed, so the offending object can be identified in one run.

    julia --project=. probe_dynare_call.jl
"""

using CSV, DataFrames, Printf, Serialization, Dynare, LinearAlgebra

SCRIPT_DIR = @__DIR__
context = deserialize(joinpath(SCRIPT_DIR, "mod", "nk_iosoe_context.jls"))
m   = context.models[1]
res = context.results.model_results[1]

@printf "\n== dimensions ==\n"
@printf "  endogenous_nbr        %d\n" m.endogenous_nbr
@printf "  exogenous_nbr         %d\n" m.exogenous_nbr
@printf "  steady state length   %d\n" length(res.trends.endogenous_steady_state)
@printf "  work.params length    %d\n" length(context.work.params)
for f in (:residuals, :dynamic_variables, :exogenous_variables, :jacobian)
    v = try getfield(context.work, f) catch; nothing end
    @printf "  work.%-20s %s\n" String(f) (v === nothing ? "unreadable" :
        (v isa AbstractArray ? string(size(v)) : string(typeof(v))))
end
for f in (:dynamic_tmp_nbr, :static_tmp_nbr, :orig_endo_nbr, :n_bkwrd, :n_fwrd, :n_both, :n_static)
    v = try getfield(m, f) catch; nothing end
    v === nothing || @printf "  model.%-19s %s\n" String(f) string(v)
end
@printf "  workspaces            %s\n" string(try keys(context.workspaces) catch; "unreadable" end)

@printf "\n== constructing workspaces ==\n"
ws     = Dynare.DynamicWs(context; order = 1);         @printf "  DynamicWs           ok\n"
css_ws = Dynare.ComputeStochSimulWs(context);          @printf "  ComputeStochSimulWs ok\n"
opts   = Dynare.StochSimulOptions(Dict{String,Any}("order" => 1, "irf" => 0))
@printf "  StochSimulOptions   ok\n"
for f in (:residuals, :temporary_values, :dynamic_variables)
    v = try getfield(ws, f) catch; nothing end
    v === nothing || @printf "  ws.%-20s %s\n" String(f) (v isa AbstractArray ? string(size(v)) : string(typeof(v)))
end

endo_ss = res.trends.endogenous_steady_state
exo_ss  = try res.trends.exogenous_steady_state catch; zeros(Float64, m.exogenous_nbr) end
@printf "\n== calling compute_first_order_solution! ==\n"
@printf "  endo_ss %d  exo_ss %d  params %d\n" length(endo_ss) length(exo_ss) length(context.work.params)
try
    Dynare.compute_first_order_solution!(
        context, endo_ss, exo_ss, endo_ss, context.work.params, m,
        ws, css_ws, opts; variance_decomposition = false)
    @printf "  SUCCEEDED\n"
catch err
    @printf "\n---- FULL BACKTRACE ----\n"
    showerror(stdout, err, catch_backtrace())
    println()
end
