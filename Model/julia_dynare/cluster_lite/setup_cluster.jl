#=
setup_cluster.jl
================
One-time package installation for the NK-IOSOE pipeline (estimation + shocks).

Run ONCE on the cluster login node before submitting any SLURM job:

    julia --project=. cluster/setup_cluster.jl

Installs + precompiles every package needed for BOTH the SMM estimation
and the shock-plot generation (Plots / StatsPlots included).

Expected runtime: 8-20 min on first run (downloads + precompiles, incl. the
Dynare.jl preprocessor binary). Subsequent runs are instant.
=#

import Pkg

println("="^62)
println("  NK-IOSOE — Cluster package setup (estimation + shocks)")
println("="^62)
println("Julia version: $(VERSION)")
println("Project:       $(Pkg.project().path)")
println("Depot:         $(first(Base.DEPOT_PATH))")
println()

println("--- Step 1: Instantiate ---")
Pkg.instantiate()
println("  Done.\n")

println("--- Step 2: Precompile ---")
Pkg.precompile()
println("  Done.\n")

println("--- Step 3: Import test ---")
import_tests = [
    "using LinearAlgebra, Statistics, Printf, Logging, Serialization, SparseArrays",
    "using CSV, DataFrames, XLSX",
    "using NLsolve",
    "using Random",
    "using StatsBase",
    "using CMAEvolutionStrategy",
    "using Plots, StatsPlots",
    "using Dynare",
]
allok = true
for t in import_tests
    try
        eval(Meta.parse(t)); println("  OK:   $t")
    catch e
        allok = false; println("  FAIL: $t\n        $(e)")
    end
end

println()
println("="^62)
if allok
    println("  Setup complete. Next:")
    println("    sbatch cluster/run_full_pipeline.sh     # moments+compile+estimate+shocks")
    println("  or step by step:")
    println("    sbatch cluster/run_estimation.sh")
    println("    sbatch cluster/run_shocks.sh            # after estimation finishes")
else
    println("  Setup FAILED for one or more packages — fix before submitting.")
end
println("="^62)
