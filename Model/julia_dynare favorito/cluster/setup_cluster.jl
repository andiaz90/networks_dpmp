"""
setup_cluster.jl
================
One-time package installation for the NK-IOSOE SMM estimation.

Run this ONCE on the cluster login node (NOT as a SLURM job) before
submitting run_smm_hpc.sh:

    julia --project=. cluster/setup_cluster.jl

This installs all packages into the depot specified by JULIA_DEPOT_PATH
(or the default ~/.julia if the variable is not set).

Expected runtime: 5-15 minutes on first run (downloads + precompiles).
Subsequent runs are instant (packages already cached).
"""

import Pkg

println("="^60)
println("  NK-IOSOE SMM — Cluster Package Setup")
println("="^60)
println()
println("Julia version: $(VERSION)")
println("Project:       $(Pkg.project().path)")
println("Depot:         $(first(Base.DEPOT_PATH))")
println()

# Instantiate: installs all packages listed in Project.toml + Manifest.toml.
# If Manifest.toml is present (committed from a local machine), this pins
# exact versions. If absent, resolves from Project.toml constraints.
println("--- Step 1: Instantiating packages ---")
Pkg.instantiate()
println("  Done.\n")

# Precompile all packages so the first SLURM job doesn't pay compile time.
println("--- Step 2: Precompiling packages ---")
Pkg.precompile()
println("  Done.\n")

# Quick import test — catches broken installs before job submission.
println("--- Step 3: Import test ---")
import_tests = [
    "using LinearAlgebra, Statistics, Printf",
    "using CSV, DataFrames",
    "using NLsolve",
    "using StatsBase",
    "using CMAEvolutionStrategy",
    "using Dynare",
]
for t in import_tests
    try
        eval(Meta.parse(t))
        println("  OK: $t")
    catch e
        println("  FAIL: $t")
        println("       $(e)")
    end
end

println()
println("="^60)
println("  Setup complete. You can now submit:")
println("    sbatch cluster/run_smm_hpc.sh")
println("="^60)
