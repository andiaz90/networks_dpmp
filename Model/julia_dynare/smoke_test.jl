"""
smoke_test.jl
=============
Cheap, fast pre-flight check — run it BEFORE submitting a multi-hour calibration
so a broken build / bad edit / missing input is caught in seconds, not on day 2.

It verifies, in order:
  1. core Julia packages load,
  2. every SMM source file PARSES and includes without error
     (this is what catches a syntax mistake in an edit to smm_estimation.jl etc.),
  3. the setup contract holds (params/moments/bounds dimensions in sync),
  4. the required raw input files are present.

Exits 0 on success, 1 on any failure, so the SLURM script can abort the run:

    julia --project=. smoke_test.jl
"""

using Printf

# Smoke mode: lets smm_estimation.jl include with placeholder data moments on a
# fresh bundle (CSVs are built later in the pipeline). Real estimation runs do
# NOT set this and hard-error on missing/stale moment files (2026-07-08).
ENV["SMM_SMOKE"] = "1"

const JD   = @__DIR__
const REPO = abspath(joinpath(JD, "..", ".."))
const DATA = joinpath(REPO, "Data")

const FAILS = String[]
function check(name::AbstractString, ok::Bool, msg::AbstractString="")
    if ok
        @printf "  OK    %s\n" name
    else
        push!(FAILS, name)
        @printf "  FAIL  %s  %s\n" name msg
    end
    return ok
end

println("="^60); println("  NK-IOSOE smoke test"); println("="^60)
@printf "Julia %s | %d thread(s) | dir %s\n\n" string(VERSION) Threads.nthreads() JD

# --- 1. packages ----------------------------------------------------------- #
println("--- packages ---")
for p in ["LinearAlgebra","Statistics","StatsBase","Printf","Random","NLsolve",
          "CSV","DataFrames","SparseArrays","CMAEvolutionStrategy","Dynare"]
    try
        @eval using $(Symbol(p))
        check(p, true)
    catch e
        check(p, false, sprint(showerror, e))
    end
end

# --- 2. source files parse + include (syntax check) ------------------------ #
# smm_model_moments.jl already includes the steady-state + utils files; including
# smm_estimation.jl runs its top-level data-moment load, which safely falls back
# to placeholders when the moment CSVs do not exist yet (pre-compile stage).
println("\n--- source includes (syntax check) ---")
cd(JD) do
    for s in ["steady_ntwsoe_system.jl","steady_ntwsoe.jl","utils.jl",
              "smm_model_moments.jl","smm_estimation.jl"]
        try
            Base.include(Main, joinpath(JD, s))
            check(s, true)
        catch e
            check(s, false, sprint(showerror, e))
        end
    end
end

# --- 3. setup contract ----------------------------------------------------- #
println("\n--- setup contract ---")
try
    Base.invokelatest(getfield(Main, :validate_smm_setup), getfield(Main, :data_moments))
    check("validate_smm_setup", true)
catch e
    check("validate_smm_setup", false, sprint(showerror, e))
end

# --- 4. required input files ----------------------------------------------- #
println("\n--- input files ---")
for f in ["sector_calibration.csv", "IO_2021_chile.csv"]
    check("Data/$f", isfile(joinpath(DATA, f)))
end

println("\n" * "="^60)
if isempty(FAILS)
    println("  SMOKE TEST PASSED"); exit(0)
else
    @printf "  SMOKE TEST FAILED (%d): %s\n" length(FAILS) join(FAILS, ", ")
    exit(1)
end
