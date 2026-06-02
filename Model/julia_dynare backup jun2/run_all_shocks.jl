#!/usr/bin/env julia
"""
run_all_shocks.jl
=================
Runs the four shock-analysis scripts — oil price, agriculture / mining /
manufacturing TFP — one after another, each in its OWN Julia subprocess.

Running each script as a separate process (rather than `include`-ing them) keeps
their global state and Dynare model cache isolated, exactly as if you launched
them by hand, but with a single command and a summary at the end.

Usage:
  julia --project=. run_all_shocks.jl            # run all four
  julia --project=. run_all_shocks.jl oil mfg    # run only the named ones

Tags: oil, agr, min, mfg
"""

using Printf

const SCRIPT_DIR = @__DIR__

# (tag, label, script file) — oil first so its IRF CSV exists for anything downstream.
const ALL_SHOCKS = [
    ("oil", "Petróleo",     "oil_shock_analysis.jl"),
    ("agr", "Agricultura",  "agr_shock_analysis.jl"),
    ("min", "Minería",      "min_shock_analysis.jl"),
    ("mfg", "Manufactura",  "mfg_shock_analysis.jl"),
]

# Optional CLI filter: keep only the shocks whose tag was passed as an argument.
selected = isempty(ARGS) ? ALL_SHOCKS : filter(s -> s[1] in ARGS, ALL_SHOCKS)
if isempty(selected)
    println("No matching shock tags in $(ARGS). Valid tags: ", join(first.(ALL_SHOCKS), ", "))
    exit(1)
end

julia_exe   = joinpath(Sys.BINDIR, "julia")
project_dir = dirname(Base.active_project())

@printf "\n%s\n  NK-IOSOE — Ejecutando %d análisis de shock en serie\n%s\n" repeat("=",70) length(selected) repeat("=",70)

results = Tuple{String,Bool,Float64}[]
for (tag, label, script) in selected
    path = joinpath(SCRIPT_DIR, script)
    @printf "\n%s\n  ▶ Shock: %s  (%s)\n%s\n" repeat("─",70) label script repeat("─",70)
    t0 = time()
    proc = run(ignorestatus(`$julia_exe --project=$project_dir $path`); wait=true)
    dt = time() - t0
    ok = success(proc)
    push!(results, (label, ok, dt))
    @printf "\n  %s  —  %s  (%.1f s)\n" (ok ? "✓ OK" : "✗ FALLÓ") label dt
end

@printf "\n%s\n  RESUMEN\n%s\n" repeat("=",70) repeat("=",70)
for (label, ok, dt) in results
    @printf "  %-14s %s   (%.1f s)\n" label (ok ? "✓" : "✗") dt
end
@printf "%s\n" repeat("=",70)

nfail = count(r -> !r[2], results)
if nfail == 0
    @printf "  Todos los shocks se ejecutaron correctamente.\n\n"
else
    @printf "  %d de %d shocks fallaron — revisar la salida arriba.\n\n" nfail length(selected)
end
exit(nfail == 0 ? 0 : 1)
