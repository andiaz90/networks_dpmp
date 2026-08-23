"""
plot_chiib_comparison.jl
========================
Overlay oil-shock IRFs for different values of the debt-elastic risk-premium
parameter phi_b (SOE closure sensitivity).

Prerequisite runs (each writes tables/<tag>_shock_irfs.csv):
    julia --project=. oil_shock_analysis.jl                       # phi_b = 0.001 → tag "oil"
    PHIB_OVERRIDE=0.01 julia --project=. oil_shock_analysis.jl   # tag "oil_chiib0p01"
    PHIB_OVERRIDE=0.05 julia --project=. oil_shock_analysis.jl   # tag "oil_chiib0p05"

Then:
    julia --project=. plot_chiib_comparison.jl

Output: figures/oil_shock/irf_chiib_comparison_oil_shock.png
"""

using CSV, DataFrames, Plots, Printf

const JD  = @__DIR__
const TAB = joinpath(JD, "tables")
const OUT = joinpath(JD, "figures", "oil_shock")
isdir(OUT) || mkpath(OUT)

# (tag, phi_b, color)
runs = [
    ("oil",           0.001, :black),
    ("oil_chiib0p01", 0.01,  :steelblue),
    ("oil_chiib0p05", 0.05,  :firebrick),
]

# variable => (title, ylabel)
panels = [
    ("GDP", "PIB",                  "% desv. EE"),
    ("pi",  "Inflación",            "pp"),
    ("Q",   "TCR",                  "% desv. EE"),
    ("TB",  "Balanza comercial",    "% desv. EE"),
    ("C",   "Consumo",              "% desv. EE"),
    ("r",   "Tasa de política",     "pp"),
]

loaded = Tuple{Float64,Symbol,DataFrame}[]
for (tag, cb, col) in runs
    f = joinpath(TAB, "$(tag)_shock_irfs.csv")
    if isfile(f)
        df = CSV.read(f, DataFrame)
        push!(loaded, (cb, col, df))
    else
        @printf "  [skip] %s not found — run oil_shock_analysis.jl with PHIB_OVERRIDE=%g first\n" f cb
    end
end
isempty(loaded) && error("No IRF csv files found in $(TAB).")

plt = plot(layout=(2,3), size=(1350, 750), left_margin=8Plots.mm,
           bottom_margin=6Plots.mm, legend=:outerbottom, legend_columns=length(loaded))

for (k, (vname, ttl, ylab)) in enumerate(panels)
    for (cb, col, df) in loaded
        sub = sort(df[df.variable .== vname, :], :period)
        isempty(sub) && continue
        plot!(plt[k], sub.period, sub.value, lw=2.5, color=col,
              label=(k == 1 ? @sprintf("\\chi_b = %g%s", cb, cb == 0.001 ? " (base)" : "") : ""))
    end
    hline!(plt[k], [0.0], color=:gray, ls=:dash, lw=1, label="")
    plot!(plt[k], title=ttl, ylabel=ylab, xlabel="Trimestres")
end

outfile = joinpath(OUT, "irf_chiib_comparison_oil_shock.png")
savefig(plt, outfile)
@printf "\nSaved: %s\n" outfile

# Console diagnostic: tail decay ratio (dominant eigenvalue proxy) per run
@printf "\n%-12s %10s %10s %12s\n" "phi_b" "Q(h40)" "TB(h40)" "Q tail-ratio"
for (cb, _, df) in loaded
    q  = sort(df[df.variable .== "Q",  :], :period).value
    tb = sort(df[df.variable .== "TB", :], :period).value
    H  = length(q)
    ratios = [q[h+1]/q[h] for h in H-6:H-1 if abs(q[h]) > 1e-9]
    lam = isempty(ratios) ? NaN : sum(ratios)/length(ratios)
    @printf "%-12g %10.4f %10.4f %12.5f\n" cb q[H] tb[H] lam
end
