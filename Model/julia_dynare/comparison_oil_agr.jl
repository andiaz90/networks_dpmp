"""
comparison_oil_agr.jl
=====================
Standalone script that compares the oil price shock and agriculture TFP shock
(El Niño) IRFs side by side.  Reads pre-computed IRF CSVs produced by:
  - oil_shock_analysis.jl  →  tables/oil_shock_irfs.csv
  - agr_shock_analysis.jl  →  tables/agr_shock_irfs.csv

Generates comparison figures saved to:
  figures/agr_shock/                     (local)
  ~/…/Overleaf/Network DPMP-DME/Figures/agr_shock/  (Overleaf sync)

Usage:
  julia --project=. comparison_oil_agr.jl
"""

using CSV, DataFrames, Printf

SCRIPT_DIR = @__DIR__
TABLES_DIR = joinpath(SCRIPT_DIR, "tables")
FIGURES_DIR = joinpath(SCRIPT_DIR, "figures", "agr_shock")
mkpath(FIGURES_DIR)

# Overleaf Dropbox sync folder
OVERLEAF_DIR = get(ENV, "OVERLEAF_AGR",
    abspath(joinpath(homedir(), "Library", "CloudStorage",
        "Dropbox", "Apps", "Overleaf", "Network DPMP-DME", "Figures", "agr_shock")))
overleaf_ok = try mkpath(OVERLEAF_DIR); true catch; false end

@printf "  Local figures  → %s\n" FIGURES_DIR
overleaf_ok && @printf "  Overleaf sync  → %s\n" OVERLEAF_DIR

# =========================================================================== #
#  LOAD PLOTS                                                                  #
# =========================================================================== #

try
    @eval Main using Plots
    @eval Main gr(dpi=200)
    @eval Main Plots.default(
        titlefontsize  = 14,
        guidefontsize  = 13,
        tickfontsize   = 12,
        legendfontsize = 12,
    )
catch e
    error("Plots.jl required. Install with: using Pkg; Pkg.add(\"Plots\")\n$e")
end

using Logging

function save_fig(p, fname)
    local_path = joinpath(FIGURES_DIR, fname)
    try
        redirect_stderr(devnull) do
            with_logger(NullLogger()) do
                Plots.savefig(p, local_path)
            end
        end
    catch
        Plots.savefig(p, local_path)
    end
    @printf "  Saved: %s\n" fname
    if overleaf_ok
        ol_path = joinpath(OVERLEAF_DIR, fname)
        cp(local_path, ol_path; force=true)
        @printf "  → Overleaf: %s\n" ol_path
    end
end

# =========================================================================== #
#  LOAD IRF DATA                                                               #
# =========================================================================== #

oil_path = joinpath(TABLES_DIR, "oil_shock_irfs.csv")
agr_path = joinpath(TABLES_DIR, "agr_shock_irfs.csv")

isfile(oil_path) || error("Oil IRF data not found: $oil_path\n  Run oil_shock_analysis.jl first.")
isfile(agr_path) || error("Agriculture IRF data not found: $agr_path\n  Run agr_shock_analysis.jl first.")

df_oil = CSV.read(oil_path, DataFrame)
df_agr = CSV.read(agr_path, DataFrame)

@printf "\n  Oil IRFs: %d rows  |  Agr IRFs: %d rows\n\n" nrow(df_oil) nrow(df_agr)

# Determine IRF horizon from the data
n_irf = maximum(df_oil.period)
periods = 1:n_irf

function get_irf(df::DataFrame, varname::String)
    sub = filter(r -> String(r.variable) == varname, df)
    isempty(sub) && return zeros(n_irf)
    s = sort(sub, :period)
    out = zeros(n_irf)
    n = min(n_irf, nrow(s))
    out[1:n] .= s.value[1:n]
    return out
end

# Sector names (short, for subplot titles)
short_names = [
    "Agric. & Fishing", "Mining", "Manufacturing", "Utilities",
    "Construction", "Trade & Hotels", "Transport & ICT", "Finance",
    "Real Estate", "Business Serv.", "Personal Serv.", "Public Admin.",
]

# =========================================================================== #
#  COMPUTE SECTORAL INFLATION FROM PRICE LEVELS                                #
#                                                                              #
#  IRF CSVs store PH_i (price level % dev from SS) and pi (CPI inflation).    #
#  Sectoral home-price inflation ≈ ΔPH_i + π_agg, annualized (×4).           #
# =========================================================================== #

function compute_inflation_irf(df::DataFrame, sector::Int)
    ph = get_irf(df, "PH_$(sector)")
    pi_cpi = get_irf(df, "pi")
    dph = similar(ph)
    dph[1] = ph[1] + pi_cpi[1]
    for t in 2:n_irf
        dph[t] = ph[t] - ph[t-1] + pi_cpi[t]
    end
    return dph .* 4   # annualize
end


# =========================================================================== #
#  FIGURE 1: Aggregate Comparison (2×2)                                        #
#  GDP, CPI Inflation, Real Exchange Rate, Policy Rate                        #
# =========================================================================== #

@printf "--- Generating comparison figures ---\n"

p_comp = Plots.plot(layout=(2, 2), size=(1000, 700),
    plot_title="Oil Price Shock vs. Agriculture TFP Shock (El Niño)",
    titlefontsize=10, margin=5Plots.mm)

# GDP
Plots.plot!(p_comp, periods, get_irf(df_oil, "GDP"), subplot=1,
    label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
Plots.plot!(p_comp, periods, get_irf(df_agr, "GDP"), subplot=1,
    label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
    ylabel="% dev.", title="GDP")
Plots.hline!(p_comp, [0.0], subplot=1, color=:black, lw=0.5, ls=:dot, label="")

# CPI Inflation (annualized)
Plots.plot!(p_comp, periods, get_irf(df_oil, "pi") .* 4, subplot=2,
    label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
Plots.plot!(p_comp, periods, get_irf(df_agr, "pi") .* 4, subplot=2,
    label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
    title="CPI Inflation (ann. pp)")
Plots.hline!(p_comp, [0.0], subplot=2, color=:black, lw=0.5, ls=:dot, label="")

# Real Exchange Rate
Plots.plot!(p_comp, periods, get_irf(df_oil, "Q"), subplot=3,
    label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
Plots.plot!(p_comp, periods, get_irf(df_agr, "Q"), subplot=3,
    label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
    xlabel="Quarters", ylabel="% dev.", title="Real Exchange Rate (Q)")
Plots.hline!(p_comp, [0.0], subplot=3, color=:black, lw=0.5, ls=:dot, label="")

# Policy Rate (annualized)
Plots.plot!(p_comp, periods, get_irf(df_oil, "r") .* 4, subplot=4,
    label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
Plots.plot!(p_comp, periods, get_irf(df_agr, "r") .* 4, subplot=4,
    label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
    xlabel="Quarters", title="Policy Rate (ann. pp)")
Plots.hline!(p_comp, [0.0], subplot=4, color=:black, lw=0.5, ls=:dot, label="")

save_fig(p_comp, "irf_comparison_oil_agr.png")


# =========================================================================== #
#  FIGURE 2: Sectoral Inflation Comparison (4×3)                               #
# =========================================================================== #

p_comp_sec = Plots.plot(layout=(4, 3), size=(1200, 900),
    plot_title="Sectoral Inflation: Oil vs. Agriculture TFP Shock (El Niño)",
    titlefontsize=8)

for i in 1:12
    pi_oil_i = compute_inflation_irf(df_oil, i)
    pi_agr_i = compute_inflation_irf(df_agr, i)
    Plots.plot!(p_comp_sec, periods, pi_oil_i, subplot=i,
        label=(i == 1 ? "Oil (10%)" : ""), color=:darkorange, lw=1.8, ls=:solid)
    Plots.plot!(p_comp_sec, periods, pi_agr_i, subplot=i,
        label=(i == 1 ? "Agr TFP (−10%)" : ""), color=:steelblue, lw=1.8, ls=:dash,
        title=short_names[i], titlefontsize=7,
        ylabel=(i % 3 == 1 ? "ann. pp" : ""))
    Plots.hline!(p_comp_sec, [0.0], subplot=i, color=:black, lw=0.5, ls=:dot, label="")
end

save_fig(p_comp_sec, "irf_comparison_sectoral_inflation.png")


# =========================================================================== #
#  FIGURE 3: GDP Gap Comparison (single panel)                                 #
# =========================================================================== #

# Check if GDPgap is available in the IRF data
_gdpgap_oil = get_irf(df_oil, "GDPgap")
_gdpgap_agr = get_irf(df_agr, "GDPgap")

if any(abs.(_gdpgap_oil) .> 1e-12) || any(abs.(_gdpgap_agr) .> 1e-12)
    p_gap = Plots.plot(size=(700, 400),
        title="GDP Gap: Oil vs. Agriculture TFP Shock (El Niño)",
        titlefontsize=11, margin=5Plots.mm,
        xlabel="Quarters", ylabel="% dev. from SS",
        legend=:topright)
    Plots.plot!(p_gap, periods, _gdpgap_oil,
        label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
    Plots.plot!(p_gap, periods, _gdpgap_agr,
        label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash)
    Plots.hline!(p_gap, [0.0], color=:black, lw=0.6, ls=:dash, label="")
    save_fig(p_gap, "irf_comparison_gdpgap.png")
else
    @printf "  [GDPgap not found in IRF data — skipping GDP gap comparison]\n"
end


# =========================================================================== #
#  FIGURE 4: Sectoral Output Comparison (4×3)                                  #
# =========================================================================== #

p_comp_y = Plots.plot(layout=(4, 3), size=(1200, 900),
    plot_title="Sectoral Output: Oil vs. Agriculture TFP Shock (El Niño)",
    titlefontsize=8)

for i in 1:12
    y_oil_i = get_irf(df_oil, "Y_$(i)")
    y_agr_i = get_irf(df_agr, "Y_$(i)")
    Plots.plot!(p_comp_y, periods, y_oil_i, subplot=i,
        label=(i == 1 ? "Oil (10%)" : ""), color=:darkorange, lw=1.8, ls=:solid)
    Plots.plot!(p_comp_y, periods, y_agr_i, subplot=i,
        label=(i == 1 ? "Agr TFP (−10%)" : ""), color=:steelblue, lw=1.8, ls=:dash,
        title=short_names[i], titlefontsize=7,
        ylabel=(i % 3 == 1 ? "% dev." : ""))
    Plots.hline!(p_comp_y, [0.0], subplot=i, color=:black, lw=0.5, ls=:dot, label="")
end

save_fig(p_comp_y, "irf_comparison_sectoral_output.png")


# =========================================================================== #
#  DONE                                                                        #
# =========================================================================== #

@printf "\n%s\n" repeat("=", 60)
@printf "  Comparison figures complete.\n"
@printf "  Local   → %s\n" FIGURES_DIR
overleaf_ok && @printf "  Overleaf → %s\n" OVERLEAF_DIR
@printf "%s\n\n" repeat("=", 60)
