"""
oil_shock_analysis.jl
=====================
NK-IOSOE 12-sector model for Chile — Oil Price Shock Analysis

Computes the response of the Chilean economy to a 10% world oil price shock,
decomposes the transmission into direct cost and network amplification channels,
and generates publication-quality figures (PDF) and LaTeX tables.

Workflow:
  1.  Load sector-level data (same as main_SOE_gap.jl)
  2.  Set structural parameters (from SMM estimates or defaults)
  3.  Configure Exercise 4: oil price shock only (eps_postar)
  4.  Solve steady state and write params_jl.mod
  5.  Run Dynare (subprocess) to get decision rules
  6.  Compute IRFs to a 10% oil price shock
  7.  Decompose: direct oil cost channel vs. IO network amplification
  8.  Generate PDF figures and LaTeX tables
  9.  Run "no-network" counterfactual (diagonal IO matrix)

Usage:
  julia --project=. oil_shock_analysis.jl

REQUIREMENTS:
  Same as main_SOE_gap.jl (Julia >= 1.9, Dynare.jl, CSV, DataFrames, etc.)

OUTPUT FILES (saved to figures/oil_shock/ and tables/):
  Figures:
    irf_aggregate_oil_shock.png     — GDP, inflation, RER, trade balance
    irf_sectoral_Y_oil_shock.png    — 12-sector output IRFs
    irf_sectoral_PH_oil_shock.png   — 12-sector price IRFs
    irf_sectoral_MC_oil_shock.png   — 12-sector marginal cost IRFs
    decomposition_output_oil.png        — sectoral output and MC at impact
    decomposition_mc_inflation_oil.png     — MC decomp (direct+network, stacked) + impact inflation diamonds
    decomposition_mc_inflation_6m_oil.png  — same bars + 6-month cumulative inflation diamonds
    decomposition_mc_inflation_12m_oil.png — same bars + 12-month cumulative inflation diamonds
    oil_intensity_exposure.png          — sector oil exposure map
  Tables:
    oil_shock_decomposition.tex     — LaTeX table of direct/network/total effects
    oil_shock_aggregates.tex        — Aggregate responses summary
"""

# =========================================================================== #
#  PACKAGES                                                                    #
# =========================================================================== #

using LinearAlgebra
using Statistics
using Printf
using Logging
using NLsolve
using CSV
using DataFrames
using XLSX
using StatsBase

# =========================================================================== #
#  INCLUDE HELPERS                                                              #
# =========================================================================== #

SCRIPT_DIR = @__DIR__

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "shock_plots_common.jl"))   # shared figure/table/decomposition routines

# Plotting (graceful failure if not installed)
const _HAS_PLOTS = Ref(false)
try
    @eval Main using Plots
    @eval Main using StatsPlots
    @eval Main gr(dpi=200)
    # Global font-size defaults (applied to every plot unless overridden)
    @eval Main Plots.default(
        titlefontsize  = 14,   # subplot / figure titles
        guidefontsize  = 13,   # axis labels (ylabel, xlabel)
        tickfontsize   = 12,   # tick labels
        legendfontsize = 12,   # legend text
        annotationfontsize = 12,
    )
    _HAS_PLOTS[] = true
catch
    @warn "Plots.jl or StatsPlots.jl not available — figures will be skipped. Install with: ] add Plots StatsPlots"
end

function _main()

@printf "\n%s\n" repeat("=", 70)
@printf "  NK-SOE 12-sector model — Oil Price Shock Analysis (10%%)\n"
@printf "%s\n\n" repeat("=", 70)


# =========================================================================== #
#  PATHS                                                                       #
# =========================================================================== #

MOD_DIR    = joinpath(SCRIPT_DIR, "mod")
REPO_ROOT  = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR   = joinpath(REPO_ROOT, "Data")

FIGURES_DIR = joinpath(SCRIPT_DIR, "figures", "oil_shock")
TABLES_DIR  = joinpath(SCRIPT_DIR, "tables")
mkpath(FIGURES_DIR)
mkpath(TABLES_DIR)

# Overleaf Dropbox sync — mirrors the mfg_shock_analysis.jl pattern.
# Figures go to Figures/oil_shock/ and tables to tables/ inside the Overleaf project.
OVERLEAF_ROOT    = get(ENV, "OVERLEAF_ROOT",
    abspath(joinpath(homedir(), "Library", "CloudStorage",
        "Dropbox", "Apps", "Overleaf", "Network DPMP-DME")))
OVERLEAF_FIG_DIR = joinpath(OVERLEAF_ROOT, "Figures", "oil_shock")
OVERLEAF_TAB_DIR = joinpath(OVERLEAF_ROOT, "tables")
overleaf_ok = try
    mkpath(OVERLEAF_FIG_DIR)
    mkpath(OVERLEAF_TAB_DIR)
    true
catch
    false
end
if overleaf_ok
    @printf "  Overleaf figures : %s\n" OVERLEAF_FIG_DIR
    @printf "  Overleaf tables  : %s\n\n" OVERLEAF_TAB_DIR
else
    @printf "  [Overleaf folder not found — saving locally only]\n\n"
end


# =========================================================================== #
#  THE MODEL IS READ, NOT REBUILT  (refactored 2026-08-20)                     #
# =========================================================================== #
#
# This script used to rebuild the whole calibration and steady state, overwrite
# mod/params_jl.mod, and run Dynare a second time — one of five near-identical
# copies of the same 400 lines, each drifting independently from the model
# main_SOE_gap.jl actually solves. See the header of read_model_params() in
# shock_plots_common.jl for what that drift had reached.
#
# It now READS mod/params_jl.mod and the decision rules main just wrote. IRFs
# are linear and shock-specific, so the relevant ghu column of main's decision
# rule IS this shock's impulse response: no steady state is solved here and
# Dynare is not re-run. A shock script may SELECT and SCALE a shock; it may not
# change a parameter. Anything that changes a parameter is a different model and
# belongs in main_SOE_gap.jl.
#
# Run main_SOE_gap.jl first.
# =========================================================================== #

MP = read_model_params(MOD_DIR, DATA_DIR)
report_model_source(MP)

(; nsec, names_vec, epsilon, beta_val, gamma, psi, subsMC_val, mc_over_ph,
   etastar_val, kappaw_val, phi_b_val, ombar_val,
   modalpha, modalphaV, modalphaK, modepsY, modepsM, modkappa, modalphaOil,
   modcl, goods, services, modbeta,
   isigma_tfp_val, rho_tfp1_val, sigma_om_vec,
   rho_postar_val, sigma_postar_val, epsilonV_oil_val,
   pH_ss, MCi_ss, Yi_ss, L_ss, M_ss, Vi_ss, PMi_ss, P_ss,
   C_gi_ss, C_si_ss, PIV_ss_vec, PL_ss,
   Q_ss, GDP_ss, TB_ss, C_ss, w_ss, N_ss, PV_ss, epsY_baseline) = MP


# Oil-specific: the ghu column for eps_postar is the response to a ONE-S.D.
# innovation. If the shock was off in the run that produced these rules that
# column is identically zero and every figure would be a flat line at zero.
MP.shock_eps_postar_val > 0 || error("""
    shock_eps_postar = 0 in mod/params_jl.mod: the oil shock was OFF in the run
    that produced these decision rules, so its IRF is identically zero. Turn it
    on in main_SOE_gap.jl and re-run before plotting an oil shock.
    """)
sigma_postar_val > 0 || error("sigma_postar = 0 in mod/params_jl.mod — cannot scale the oil IRF to a 10% shock.")

# =========================================================================== #
#  READ DECISION RULES                                                         #
# =========================================================================== #

df_names   = CSV.read(joinpath(MOD_DIR, "dynare_endo_names.csv"), DataFrame)
df_ss_dyn  = CSV.read(joinpath(MOD_DIR, "dynare_ss.csv"), DataFrame)
df_g1_1    = CSV.read(joinpath(MOD_DIR, "dynare_g1_1.csv"), DataFrame)
df_g1_2    = CSV.read(joinpath(MOD_DIR, "dynare_g1_2.csv"), DataFrame)
df_states  = CSV.read(joinpath(MOD_DIR, "dynare_state_rows.csv"), DataFrame)

endo_names = String.(df_names.variable)
ss_vec     = Float64.(df_ss_dyn.ss_value)
ghx        = Matrix{Float64}(df_g1_1)
ghu        = Matrix{Float64}(df_g1_2)
state_rows = Int.(df_states.state_row)

n_endo   = length(endo_names)
n_states = length(state_rows)
n_shocks = size(ghu, 2)

endo_idx = Dict(nm => i for (i, nm) in enumerate(endo_names))

@printf "  Variables: %d  States: %d  Shocks: %d\n" n_endo n_states n_shocks

# Validate decision rules — detect Dynare solver failures (zero matrices)
_nz_ghx = sum(abs.(ghx) .> 1e-12)
_nz_ghu = sum(abs.(ghu) .> 1e-12)
@printf "  Decision rule nonzeros:  ghx=%d  ghu=%d\n" _nz_ghx _nz_ghu
if _nz_ghx == 0 || _nz_ghu == 0
    @error """
    ╔══════════════════════════════════════════════════════════════╗
    ║  DECISION RULES ARE ALL ZEROS — Dynare solver failed!      ║
    ║  This parameterization likely hit the ARM/aarch64 gees bug  ║
    ║  or another Dynare.jl solver failure.                       ║
    ║  IRFs for this εY case will be meaningless (flat at zero).  ║
    ╚══════════════════════════════════════════════════════════════╝
    """
end

# Identify the eps_postar column BY NAME (robust to .mod shock-order changes;
# in the unified model eps_postar is column 4, not 5).
postar_col = exo_col("eps_postar", MOD_DIR)
@printf "  Oil shock (eps_postar) column: %d\n\n" postar_col


# =========================================================================== #
#  COMPUTE IRFs: 10% OIL PRICE SHOCK                                           #
#                                                                              #
#  The shock size is chosen so that POstar jumps by 10% on impact.             #
#  In the .mod file: log(POstar/POstar_ss) = rho*log(POstar(-1)/POstar_ss)    #
#                                            + sigma_postar * eps_postar       #
#  For a 10% shock: sigma_postar * eps_postar = log(1.10)                     #
#  Since sigma_postar = 0.02, we need eps_postar = log(1.10) / 0.02 ≈ 4.76   #
#  The IRFs are linear, so we scale: IRF(10%) = scale_factor * IRF(1 s.d.)    #
# =========================================================================== #

shock_pct = 0.10   # 10% oil price increase
scale_factor = log(1 + shock_pct) / sigma_postar_val   # units of s.d.

@printf "--- Computing IRFs (10%% oil price shock) ---\n"
@printf "  Scale factor: %.4f std devs (to get %.0f%% POstar increase)\n\n" scale_factor 100*shock_pct

n_irf = 80   # quarters (extended from 40 to assess convergence of the slow NFA mode)

# nT (the plotting horizon) MUST equal n_irf. It used to be set independently
# in the reader block and silently disagreed for the four TFP scripts
# (nT=80 vs n_irf=40), which is a BoundsError inside the figure code.
nT = n_irf

# State-space matrices
A = ghx[state_rows, :]  # n_states × n_states
B = ghu[state_rows, :]  # n_states × n_shocks

# Compute scaled IRFs
irf_mat = zeros(n_endo, n_irf)   # variable × period
x = zeros(n_states)

for h in 1:n_irf
    if h == 1
        y_h = ghu[:, postar_col] .* scale_factor
        x   = B[:, postar_col]   .* scale_factor
    else
        y_h = ghx * x
        x   = A * x
    end
    irf_mat[:, h] .= y_h
end

# Convert to percentage deviations from SS
# For variables with SS ≈ 0 (gaps, log-deviations), the linearized IRF is already
# in deviation units — multiply by 100 to get percentage points.
irf_pct = zeros(n_endo, n_irf)
for i in 1:n_endo
    ssv = abs(ss_vec[i])
    if ssv > 1e-10
        irf_pct[i, :] .= 100.0 .* irf_mat[i, :] ./ ssv
    else
        # Gap variables (Ygap, GDPgap, Ngap, etc.) have SS = 0;
        # linearized IRF is in log-deviation units → ×100 for pp.
        irf_pct[i, :] .= 100.0 .* irf_mat[i, :]
    end
end

# Helper to extract IRF for a variable
function get_irf_var(varname::String)
    idx = get(endo_idx, varname, 0)
    idx == 0 && return zeros(n_irf)
    return irf_pct[idx, :]
end

# =========================================================================== #
#  EXTRACT IRFs, DECOMPOSE, AND GENERATE OUTPUTS (shared module)               #
# =========================================================================== #

# Aggregate IRFs
gdp_irf = get_irf_var("GDP");  pi_irf = get_irf_var("pi")
q_irf   = get_irf_var("Q");    tb_irf = get_irf_var("TB")
c_irf   = get_irf_var("C");    cg_irf = get_irf_var("C_g"); cs_irf = get_irf_var("C_s")
n_irf_v = get_irf_var("N");    r_irf  = get_irf_var("r");   w_irf  = get_irf_var("w")
po_irf  = get_irf_var("PO");   pv_irf = get_irf_var("PV")
gdpgap_irf = get_irf_var("GDPgap")

# Sectoral IRF matrices (nsec × n_irf)
y_irf_mat    = vcat([permutedims(get_irf_var("Y_$(i)"))    for i in 1:nsec]...)
ph_irf_mat   = vcat([permutedims(get_irf_var("PH_$(i)"))   for i in 1:nsec]...)
mc_irf_mat   = vcat([permutedims(get_irf_var("MC_$(i)"))   for i in 1:nsec]...)
l_irf_mat    = vcat([permutedims(get_irf_var("L_$(i)"))    for i in 1:nsec]...)
ygap_irf_mat = vcat([permutedims(get_irf_var("Ygap_$(i)")) for i in 1:nsec]...)

y_irf_impact  = y_irf_mat[:, 1]
mc_irf_impact = mc_irf_mat[:, 1]
ph_irf_impact = ph_irf_mat[:, 1]

# Real sectoral VALUE ADDED (GDP-consistent). Reported as "Y_i" in the
# decomposition figure/table: gross output double-counts intermediates and can
# rise under a negative TFP shock even as value added (and GDP) falls.
va_irf_mat    = sectoral_va_irf(get_irf_var, nsec, n_irf, Yi_ss, M_ss, Vi_ss, pH_ss, PMi_ss, PV_ss)
va_irf_impact = va_irf_mat[:, 1]

# Sectoral & group home-price inflation (shared helper)
infl = compute_inflation_aggregates(ph_irf_mat, pi_irf, C_gi_ss, C_si_ss, nsec, n_irf)

# Direct vs. network MC decomposition (Leontief): per-sector oil cost push
d_log_PO  = log(1 + shock_pct) * 100
direct_mc = (modalphaV .* modalphaOil) .* d_log_PO
dec = leontief_decomp(direct_mc, modalpha, modbeta, Yi_ss)

@printf "\n  %-30s %8s %8s %8s\n" "Sector" "Direct" "Network" "Total"
for i in 1:nsec
    @printf "  %-30s %7.3f%% %7.3f%% %7.3f%%\n" names_vec[i] dec.direct_mc[i] dec.network_mc[i] dec.total_mc[i]
end
@printf "  %-30s %7.3f%% %7.3f%% %7.3f%%\n\n" "Aggregate (output-wtd)" dec.agg_direct dec.agg_network dec.agg_total

# 5-way GE marginal-cost decomposition at h = 1, 2, 4
alpha_L_vec = 1.0 .- modalpha .- modalphaV .- modalphaK
# Rental IRFs feed the capital bar of the marginal-cost decomposition.
rk_irf_mat = reduce(vcat, [get_irf_var("RK_$(i)")' for i in 1:nsec])
ge_kwargs = (po_irf=po_irf, pv_irf=pv_irf, modalphaOil=modalphaOil)
ge_h1 = ge_mc_components(:price, 1, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_h2 = ge_mc_components(:price, 2, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_h4 = ge_mc_components(:price, 4, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_fun = h -> ge_mc_components(:price, h, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_colors, ge_labels = ge_component_style(:price)

# "Affected" = sectors with above-average DIRECT oil cost share (αV·αOil).
# Used for the affected-vs-rest inflation figure. The oil shock hits every
# sector, so the inflation decompositions use a SINGLE axis (no pi_decomp_axis2).
oil_affected = findall((modalphaV .* modalphaOil) .> mean(modalphaV .* modalphaOil))

@printf "--- Aggregate IRFs (10%% oil shock, impact) ---\n"
@printf "  GDP %+.3f%%   pi %+.3f ann.pp   Q %+.3f%%   r %+.3f ann.pp\n\n" gdp_irf[1] (pi_irf[1]*4) q_irf[1] (r_irf[1]*4)

# Assemble context and generate the full figure + table set via the shared module
ctx = (
    # Tag carries the phi_b override so sensitivity runs don't overwrite the
    # baseline figures/tables (e.g. tag = "oil_chiib0p01").
    HAS_PLOTS = _HAS_PLOTS[],
    tag = "oil" *
        (haskey(ENV, "PHIB_OVERRIDE")  ? "_chiib" * replace(ENV["PHIB_OVERRIDE"],  "." => "p") : "") *
        (haskey(ENV, "KAPPAW_OVERRIDE") ? "_kw"    * replace(ENV["KAPPAW_OVERRIDE"], "." => "p") : ""),
    title_long = "Shock al Precio del Petróleo (+10%)", title_short = "Shock Petróleo",
    nsec = nsec, nT = nT, names_vec = names_vec, goods = goods,
    TB_ss = TB_ss, GDP_ss = GDP_ss,
    FIGURES_DIR = FIGURES_DIR, TABLES_DIR = TABLES_DIR,
    OVERLEAF_FIG_DIR = OVERLEAF_FIG_DIR, OVERLEAF_TAB_DIR = OVERLEAF_TAB_DIR,
    overleaf_ok = overleaf_ok, fignames = oil_fignames(),
    epsY_baseline = epsY_baseline,
    baseline_label = "Base (εY=$(round(epsY_baseline,digits=2)))",
    baseline_color = IPOM_NAVY, baseline_lw = 2.5,   # IPoM palette (shock_plots_common.jl)
    gdp_irf=gdp_irf, pi_irf=pi_irf, q_irf=q_irf, tb_irf=tb_irf, c_irf=c_irf,
    cg_irf=cg_irf, cs_irf=cs_irf, n_irf_v=n_irf_v, r_irf=r_irf, w_irf=w_irf,
    gdpgap_irf=gdpgap_irf,
    y_irf_mat=y_irf_mat, ph_irf_mat=ph_irf_mat, mc_irf_mat=mc_irf_mat,
    l_irf_mat=l_irf_mat, ygap_irf_mat=ygap_irf_mat,
    y_irf_impact=y_irf_impact, mc_irf_impact=mc_irf_impact, ph_irf_impact=ph_irf_impact,
    va_irf_impact=va_irf_impact, va_irf_mat=va_irf_mat,
    pi_sec_mat=infl.pi_sec_mat, pi_agg_irf=infl.pi_agg_irf,
    pi_goods_irf=infl.pi_goods_irf, pi_serv_irf=infl.pi_serv_irf,
    cons_ss_all=infl.cons_ss_all,
    affected_sectors=oil_affected,
    affected_label="Sectores intensivos en petróleo",
    infl_irf_impact=infl.infl_irf_impact, infl_6m=infl.infl_6m, infl_12m=infl.infl_12m,
    direct_mc=dec.direct_mc, network_mc=dec.network_mc, total_mc=dec.total_mc,
    amp_ratio=dec.amp_ratio, agg_direct=dec.agg_direct, agg_network=dec.agg_network,
    agg_total=dec.agg_total, agg_amp=dec.agg_amp,
    ge_h1=ge_h1, ge_h2=ge_h2, ge_h4=ge_h4, ge_fun=ge_fun, ge_colors=ge_colors, ge_labels=ge_labels,
    exposure_vec = (modalphaV .* modalphaOil) .* 100,
    exposure_label = "Participación del petróleo en insumos totales (%)",
    exposure_ylabel = "Participación de costo de petróleo (%)",
    exposure_title = "Participación del Petróleo en Insumos por Sector",
    share_vec = modalphaOil, modalphaV = modalphaV, Yi_ss = Yi_ss,
    decomp_share_head = "Petróleo", decomp_share_sub = "Part. (\\%)",
    decomp_caption = "Descomposición Sectorial de un Shock de 10\\% al Precio del Petróleo",
    decomp_notes = "La participación del petróleo es \$\\alpha^{\\text{Oil}}_i\$, la fracción de las importaciones intermedias del sector \$i\$ que corresponde a petróleo/combustibles (matriz IP de Chile 2021). El CM directo es \$\\alpha_{Vi}\\,\\alpha^{\\text{Oil}}_i\\,\\Delta\\log P^O\$. La red es el costo adicional vía encadenamientos IP bajo la aproximación de Leontief de equilibrio parcial. Total = Directo + Red. Razón de amplificación = Total/Directo. La IRF de producto es la respuesta sectorial en el período de impacto. B\\,=\\,Bienes, S\\,=\\,Servicios.",
    agg_caption = "Respuestas Agregadas a un Shock de 10\\% al Precio del Petróleo",
    agg_notes = "PIB, consumo, empleo y TCR en \\% de desviación del estado estacionario; inflación y tasa de política en pp anualizados; balanza comercial en pp del PIB. El shock es un aumento único de 10\\% en el precio mundial del petróleo \$P^{O*}_t\$ con persistencia \$\\rho = $(rho_postar_val)\$. Impacto = respuesta en el trimestre 1. Mín/Máx = respuesta extrema en 40 trimestres.",
    get_irf = get_irf_var,
)

generate_shock_outputs(ctx)



# =========================================================================== #
#  DONE                                                                        #
# =========================================================================== #

@printf "\n%s\n" repeat("=", 70)
@printf "  Oil shock analysis complete.\n"
@printf "  Figures (local)   → %s\n" FIGURES_DIR
@printf "  Tables  (local)   → %s\n" TABLES_DIR
if overleaf_ok
    @printf "  Figures (Overleaf) → %s\n" OVERLEAF_FIG_DIR
    @printf "  Tables  (Overleaf) → %s\n" OVERLEAF_TAB_DIR
else
    @printf "  Overleaf sync     → not available (folder not found)\n"
end
@printf "%s\n\n" repeat("=", 70)

end  # function _main()

Base.invokelatest(_main)
