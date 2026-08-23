"""
agrmin_shock_analysis.jl
========================
NK-IOSOE 12-sector model for Chile — Joint Agriculture + Mining TFP Shock Analysis

Applies a SIMULTANEOUS -1% total-factor-productivity shock to both Agriculture
(sector 1) and Mining (sector 2) and traces the joint transmission through the
input-output network, mirroring the single-sector TFP shock scripts.

Computes the response of the Chilean economy to a SIMULTANEOUS -1% TFP shock
to Agriculture (sector 1) and Mining (sector 2), decomposes the transmission
into direct cost and network amplification channels, and generates
publication-quality figures (PDF) and LaTeX tables.

Workflow:
  1.  Load sector-level data (same as main_SOE_gap.jl)
  2.  Set structural parameters (from SMM estimates or defaults)
  3.  Turn on epsA_1 and epsA_2 (Agriculture + Mining TFP) simultaneously
  4.  Solve steady state and write params_jl.mod
  5.  Run Dynare (subprocess) to get decision rules
  6.  Compute IRFs to the joint -1% TFP shock (superpose the two scaled impulses)
  7.  Decompose: direct own-cost push vs. IO network amplification (Leontief)
  8.  Generate PDF figures and LaTeX tables (shared shock_plots_common routines)

Usage:
  julia --project=. agrmin_shock_analysis.jl

REQUIREMENTS:
  Same as main_SOE_gap.jl (Julia >= 1.9, Dynare.jl, CSV, DataFrames, etc.)

OUTPUT FILES (saved to figures/agrmin_shock/ and tables/):
  Figures (tag = agrmin):
    irf_aggregate_agrmin_shock.png        — GDP, inflation, RER, trade balance
    irf_sectoral_{Y,PH,MC,L}_agrmin_shock.png — 12-sector IRFs
    decomposition_mc_inflation_{,6m_,12m_}agrmin_baseline.png  — 5-way GE MC decomp
    decomposition_mc3_inflation_{,6m_,12m_}agrmin_baseline.png — 3-way direct/indirect/others
    decomposition_output_agrmin_baseline.png — sectoral value added + inflation at impact
    agrmin_intensity_exposure.png         — combined Agriculture+Mining input exposure
  Tables:
    agrmin_shock_decomposition.tex        — direct/network/total effects
    agrmin_shock_aggregates.tex           — aggregate responses summary
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
@printf "  NK-SOE 12-sector model — Joint Agriculture + Mining TFP Shock Analysis (-1%%)\n"
@printf "%s\n\n" repeat("=", 70)


# =========================================================================== #
#  PATHS                                                                       #
# =========================================================================== #

MOD_DIR    = joinpath(SCRIPT_DIR, "mod")
REPO_ROOT  = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR   = joinpath(REPO_ROOT, "Data")

FIGURES_DIR = joinpath(SCRIPT_DIR, "figures", "agrmin_shock")
TABLES_DIR  = joinpath(SCRIPT_DIR, "tables")
mkpath(FIGURES_DIR)
mkpath(TABLES_DIR)

# Overleaf Dropbox sync — mirrors the mfg_shock_analysis.jl pattern.
# Figures go to Figures/oil_shock/ and tables to tables/ inside the Overleaf project.
OVERLEAF_ROOT    = get(ENV, "OVERLEAF_ROOT",
    abspath(joinpath(homedir(), "Library", "CloudStorage",
        "Dropbox", "Apps", "Overleaf", "Network DPMP-DME")))
OVERLEAF_FIG_DIR = joinpath(OVERLEAF_ROOT, "Figures", "agrmin_shock")
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


# --- shock selection (the ONLY thing this script chooses) ---
shock_sectors = [1, 2]          # both shocked simultaneously
shock_sector  = shock_sectors[1]  # kept for legacy single-sector references

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

# Identify the epsA columns for the shocked sectors BY NAME (robust to .mod
# shock-order changes).
shock_cols = [exo_col("epsA_$(s)", MOD_DIR) for s in shock_sectors]
@printf "  TFP shock columns: %s (sectors %s)\n\n" string(shock_cols) string(shock_sectors)


# =========================================================================== #
#  COMPUTE IRFs: JOINT -1% TFP SHOCK (Agriculture + Mining)                   #
#                                                                              #
#  Each shocked sector's own TFP is driven down by exactly 1% on impact.      #
#  The model is linear, so the joint response is the SUM of the two scaled     #
#  single-sector impulses applied in the SAME period.  For each sector k we    #
#  pick the scale so that A_k drops 1% on impact, then superpose.             #
# =========================================================================== #

shock_pct = -0.01   # 1% negative productivity (TFP) shock per shocked sector
# Per-sector scale: drive each shocked sector's own A down exactly 1% on impact.
scale_factors = Vector{Float64}(undef, length(shock_sectors))
for (idx, s) in enumerate(shock_sectors)
    _a_unit = ghu[get(endo_idx, "A_$(s)", 0), shock_cols[idx]]
    scale_factors[idx] = abs(_a_unit) > 1e-12 ? shock_pct / _a_unit : shock_pct / isigma_tfp_val[s]
    @printf "  A_%d unit response = %.5f -> scale factor = %.4f\n" s _a_unit scale_factors[idx]
end

@printf "\n--- Computing IRFs (joint -1%% TFP shock, sectors %s) ---\n\n" string(shock_sectors)

n_irf = 40   # quarters

# nT (the plotting horizon) MUST equal n_irf. It used to be set independently
# in the reader block and silently disagreed for the four TFP scripts
# (nT=80 vs n_irf=40), which is a BoundsError inside the figure code.
nT = n_irf

# State-space matrices
A = ghx[state_rows, :]  # n_states × n_states
B = ghu[state_rows, :]  # n_states × n_shocks

# Compute scaled IRFs (superpose the two simultaneous impulses at h = 1)
irf_mat = zeros(n_endo, n_irf)   # variable × period
x = zeros(n_states)

for h in 1:n_irf
    if h == 1
        y_h = zeros(n_endo)
        x   = zeros(n_states)
        for (idx, col) in enumerate(shock_cols)
            y_h .+= ghu[:, col] .* scale_factors[idx]
            x   .+= B[:, col]   .* scale_factors[idx]
        end
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

k = shock_sector

# Aggregate IRFs
gdp_irf = get_irf_var("GDP");  pi_irf = get_irf_var("pi")
q_irf   = get_irf_var("Q");    tb_irf = get_irf_var("TB")
c_irf   = get_irf_var("C");    cg_irf = get_irf_var("C_g"); cs_irf = get_irf_var("C_s")
n_irf_v = get_irf_var("N");    r_irf  = get_irf_var("r");   w_irf  = get_irf_var("w")
pv_irf  = get_irf_var("PV")
gdpgap_irf = get_irf_var("GDPgap")

# Sectoral IRF matrices (nsec x n_irf)
y_irf_mat    = vcat([permutedims(get_irf_var("Y_$(i)"))    for i in 1:nsec]...)
ph_irf_mat   = vcat([permutedims(get_irf_var("PH_$(i)"))   for i in 1:nsec]...)
mc_irf_mat   = vcat([permutedims(get_irf_var("MC_$(i)"))   for i in 1:nsec]...)
l_irf_mat    = vcat([permutedims(get_irf_var("L_$(i)"))    for i in 1:nsec]...)
ygap_irf_mat = vcat([permutedims(get_irf_var("Ygap_$(i)")) for i in 1:nsec]...)
a_irf_mat    = vcat([permutedims(get_irf_var("A_$(i)"))    for i in 1:nsec]...)

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

# Direct vs. network MC decomposition (Leontief):
# the joint -1% productivity shock raises BOTH shocked sectors' own marginal
# cost by +1% on impact (= -Ahat_k); zero direct push elsewhere.  Network
# propagation of those cost increases downstream is the Leontief amplification.
direct_mc = -a_irf_mat[:, 1]                 # own-sector cost push (%, impact)
dec = leontief_decomp(direct_mc, modalpha, modbeta, Yi_ss)

@printf "
  %-30s %8s %8s %8s
" "Sector" "Direct" "Network" "Total"
for i in 1:nsec
    @printf "  %-30s %7.3f%% %7.3f%% %7.3f%%
" names_vec[i] dec.direct_mc[i] dec.network_mc[i] dec.total_mc[i]
end
@printf "
"

# Sanity check: confirm both shocked sectors' TFP and output signs
for s in shock_sectors
    @printf "  Verificación shock: A_%d impacto = %+.2f%% (debe ser ≈ -1%%),  VA_%d impacto = %+.2f%%\n" s a_irf_mat[s,1] s va_irf_impact[s]
end

# 5-way GE marginal-cost decomposition at h = 1, 2, 4
alpha_L_vec = 1.0 .- modalpha .- modalphaV .- modalphaK
# Rental IRFs feed the capital bar of the marginal-cost decomposition.
rk_irf_mat = reduce(vcat, [get_irf_var("RK_$(i)")' for i in 1:nsec])
ge_kwargs = (pv_irf=pv_irf, a_irf_mat=a_irf_mat)
ge_h1 = ge_mc_components(:tfp, 1, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_h2 = ge_mc_components(:tfp, 2, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_h4 = ge_mc_components(:tfp, 4, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_fun = h -> ge_mc_components(:tfp, h, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; modalphaK=modalphaK, rk_irf_mat=rk_irf_mat, ge_kwargs...)
ge_colors, ge_labels = ge_component_style(:tfp)

@printf "--- Aggregate IRFs (joint -1%% Agriculture + Mining TFP shock, impact) ---
"
@printf "  GDP %+.3f%%   pi %+.3f ann.pp   Q %+.3f%%   r %+.3f ann.pp

" gdp_irf[1] (pi_irf[1]*4) q_irf[1] (r_irf[1]*4)

# Exposure of each sector to the shocked sectors' output (Agriculture + Mining):
# direct input share and input-cost share summed across both shocked sectors.
exposure_share = vec(sum(modbeta[:, shock_sectors], dims=2))           # direct input share from shocked sectors
exposure_cost  = (modalpha .* vec(sum(modbeta[:, shock_sectors], dims=2))) .* 100   # share of total input cost

ctx = (
    HAS_PLOTS = _HAS_PLOTS[], tag = "agrmin",
    title_long = "Shock Conjunto de Productividad en Agricultura y Minería (-1%)", title_short = "Shock PTF Agric.+Minería",
    nsec = nsec, nT = nT, names_vec = names_vec, goods = goods,
    TB_ss = TB_ss, GDP_ss = GDP_ss,
    FIGURES_DIR = FIGURES_DIR, TABLES_DIR = TABLES_DIR,
    OVERLEAF_FIG_DIR = OVERLEAF_FIG_DIR, OVERLEAF_TAB_DIR = OVERLEAF_TAB_DIR,
    overleaf_ok = overleaf_ok, fignames = standard_fignames("agrmin"),
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
    cons_ss_all=infl.cons_ss_all, affected_sectors=shock_sectors,
    pi_decomp_axis2=shock_sectors,   # shocked sectors on a right-hand axis in the inflation decomps
    affected_label="Agricultura y Minería",
    infl_irf_impact=infl.infl_irf_impact, infl_6m=infl.infl_6m, infl_12m=infl.infl_12m,
    direct_mc=dec.direct_mc, network_mc=dec.network_mc, total_mc=dec.total_mc,
    amp_ratio=dec.amp_ratio, agg_direct=dec.agg_direct, agg_network=dec.agg_network,
    agg_total=dec.agg_total, agg_amp=dec.agg_amp,
    ge_h1=ge_h1, ge_h2=ge_h2, ge_h4=ge_h4, ge_fun=ge_fun, ge_colors=ge_colors, ge_labels=ge_labels,
    exposure_vec = exposure_cost,
    exposure_label = "Participación de costo de insumos desde Agricultura + Minería (%)",
    exposure_ylabel = "Participación de costo (%)",
    exposure_title = "Participación de Insumos desde Agricultura y Minería por Sector",
    share_vec = exposure_share, modalphaV = modalphaV, Yi_ss = Yi_ss,
    decomp_share_head = "Agr+Min", decomp_share_sub = "Insumo (\\%)",
    decomp_caption = "Descomposición Sectorial de un Shock Conjunto Negativo de $(round(abs(shock_pct)*100, digits=1))\\% a la PTF en Agricultura y Minería",
    decomp_notes = "La participación de insumos es la fracción de los insumos intermedios del sector \$i\$ provenientes de Agricultura y Minería combinadas (matriz IP de Chile 2021). El CM directo es el empuje de costo por productividad propia \$-\\hat A_i\$, distinto de cero solo para los sectores afectados (Agricultura y Minería). La red es el costo adicional propagado por encadenamientos IP bajo la aproximación de Leontief de equilibrio parcial. Total = Directo + Red. Razón de amplificación = Total/Directo (mostrada solo para los sectores afectados). La IRF de producto es la respuesta sectorial en el período de impacto. B\\,=\\,Bienes, S\\,=\\,Servicios.",
    agg_caption = "Respuestas Agregadas a un Shock Conjunto Negativo de $(round(abs(shock_pct)*100, digits=1))\\% a la PTF en Agricultura y Minería",
    agg_notes = "PIB, consumo, empleo y TCR en \\% de desviación del estado estacionario; inflación y tasa de política en pp anualizados; balanza comercial en pp del PIB. El shock es una reducción única y simultánea de $(round(abs(shock_pct)*100, digits=1))\\% en la productividad total de factores de los sectores Agricultura y Minería, con persistencia \$\\rho_{Agr} = $(round(rho_tfp1_val[1],digits=3))\$ y \$\\rho_{Min} = $(round(rho_tfp1_val[2],digits=3))\$. Impacto = respuesta en el trimestre 1. Mín/Máx = respuesta extrema en 40 trimestres.",
    get_irf = get_irf_var,
)

generate_shock_outputs(ctx)



# =========================================================================== #
#  DONE                                                                        #
# =========================================================================== #

@printf "\n%s\n" repeat("=", 70)
@printf "  Joint Agriculture + Mining TFP shock analysis complete.\n"
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
