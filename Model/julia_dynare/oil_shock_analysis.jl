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
#  READ DATA (identical to main_SOE_gap.jl)                                    #
# =========================================================================== #

nsec = 12

DATA_CANDIDATES = filter!(!isempty, [DATA_DIR, SCRIPT_DIR])
function find_file(candidates, fnames...)
    for d in candidates, fname in fnames
        p = joinpath(d, fname); isfile(p) && return p
    end; return ""
end

path_cal     = find_file(DATA_CANDIDATES, "sector_calibration.csv")
path_io      = find_file(DATA_CANDIDATES, "IO_2021_chile.csv")
path_fpa     = find_file(DATA_CANDIDATES, "fpa_vector_few_industries_chile.csv")
path_sec_mom = find_file(DATA_CANDIDATES, "sectoral_moments.csv")

for p in [path_cal, path_io, path_fpa]
    isempty(p) && error("Required data file not found. Searched: $(join(DATA_CANDIDATES, ", "))")
end

cal_df     = CSV.read(path_cal, DataFrame)
names_vec  = String.(cal_df.name)
alpha      = Float64.(cal_df.alpha)
alpha_V    = Float64.(cal_df.alpha_V)
var_rho    = Float64.(cal_df.var_rho)
spend_good = Float64.(cal_df.spend_good)
spend_serv = Float64.(cal_df.spend_serv)

if !isempty(path_sec_mom)
    sec_mom = CSV.read(path_sec_mom, DataFrame)
    y_d = Float64.(sec_mom.std_Y); p_d = Float64.(sec_mom.std_PH); l_d = Float64.(sec_mom.std_L)
else
    y_d = zeros(nsec); p_d = zeros(nsec); l_d = zeros(nsec)
end

betaio_df = CSV.read(path_io, DataFrame, header=false)
betaio    = Matrix{Float64}(betaio_df[1:nsec, 1:nsec])
col_sums  = sum(betaio, dims=1)
betax     = betaio ./ col_sums
modbeta   = Matrix(betax')   # modbeta[i,j] = share of inputs sector i gets from j

fpa_df    = CSV.read(path_fpa, DataFrame, header=false)
theta_vec = vec(Matrix{Float64}(fpa_df))

modalpha  = alpha
modalphaV = alpha_V

@printf "  Loaded %d sectors.\n\n" nsec


# =========================================================================== #
#  STRUCTURAL PARAMETERS                                                       #
# =========================================================================== #

beta_val = 0.986;  epsilon = 10.0;  gamma = 2.0;  chi = 1.0
psi = haskey(ENV, "PSI_OVERRIDE") ? parse(Float64, ENV["PSI_OVERRIDE"]) : 0.5   # inverse Frisch; lower = flatter labor supply, wage reacts less
# theta_vec is the FREQUENCY of price adjustment (fraction of firms that reset each
# quarter). Calvo stickiness = probability of NOT adjusting = 1 - theta_vec.
stick    = 1 .- theta_vec
modkappa = stick .* (epsilon - 1) ./ ((1 .- stick) .* (1 .- stick .* beta_val))

goods    = spend_good .> spend_serv
services = spend_serv .> spend_good
modgammag = spend_good ./ sum(spend_good)
modgammas = spend_serv ./ sum(spend_serv)

modcl = zeros(nsec); modclneg = zeros(nsec); modcm = zeros(nsec)
modepsM  = fill(0.1, nsec)
modepsY  = fill(0.8, nsec)
modpsil = fill(-1000.0, nsec); modpsim = fill(-1000.0, nsec)

phi_val = 2.5; rhoi_val = 0.6; rhoirule_val = 0.74
gammaind_val = 0.0; ilabcosts_val = 0.1; ombar_val = 0.57

Pistar_ss = 1.00; Rworld_ss = Pistar_ss / beta_val
kappaV_val = 1e13; epsilonV_val = 1e13; epsilonX_val = 1.0; omegaX_val = 1.0
# Debt-elastic risk premium: chii_b = 0.001 implies NFA eigenvalue ≈ 0.999
# (Q half-life ~600 quarters) — IRFs look non-convergent on a 40q window.
# Override per-run with ENV["CHIIB_OVERRIDE"], e.g. `CHIIB_OVERRIDE=0.01 julia …`.
# chii_b does NOT affect the steady state (premium term is zero at SS).
chii_b_val = haskey(ENV, "CHIIB_OVERRIDE") ? parse(Float64, ENV["CHIIB_OVERRIDE"]) : 0.0024
# Default = XMAS posterior mean (100ψ = 0.24 [0.18,0.30], Bayesian with Chile
# EMBIG observable; NFA/quarterly-GDP units — verified same scaling as here).
etastar_val = 3.5
# Sticky wages (Rotemberg wage PC): kappaw = 0 → flexible wages (pre-existing
# model); default 115 ≈ Calvo wage duration of 4 quarters. Override per-run
# with ENV["KAPPAW_OVERRIDE"], e.g. `KAPPAW_OVERRIDE=0 julia …`.
# Steady state is unchanged by kappaw (SS wage markdown offset by subsidy).
epsw_val   = 10.0
kappaw_val = haskey(ENV, "KAPPAW_OVERRIDE") ? parse(Float64, ENV["KAPPAW_OVERRIDE"]) : 115.0
xi_rstar_val = 0.2; ystar_ss_val = 1.0; PVstar_ss = 1.0; sigmaH_val = 0.999

modchiX   = let   # sectoral export shares chi_i^X from Chilean 2021 supply-use table (Data/computed)
    _f = joinpath(DATA_DIR, "computed", "export_shares_chile.csv")
    isfile(_f) ? (v = Float64.(CSV.read(_f, DataFrame).chi_x); v ./ sum(v)) : fill(1/nsec, nsec)
end
modvarrho = var_rho
modA      = ones(nsec)

# ---- Oil price shock parameters ----
modalphaOil = [0.1635, 0.2164, 0.1890, 0.0871, 0.0417,
               0.0793, 0.3734, 0.0043, 0.0449, 0.0633, 0.0391, 0.0447]
epsilonV_oil_val  = haskey(ENV, "OILSUB_OVERRIDE") ? parse(Float64, ENV["OILSUB_OVERRIDE"]) : 0.5   # oil↔non-oil import substitution (default 0.5: model solves)
rho_postar_val    = 0.9
sigma_postar_val  = 0.02
POstar_ss_val     = 1.0

# ---- Exercise 4: OIL SHOCK ONLY ----
# All other shocks OFF; oil shock ON with 10% innovation
sigma_i_val     = 0.0
rho_om1_val     = 0.1
sigma_om_vec    = zeros(nsec)
rho_tfp1_val    = 0.5; rho_tfp2_val = 0.0
isigma_tfp_val  = zeros(nsec)
rho_val         = 0.1
sigma_L_agg_val = 0.0
rho_pvstar_val  = 0.9; sigma_pvstar_val = 0.0   # import price shock OFF
rho_xi_val      = 0.80; sigma_xi_val = 0.0      # preference shock OFF

# Shock activation flags
shock_eps_om_vec     = zeros(nsec)
shock_eps_i_val      = 0.0
shock_eps_pvstar_val = 0.0
shock_eps_xi_val     = 0.0
shock_eps_postar_val = 1.0   # << OIL SHOCK ON
shock_epsA_val       = ones(nsec)

# Load SMM estimates if available (override structural params)
smm_est_file = joinpath(DATA_DIR, "smm_estimates.csv")
if isfile(smm_est_file)
    est_df = CSV.read(smm_est_file, DataFrame)
    est    = Dict(String(r.param) => Float64(r.value) for r in eachrow(est_df))
    ilabcosts_val    = est["ilabcosts"]
    modepsY          = fill(est["epsY"], nsec)
    modepsM          = fill(est["epsM"], nsec)
    kappaV_val       = exp(est["log_kappaV"])
    rho_om1_val      = est["rho_om"]
    rho_tfp1_val     = est["rho_A"]
    isigma_tfp_val   = [est["isigma_tfp_$(i)"] for i in 1:nsec]
    haskey(est, "etastar") && (etastar_val = est["etastar"])
    @printf "  Loaded SMM estimates from %s\n" smm_est_file
end

# ---- εY (production-input elasticity) — pinned default ------------------- #
# Pinned to 0.5 for ALL runs, overriding the SMM estimate (εY=1.48). The SMM
# value implied input SUBSTITUTABILITY (εY>1), which produced counterintuitive
# supply-shock responses (rising employment / output after adverse shocks).
# εY<1 = complements (gross-complementarity, the network-amplification regime).
# Override per-run with ENV["EPSY_OVERRIDE"], e.g. `EPSY_OVERRIDE=0.8 julia …`.
# Propagates to subprocesses launched by run_all_shocks.jl.
_epsY_set = haskey(ENV, "EPSY_OVERRIDE") ? parse(Float64, ENV["EPSY_OVERRIDE"]) : 0.5
modepsY   = fill(_epsY_set, nsec)
@printf "  εY pinned: εY = %.4f (overrides SMM estimate)\n" _epsY_set


# =========================================================================== #
#  THREE ε_Y CASES: NEAR-LEONTIEF / BASELINE / HIGH SUBSTITUTION              #
# =========================================================================== #

# Use baseline εY from SMM estimates (no εY sensitivity loop)
epsY_baseline = modepsY[1]
nT = 80   # IRF horizon (extended from 40 to assess convergence of the slow NFA mode)

@printf "\n%s\n  Baseline εY = %.4f (from SMM estimates)\n%s\n" repeat("─",60) epsY_baseline repeat("─",60)

# εY is already set in modepsY from the SMM parameter load — no override needed


# =========================================================================== #
#  STEADY STATE (same as main_SOE_gap.jl)                                      #
# =========================================================================== #

sigmaH     = sigmaH_val
varrho_val = modvarrho
gammag_vec = modgammag; gammas_vec = modgammas
om_g = ombar_val; om_s = 1 - ombar_val
chiX_vec = modchiX; omegaX = omegaX_val; etastar = etastar_val
Ystar = ystar_ss_val
alpha_vec = modalpha; alphaV_vec = modalphaV
beta_mat = modbeta; epsY_vec = modepsY; epsM_vec = modepsM; A_vec = modA

tb_target = let
    _tb = 0.02
    agg_mom_tb = joinpath(DATA_DIR, "aggregate_moments.csv")
    if isfile(agg_mom_tb)
        try
            agg_tb = CSV.read(agg_mom_tb, DataFrame)
            row    = filter(r -> String(r.moment) == "TBGDP", agg_tb)
            !isempty(row) && isfinite(row[1, :value]) && (_tb = Float64(row[1, :value]))
        catch; end
    end; _tb
end

x_guess = [ones(nsec); 1.0; 1.0; 1.0]
ss_result = nlsolve(
    (F, x) -> F .= steady_ntwsoe(
        x, PVstar_ss, epsilon, varrho_val, sigmaH,
        gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar,
        alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec,
        gamma, chi, psi, A_vec, tb_target
    ),
    x_guess; ftol=1e-14, show_trace=false, method=:trust_region)

!converged(ss_result) && @warn "SS solver did not converge (res=$(ss_result.residual_norm))"

pH_ss = ss_result.zero[1:nsec]
w_ss  = ss_result.zero[nsec+1]
Q_ss  = ss_result.zero[nsec+2]
C_ss  = ss_result.zero[nsec+3]

# Full SS evaluation (replicates main_SOE_gap.jl)
PL_ss  = fill(w_ss, nsec)
PV_ss  = Q_ss * PVstar_ss
MCi_ss = (epsilon-1)/epsilon .* pH_ss
PMi_ss = (beta_mat * (pH_ss .^ (1 .- epsM_vec))) .^ (1 ./ (1 .- epsM_vec))
P_ss   = (varrho_val .^ sigmaH .* pH_ss .^ (1-sigmaH)
         .+ (1 .- varrho_val) .^ sigmaH .* PV_ss .^ (1-sigmaH)) .^ (1/(1-sigmaH))
p_g_ss = prod(P_ss .^ gammag_vec); p_s_ss = prod(P_ss .^ gammas_vec)
C_g_ss = om_g * C_ss / p_g_ss; C_s_ss = om_s * C_ss / p_s_ss
C_gi_ss = gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss
C_si_ss = gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss
CHg_ss = varrho_val .^ sigmaH .* (pH_ss ./ P_ss) .^ (-sigmaH) .* C_gi_ss
CHs_ss = varrho_val .^ sigmaH .* (pH_ss ./ P_ss) .^ (-sigmaH) .* C_si_ss
CFg_ss = (1 .- varrho_val) .^ sigmaH .* (PV_ss ./ P_ss) .^ (-sigmaH) .* C_gi_ss
CFs_ss = (1 .- varrho_val) .^ sigmaH .* (PV_ss ./ P_ss) .^ (-sigmaH) .* C_si_ss
CHi_ss = CHg_ss .+ CHs_ss; CFi_ss = CFg_ss .+ CFs_ss
PX_ss  = prod(pH_ss .^ chiX_vec)
X_ss   = omegaX * (PX_ss / Q_ss)^(-etastar) * Ystar
Xi_ss  = chiX_vec .* X_ss .* PX_ss ./ pH_ss

ig_ss = zeros(nsec)
for i in 1:nsec, j in 1:nsec; ig_ss[i] += beta_mat[j, i]; end
ig_ss .*= mean(CHi_ss .+ Xi_ss)

M_init  = (MCi_ss ./ PMi_ss) .^ epsY_vec .* alpha_vec .* (CHi_ss .+ Xi_ss .+ ig_ss)
L_init  = (MCi_ss ./ PL_ss)  .^ epsY_vec .* (1 .- alpha_vec .- alphaV_vec) .* (CHi_ss .+ Xi_ss .+ ig_ss)
Vi_init = (MCi_ss ./ PV_ss)  .^ epsY_vec .* alphaV_vec .* (CHi_ss .+ Xi_ss .+ ig_ss)
Yi_init = A_vec .* (
    alpha_vec .^ (1 ./ epsY_vec) .* max.(M_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
  .+ alphaV_vec .^ (1 ./ epsY_vec) .* max.(Vi_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
  .+ (1 .- alphaV_vec .- alpha_vec) .^ (1 ./ epsY_vec) .* max.(L_init, 1e-20) .^ ((epsY_vec .- 1) ./ epsY_vec)
) .^ (epsY_vec ./ (epsY_vec .- 1))

inner_sol = nlsolve(
    (F, x) -> steady_ntwsoe_system!(F, x, alpha_vec, alphaV_vec, beta_mat,
        MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss, epsY_vec, epsM_vec, A_vec, pH_ss),
    [M_init; L_init; Vi_init; Yi_init]; ftol=1e-10, show_trace=false, method=:trust_region)

M_ss  = inner_sol.zero[1:nsec]
L_ss  = inner_sol.zero[nsec+1:2*nsec]
Vi_ss = inner_sol.zero[2*nsec+1:3*nsec]
Yi_ss = inner_sol.zero[3*nsec+1:4*nsec]

V_ss       = sum(Vi_ss); CF_ss = sum(CFi_ss)
mkupV      = 1.0; IMP_tot_ss = mkupV * (V_ss + CF_ss)
TB_ss      = PX_ss * X_ss - PV_ss * IMP_tot_ss
GDP_ss     = C_ss + TB_ss; N_ss = sum(L_ss)
r_star_ss  = Rworld_ss; Bstar_ss = -TB_ss / (Q_ss * (1 - r_star_ss / Pistar_ss))
bbar_val   = Q_ss * Bstar_ss / GDP_ss

PO_ss     = Q_ss * POstar_ss_val
PIV_ss_vec = [
    (modalphaOil[i] * PO_ss^(1-epsilonV_oil_val)
     + (1-modalphaOil[i]) * PV_ss^(1-epsilonV_oil_val))^(1/(1-epsilonV_oil_val))
    for i in 1:nsec
]

CFg_total_ss = sum(CFg_ss); CFs_total_ss = sum(CFs_ss)

@printf "  SS: GDP=%.4f  TB/GDP=%.4f  C=%.4f\n\n" GDP_ss (TB_ss/GDP_ss) C_ss


# =========================================================================== #
#  WRITE params_jl.mod AND RUN DYNARE                                          #
# =========================================================================== #

params_nt = (
    nsec=nsec, sigma_i_val=sigma_i_val, sigma_L_agg_val=sigma_L_agg_val,
    ilabcosts_val=ilabcosts_val, gamma_val=gamma, beta_val=beta_val,
    phi_val=phi_val, rho_val=rho_val, rho_om1_val=rho_om1_val,
    rho_tfp1_val=rho_tfp1_val, rho_tfp2_val=rho_tfp2_val,
    rhoi_val=rhoi_val, rhoirule_val=rhoirule_val, ombar_val=ombar_val,
    chii_b_val=chii_b_val, bbar_val=bbar_val,
    epsw_val=epsw_val, kappaw_val=kappaw_val,
    epsilonX_val=epsilonX_val, omegaX_val=omegaX_val,
    ystar_ss_val=ystar_ss_val, etastar_val=etastar_val,
    epsilonV_val=epsilonV_val, kappaV_val=kappaV_val, sigmaH_val=sigmaH_val,
    Pistar_ss_val=Pistar_ss, PVstar_ss_val=PVstar_ss,
    Rworld_ss_val=Rworld_ss,
    rho_pvstar_val=rho_pvstar_val, sigma_pvstar_val=sigma_pvstar_val,
    rho_xi_val=rho_xi_val, sigma_xi_val=sigma_xi_val,
    # Oil sector
    modalphaOil=modalphaOil, epsilonV_oil_val=epsilonV_oil_val,
    rho_postar_val=rho_postar_val, sigma_postar_val=sigma_postar_val,
    POstar_ss_val=POstar_ss_val, shock_eps_postar_val=shock_eps_postar_val,
    PIV_ss_vec=PIV_ss_vec,
    # Steady state
    w_ss=w_ss, C_ss=C_ss, GDP_ss=GDP_ss, N_ss=N_ss,
    p_s_ss=p_s_ss, p_g_ss=p_g_ss, C_s_ss=C_s_ss, C_g_ss=C_g_ss,
    Bstar_ss=Bstar_ss, Q_ss=Q_ss, TB_ss=TB_ss, PX_ss=PX_ss,
    V_ss=V_ss, CF_ss=CF_ss, CFg_total_ss=CFg_total_ss, CFs_total_ss=CFs_total_ss,
    IMP_ss_val=IMP_tot_ss, Ctot_ss_val=sum(gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss) + sum(gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss),
    Ctotg_ss_val=sum(gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss),
    Ctots_ss_val=sum(gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss),
    VA_ss_val=sum(pH_ss .* Yi_ss .- PMi_ss .* M_ss .- PV_ss .* Vi_ss), M_tot_ss=sum(M_ss), Y_tot_ss=sum(pH_ss .* Yi_ss),
    # Shock flags
    shock_eps_i_val=shock_eps_i_val, shock_eps_pvstar_val=shock_eps_pvstar_val,
    shock_eps_xi_val=shock_eps_xi_val,
    shock_epsA_val=shock_epsA_val,
    # Option-A
    sigma_om_vec=sigma_om_vec, shock_eps_om_vec=shock_eps_om_vec,
    # Sectoral vectors
    modgammag=modgammag, modgammas=modgammas, modalpha=modalpha,
    modalphaV=modalphaV, modepsY=modepsY, modepsM=modepsM,
    modkappa=modkappa, goods=goods, services=services,
    modcl=modcl, modclneg=modclneg, modcm=modcm,
    modchiX=modchiX, modvarrho=modvarrho,
    isigma_tfp_val=isigma_tfp_val, PL_ss=PL_ss,
    CFg_ss=CFg_ss, CFs_ss=CFs_ss, CHg_ss=CHg_ss, CHs_ss=CHs_ss,
    Vi_ss=Vi_ss, pH_ss=pH_ss, MCi_ss=MCi_ss, Yi_ss=Yi_ss,
    L_ss=L_ss, C_gi_ss=C_gi_ss, C_si_ss=C_si_ss, P_ss=P_ss,
    PMi_ss=PMi_ss, M_ss=M_ss, modbeta=modbeta,
)

params_mod_path = write_params_mod(MOD_DIR, params_nt)
@printf "--- Wrote params_jl.mod ---\n\n"

# Run Dynare subprocess
dynare_script = joinpath(SCRIPT_DIR, "run_dynare_subprocess.jl")
julia_exe     = joinpath(Sys.BINDIR, "julia")
project_dir   = dirname(Base.active_project())

@printf "--- Running Dynare.jl (subprocess) ---\n"
dynare_out = IOBuffer()
dynare_cmd = ignorestatus(`$julia_exe --project=$project_dir $dynare_script $MOD_DIR`)
dynare_proc = run(pipeline(dynare_cmd, stdout=dynare_out, stderr=stderr), wait=true)

dynare_stdout = String(take!(dynare_out))
print(dynare_stdout)

if !success(dynare_proc)
    error("Dynare subprocess crashed (exit code $(dynare_proc.exitcode)). See stderr above.")
end

if !occursin("DYNARE_SUCCESS", dynare_stdout)
    error("Dynare subprocess did not report success. See output above.")
end

@printf "\n--- Dynare completed ---\n\n"


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
alpha_L_vec = 1.0 .- modalpha .- modalphaV
ge_kwargs = (po_irf=po_irf, pv_irf=pv_irf, modalphaOil=modalphaOil)
ge_h1 = ge_mc_components(:price, 1, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; ge_kwargs...)
ge_h2 = ge_mc_components(:price, 2, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; ge_kwargs...)
ge_h4 = ge_mc_components(:price, 4, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; ge_kwargs...)
ge_fun = h -> ge_mc_components(:price, h, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha, alpha_L_vec, modbeta; ge_kwargs...)
ge_colors, ge_labels = ge_component_style(:price)

# "Affected" = sectors with above-average DIRECT oil cost share (αV·αOil).
# Used for the affected-vs-rest inflation figure. The oil shock hits every
# sector, so the inflation decompositions use a SINGLE axis (no pi_decomp_axis2).
oil_affected = findall((modalphaV .* modalphaOil) .> mean(modalphaV .* modalphaOil))

@printf "--- Aggregate IRFs (10%% oil shock, impact) ---\n"
@printf "  GDP %+.3f%%   pi %+.3f ann.pp   Q %+.3f%%   r %+.3f ann.pp\n\n" gdp_irf[1] (pi_irf[1]*4) q_irf[1] (r_irf[1]*4)

# Assemble context and generate the full figure + table set via the shared module
ctx = (
    # Tag carries the chii_b override so sensitivity runs don't overwrite the
    # baseline figures/tables (e.g. tag = "oil_chiib0p01").
    HAS_PLOTS = _HAS_PLOTS[],
    tag = "oil" *
        (haskey(ENV, "CHIIB_OVERRIDE")  ? "_chiib" * replace(ENV["CHIIB_OVERRIDE"],  "." => "p") : "") *
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
