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
    irf_aggregate_oil_shock.pdf     — GDP, inflation, RER, trade balance
    irf_sectoral_Y_oil_shock.pdf    — 12-sector output IRFs
    irf_sectoral_PH_oil_shock.pdf   — 12-sector price IRFs
    irf_sectoral_MC_oil_shock.pdf   — 12-sector marginal cost IRFs
    decomposition_output_oil.pdf        — sectoral output and MC at impact
    decomposition_mc_inflation_oil.pdf     — MC decomp (direct+network, stacked) + impact inflation diamonds
    decomposition_mc_inflation_6m_oil.pdf  — same bars + 6-month cumulative inflation diamonds
    decomposition_mc_inflation_12m_oil.pdf — same bars + 12-month cumulative inflation diamonds
    oil_intensity_exposure.pdf          — sector oil exposure map
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
using StatsBase

# =========================================================================== #
#  INCLUDE HELPERS                                                              #
# =========================================================================== #

SCRIPT_DIR = @__DIR__

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))

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

beta_val = 0.986;  epsilon = 10.0;  gamma = 2.0;  psi = 1.0;  chi = 1.0
modkappa = theta_vec .* (epsilon-1) ./ ((1 .- theta_vec) .* (1 .- theta_vec .* beta_val))

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
chii_b_val = 0.001; etastar_val = 3.5
xi_rstar_val = 0.2; ystar_ss_val = 1.0; PVstar_ss = 1.0; sigmaH_val = 0.999

modchiX   = fill(1/nsec, nsec)
modvarrho = var_rho
modA      = ones(nsec)

# ---- Oil price shock parameters ----
modalphaOil = [0.1635, 0.2164, 0.1890, 0.0871, 0.0417,
               0.0793, 0.3734, 0.0043, 0.0449, 0.0633, 0.0391, 0.0447]
epsilonV_oil_val  = 0.5
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
smm_est_file = joinpath(DATA_DIR, "computed", "smm_estimates.csv")
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


# =========================================================================== #
#  THREE ε_Y CASES: NEAR-LEONTIEF / BASELINE / HIGH SUBSTITUTION              #
# =========================================================================== #

# Use baseline εY from SMM estimates (no εY sensitivity loop)
epsY_baseline = modepsY[1]
nT = 40   # IRF horizon

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
    VA_ss_val=sum(Yi_ss .- M_ss), M_tot_ss=sum(M_ss), Y_tot_ss=sum(Yi_ss),
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

# Identify eps_postar column
# varexo order: eps_om eps_i epschi eps_pvstar eps_postar epsA_1..12 eps_xi
# => eps_postar should be column 5
postar_col = 5
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

n_irf = 40   # quarters

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
#  DECOMPOSITION: DIRECT OIL COST vs. NETWORK AMPLIFICATION                   #
#                                                                              #
#  We use the Leontief inverse to separate the first-round cost push (direct  #
#  oil price increase on each sector's import costs) from the IO-network      #
#  amplification (cost increases propagated through intermediate input trade). #
#                                                                              #
#  Direct effect on sector i's marginal cost:                                  #
#    d_direct_i ≈ alphaV_i × alphaOilShare_i × d(log PO)                     #
#  This is the cost increase from oil imports alone, ignoring IO linkages.     #
#                                                                              #
#  Total effect includes network propagation via the Leontief inverse:         #
#    d_total_i = [(I - diag(alpha) × Gamma)^{-1} × d_direct]_i               #
#  where Gamma is the IO matrix and alpha is the materials share vector.       #
#                                                                              #
#  Network amplification = total - direct                                      #
# =========================================================================== #

@printf "--- Decomposing direct vs. network effects ---\n"

# Direct cost exposure: share of total cost that is oil
# In the CES production function, the cost share of imported inputs in MC is
# approximately alphaV_i (calibrated share), and of that, alphaOilShare_i is oil.
# First-order log-linear approximation:
# d(log MC_i) ≈ s_oil_i × d(log PO) + s_M_i × d(log PM_i) + s_L_i × d(log w)
# where s_oil_i = alphaV_i × alphaOilShare_i / MC normalization

direct_oil_exposure = modalphaV .* modalphaOil   # direct oil cost share per sector

# Leontief inverse for network amplification
# The IO matrix Gamma[i,j] = share of sector i's material inputs from j
# Materials cost share in production: alpha_i
# The amplification matrix is (I - diag(alpha) × Gamma)^{-1}
Gamma_io = modbeta   # purchaser × supplier
Leontief = inv(I(nsec) - Diagonal(modalpha) * Gamma_io)

# Direct effect vector (on MC, in % for 10% oil shock)
d_log_PO = log(1 + shock_pct) * 100   # ~9.53%
direct_mc = direct_oil_exposure .* d_log_PO   # direct cost push per sector

# Total effect via network (first-order approximation)
total_mc_approx = Leontief * direct_mc
network_mc = total_mc_approx .- direct_mc

# Amplification ratio
amp_ratio = total_mc_approx ./ max.(direct_mc, 1e-10)

# Also extract model-based IRFs for comparison
mc_irf_impact = [get_irf_var("MC_$(i)")[1] for i in 1:nsec]
y_irf_impact  = [get_irf_var("Y_$(i)")[1] for i in 1:nsec]
ph_irf_impact = [get_irf_var("PH_$(i)")[1] for i in 1:nsec]

# Aggregate effects (output-weighted)
yi_weights = Yi_ss ./ sum(Yi_ss)
agg_direct  = sum(yi_weights .* direct_mc)
agg_total   = sum(yi_weights .* total_mc_approx)
agg_network = agg_total - agg_direct
agg_amp     = agg_total / agg_direct

@printf "\n  %-45s  %8s  %8s  %8s  %8s\n" "Sector" "Direct" "Network" "Total" "Amplif."
@printf "  %s\n" repeat("-", 85)
for i in 1:nsec
    @printf "  %-45s  %7.3f%%  %7.3f%%  %7.3f%%  %7.2fx\n" names_vec[i] direct_mc[i] network_mc[i] total_mc_approx[i] amp_ratio[i]
end
@printf "  %s\n" repeat("-", 85)
@printf "  %-45s  %7.3f%%  %7.3f%%  %7.3f%%  %7.2fx\n\n" "Aggregate (output-weighted)" agg_direct agg_network agg_total agg_amp


# =========================================================================== #
#  AGGREGATE IRF RESULTS                                                       #
# =========================================================================== #

gdp_irf = get_irf_var("GDP")
pi_irf  = get_irf_var("pi")
q_irf   = get_irf_var("Q")
tb_irf  = get_irf_var("TB")
c_irf   = get_irf_var("C")
cg_irf  = get_irf_var("C_g")   # goods consumption
cs_irf  = get_irf_var("C_s")   # services consumption
n_irf_v = get_irf_var("N")
r_irf   = get_irf_var("r")
w_irf   = get_irf_var("w")
po_irf  = get_irf_var("PO")    # domestic oil price (includes exchange rate)
pv_irf  = get_irf_var("PV")    # composite non-oil import price (exchange-rate driven)
alpha_L = 1.0 .- modalpha .- modalphaV   # labor cost share per sector

@printf "--- Aggregate IRFs (10%% oil shock, impact period) ---\n"
@printf "  GDP    : %+.3f%% dev.\n" gdp_irf[1]
@printf "  CPI (π): %+.3f ann. pp\n" (pi_irf[1]*4)
@printf "  RER (Q): %+.3f%% dev.\n" q_irf[1]
@printf "  TB     : %+.4f pp of GDP\n" (tb_irf[1]*TB_ss/GDP_ss)
@printf "  C      : %+.3f%% dev.\n" c_irf[1]
@printf "  N      : %+.3f%% dev.\n" n_irf_v[1]
@printf "  r      : %+.3f ann. pp\n" (r_irf[1]*4)
@printf "\n"

# Peak effects
gdp_trough = minimum(gdp_irf)
gdp_trough_q = argmin(gdp_irf)
pi_peak = maximum(pi_irf)
pi_peak_q = argmax(pi_irf)

@printf "  GDP trough : %+.3f%% at quarter %d\n" gdp_trough gdp_trough_q
@printf "  π peak     : %+.3f%% at quarter %d\n" pi_peak pi_peak_q

# ---- Compute additional matrices for this case ----
ph_irf_mat = hcat([get_irf_var("PH_$(i)") for i in 1:nsec]...)'  # nsec × n_irf
mc_irf_mat = hcat([get_irf_var("MC_$(i)") for i in 1:nsec]...)'
l_irf_mat  = hcat([get_irf_var("L_$(i)")  for i in 1:nsec]...)'
ygap_irf_mat = hcat([get_irf_var("Ygap_$(i)") for i in 1:nsec]...)'  # nsec × n_irf
ygap_irf   = get_irf_var("Ygap")     # aggregate output gap
gdpgap_irf = get_irf_var("GDPgap")   # GDP gap

# Sectoral nominal home-price inflation (annualized)
pi_sec_mat = zeros(nsec, n_irf)
for i in 1:nsec
    pi_sec_mat[i, 1] = (ph_irf_mat[i, 1] + pi_irf[1]) * 4
    for h in 2:n_irf
        pi_sec_mat[i, h] = (ph_irf_mat[i, h] - ph_irf_mat[i, h-1] + pi_irf[h]) * 4
    end
end

# Consumption-weighted inflation aggregates
cons_ss_g = C_gi_ss; cons_ss_s = C_si_ss
cons_ss_all = cons_ss_g .+ cons_ss_s; tot_cons_all = sum(cons_ss_all)
w_agg_c = cons_ss_all ./ tot_cons_all
w_g_c   = cons_ss_g ./ sum(cons_ss_g)
w_s_c   = cons_ss_s ./ sum(cons_ss_s)
pi_agg_irf   = [sum(w_agg_c .* pi_sec_mat[:, h]) for h in 1:n_irf]
pi_goods_irf = [sum(w_g_c   .* pi_sec_mat[:, h]) for h in 1:n_irf]
pi_serv_irf  = [sum(w_s_c   .* pi_sec_mat[:, h]) for h in 1:n_irf]

# Inflation at various horizons
infl_irf_impact = ph_irf_impact .+ pi_irf[1]
infl_6m  = ph_irf_mat[:, 2] .+ sum(pi_irf[1:2])
infl_12m = ph_irf_mat[:, 4] .+ sum(pi_irf[1:4])

# GE decomposition components at horizons 1, 2, 4
alpha_L_vec = 1.0 .- modalpha .- modalphaV
function _mc_decomp(h, _ph, _mc, _po, _pv, _w, _aV, _aO, _am, _aL, _beta)
    dir_h   = _aV .* _aO .* _po[h]
    net_h   = _am .* (_beta * _ph[:, h])
    noi_h   = _aV .* (1.0 .- _aO) .* _pv[h]
    lab_h   = _aL .* _w[h]
    res_h   = _mc[:, h] .- dir_h .- net_h .- noi_h .- lab_h
    return (dir=copy(dir_h), net=copy(net_h), noi=copy(noi_h),
            lab=copy(lab_h), res=copy(res_h), mc=copy(_mc[:, h]))
end
ge_h1 = _mc_decomp(1, ph_irf_mat, mc_irf_mat, po_irf, pv_irf, w_irf,
                    modalphaV, modalphaOil, modalpha, alpha_L_vec, modbeta)
ge_h2 = _mc_decomp(2, ph_irf_mat, mc_irf_mat, po_irf, pv_irf, w_irf,
                    modalphaV, modalphaOil, modalpha, alpha_L_vec, modbeta)
ge_h4 = _mc_decomp(4, ph_irf_mat, mc_irf_mat, po_irf, pv_irf, w_irf,
                    modalphaV, modalphaOil, modalpha, alpha_L_vec, modbeta)

# Helper: retrieve IRF for a named variable (used by plotting code)
get_irf(vn::String) = begin
    i = get(endo_idx, vn, 0)
    i == 0 ? zeros(n_irf) : irf_pct[i, :]
end

@printf "\n%s\n  Baseline case computed.\n%s\n" repeat("=",60) repeat("=",60)


# =========================================================================== #
#  GENERATE FIGURES                                                            #
# =========================================================================== #

# Suppress PlotUtils "No strict ticks" Julia warnings and GR/Qt C-level stderr
# (QPainterPath NaN messages) that fire during PDF rendering.
function savefig_quiet(p, path)
    try
        redirect_stderr(devnull) do
            with_logger(NullLogger()) do
                Plots.savefig(p, path)
            end
        end
    catch e
        @warn "savefig failed for $path" exception=(e, catch_backtrace())
        # Retry without stderr suppression so GR errors are visible
        Plots.savefig(p, path)
    end
end

# Save figure locally and sync to Overleaf Dropbox folder if available.
# For Overleaf, write to a temp file first then copy (avoids Dropbox lock issues).
function save_fig(p, fname)
    local_path = joinpath(FIGURES_DIR, fname)
    savefig_quiet(p, local_path)
    @printf "  Saved: %s\n" fname
    if overleaf_ok
        ol_path = joinpath(OVERLEAF_FIG_DIR, fname)
        # Copy from the already-saved local file instead of re-rendering
        cp(local_path, ol_path; force=true)
        @printf "  → Overleaf: %s\n" ol_path
    end
end

# Copy a table file (tex/csv) to the Overleaf tables folder if available.
function sync_table(fname)
    overleaf_ok || return
    src = joinpath(TABLES_DIR, fname)
    dst = joinpath(OVERLEAF_TAB_DIR, fname)
    cp(src, dst; force=true)
    @printf "  → Overleaf: %s\n" dst
end

if _HAS_PLOTS[]
    @printf "\n--- Generating figures (baseline εY = %.2f) ---\n" epsY_baseline
    periods = 1:nT

    short_names = [s[1:min(18,length(s))] for s in names_vec]
    sector_colors = vcat(fill(:steelblue, 5), fill(:firebrick, 7))
    bar_names = [
        "Agric. & Fishing", "Mining", "Manufacturing", "Utilities",
        "Construction", "Trade & Hotels", "Transport & ICT", "Finance",
        "Real Estate", "Business Serv.", "Personal Serv.", "Public Admin.",
    ]

    baseline_label = "Baseline (εY=$(round(epsY_baseline,digits=2)))"
    baseline_color = :steelblue
    baseline_lw    = 2.5

    # ================================================================== #
    #  Figure 1: Aggregate IRFs (2×5)                                    #
    # ================================================================== #
    agg_panels = [
        (gdp_irf, "GDP", "% dev.", 1.0),
        (pi_irf,  "CPI Inflation", "ann. pp", 4.0),
        (q_irf,   "Real Exchange Rate", "% dev.", 1.0),
        (tb_irf,  "Trade Balance", "pp of GDP", TB_ss/GDP_ss),
        (cg_irf,  "Goods Consumption", "% dev.", 1.0),
        (r_irf,   "Nominal Interest Rate", "ann. pp", 4.0),
        (w_irf,   "Real Wage", "% dev.", 1.0),
        (pi_goods_irf,  "Goods Inflation", "ann. pp", 1.0),
        (pi_serv_irf,   "Services Inflation", "ann. pp", 1.0),
        (cs_irf,  "Services Consumption", "% dev.", 1.0),
    ]
    p_agg = Plots.plot(layout=(2,5), size=(2000,700),
        plot_title="10% Oil Price Shock — Aggregate Responses",
        titlefontsize=12, margin=5Plots.mm)
    for (k, (irf_v, ttl, yl, scl)) in enumerate(agg_panels)
        Plots.plot!(p_agg, periods, irf_v .* scl, subplot=k,
            label=(k==1 ? baseline_label : ""), color=baseline_color, lw=baseline_lw,
            title=ttl, ylabel=(k∈[1,6] ? yl : ""),
            xlabel=(k>5 ? "Quarters" : ""))
        Plots.hline!(p_agg, [0.0], subplot=k, color=:black, lw=0.6, ls=:dash, label="")
    end
    save_fig(p_agg, "irf_aggregate_oil_shock.pdf")


    # ================================================================== #
    #  Figure 2: Sectoral Output IRFs (4×3)                               #
    # ================================================================== #
    p_sec_y = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Output",
        titlefontsize=10)
    for i in 1:12
        Plots.plot!(p_sec_y, periods, get_irf("Y_$(i)"), subplot=i,
            label="", color=baseline_color, lw=baseline_lw,
            title=short_names[i], titlefontsize=9)
        Plots.hline!(p_sec_y, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_y, "irf_sectoral_Y_oil_shock.pdf")


    # ================================================================== #
    #  Figure 3: Sectoral Price IRFs (4×3)                                #
    # ================================================================== #
    p_sec_ph = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Home Prices",
        titlefontsize=10)
    for i in 1:12
        Plots.plot!(p_sec_ph, periods, get_irf("PH_$(i)"), subplot=i,
            label="", color=baseline_color, lw=baseline_lw,
            title=short_names[i], titlefontsize=9)
        Plots.hline!(p_sec_ph, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_ph, "irf_sectoral_PH_oil_shock.pdf")


    # ================================================================== #
    #  Figure 5b: Sectoral Inflation IRFs (4×3)                          #
    # ================================================================== #
    p_sec_pi = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Home-Price Inflation",
        titlefontsize=10)
    for i in 1:12
        Plots.plot!(p_sec_pi, periods, pi_sec_mat[i, :], subplot=i,
            label="", color=baseline_color, lw=baseline_lw,
            title=short_names[i], titlefontsize=9)
        Plots.hline!(p_sec_pi, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_pi, "irf_sectoral_inflation_oil.pdf")


    # ================================================================== #
    #  Figure 5c: Aggregate Inflation by Group                           #
    # ================================================================== #
    p_pi_agg = Plots.plot(size=(1000, 500),
        title="10% Oil Shock — Aggregate Home-Price Inflation",
        titlefontsize=11, xlabel="Quarters", ylabel="Ann. pp dev. from SS",
        legend=:topright, margin=5Plots.mm)
    Plots.plot!(p_pi_agg, periods, pi_agg_irf,
        label="Aggregate", color=baseline_color, lw=baseline_lw)
    Plots.hline!(p_pi_agg, [0.0], color=:black, lw=0.5, ls=:dash, label="")
    save_fig(p_pi_agg, "irf_aggregate_inflation_oil.pdf")


    # ================================================================== #
    #  Figure 5d: Goods vs Services Inflation                            #
    # ================================================================== #
    p_gs_pi = Plots.plot(size=(800, 500),
        title="10% Oil Shock — Home-Price Inflation: Goods vs Services",
        titlefontsize=11, xlabel="Quarters", ylabel="Ann. pp dev. from SS",
        legend=:topright, margin=5Plots.mm)
    Plots.plot!(p_gs_pi, periods, pi_agg_irf,
        label="Aggregate", color=:black, lw=2.5, ls=:solid)
    Plots.plot!(p_gs_pi, periods, pi_goods_irf,
        label="Goods (sectors 1–5)", color=:steelblue, lw=2, ls=:dash)
    Plots.plot!(p_gs_pi, periods, pi_serv_irf,
        label="Services (sectors 6–12)", color=:firebrick, lw=2, ls=:dot)
    Plots.hline!(p_gs_pi, [0.0], color=:black, lw=0.5, ls=:dash, label="")
    save_fig(p_gs_pi, "irf_goods_vs_services_inflation_oil.pdf")


    # ================================================================== #
    #  Figure 7: Sectoral MC IRFs (4×3)                                  #
    # ================================================================== #
    p_sec_mc = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Marginal Cost",
        titlefontsize=10)
    for i in 1:12
        Plots.plot!(p_sec_mc, periods, mc_irf_mat[i, :], subplot=i,
            label="", color=baseline_color, lw=baseline_lw,
            title=short_names[i], titlefontsize=9)
        Plots.hline!(p_sec_mc, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_mc, "irf_sectoral_MC_oil_shock.pdf")


    # ================================================================== #
    #  Figure 8: Aggregate Labor Market (1×2)                            #
    # ================================================================== #
    p_lab_agg = Plots.plot(layout=(1,2), size=(1000, 400),
        plot_title="10% Oil Shock — Aggregate Labor Market",
        titlefontsize=11, margin=5Plots.mm)
    Plots.plot!(p_lab_agg, periods, n_irf_v, subplot=1,
        label=baseline_label, color=baseline_color, lw=baseline_lw,
        xlabel="Quarters", ylabel="% dev. from SS", title="Aggregate Employment")
    Plots.plot!(p_lab_agg, periods, w_irf, subplot=2,
        label=baseline_label, color=baseline_color, lw=baseline_lw,
        xlabel="Quarters", ylabel="% dev. from SS", title="Real Wage")
    Plots.hline!(p_lab_agg, [0.0], subplot=1, color=:black, lw=0.6, ls=:dash, label="")
    Plots.hline!(p_lab_agg, [0.0], subplot=2, color=:black, lw=0.6, ls=:dash, label="")
    save_fig(p_lab_agg, "irf_labor_aggregate_oil_shock.pdf")


    # ================================================================== #
    #  Figure 9: Sectoral Employment (4×3)                               #
    # ================================================================== #
    p_sec_l = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Employment",
        titlefontsize=10)
    for i in 1:12
        Plots.plot!(p_sec_l, periods, l_irf_mat[i, :], subplot=i,
            label="", color=baseline_color, lw=baseline_lw,
            title=short_names[i], titlefontsize=9)
        Plots.hline!(p_sec_l, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_l, "irf_sectoral_L_oil_shock.pdf")


    # ================================================================== #
    #  Figure 10: GDP Gap                                                #
    # ================================================================== #
    p_gap_agg = Plots.plot(size=(700, 400),
        title="10% Oil Shock — GDP Gap",
        titlefontsize=11, margin=5Plots.mm,
        xlabel="Quarters", ylabel="% dev. from SS")
    Plots.plot!(p_gap_agg, periods, gdpgap_irf,
        label=baseline_label, color=baseline_color, lw=baseline_lw)
    Plots.hline!(p_gap_agg, [0.0], color=:black, lw=0.6, ls=:dash, label="")
    save_fig(p_gap_agg, "irf_gdpgap_aggregate_oil_shock.pdf")


    # ================================================================== #
    #  Figure 11: Sectoral Output Gaps (4×3)                             #
    # ================================================================== #
    p_sec_ygap = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Output Gap",
        titlefontsize=10)
    for i in 1:12
        Plots.plot!(p_sec_ygap, periods, ygap_irf_mat[i, :], subplot=i,
            label="", color=baseline_color, lw=baseline_lw,
            title=short_names[i], titlefontsize=9)
        Plots.hline!(p_sec_ygap, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_ygap, "irf_sectoral_Ygap_oil_shock.pdf")


    # ================================================================== #
    #  Figure 5: Oil Intensity Map (no εY dependence)                    #
    # ================================================================== #
    direct_oil_exposure = modalphaV .* modalphaOil   # oil share of total inputs
    p_oil = Plots.bar(1:12, direct_oil_exposure .* 100,
        xticks=(1:12, [s[1:min(12,length(s))] for s in names_vec]),
        xrotation=45, label="Oil share of total inputs (%)",
        color=:darkorange, ylabel="Oil cost share (%)",
        title="Oil Share of Total Inputs by Sector",
        titlefontsize=14, size=(900, 450), bottom_margin=10Plots.mm)
    save_fig(p_oil, "oil_intensity_exposure.pdf")


    # ================================================================== #
    #  Decomposition Figures — ONE PER εY CASE                           #
    #                                                                     #
    #  GE MC decomposition + inflation pass-through (stacked bars).      #
    #  Generates separate PDF for each εY value at each horizon.         #
    # ================================================================== #

    # Helper: solid rectangle shape for manual bar-chart plotting
    bar_rect(x, y0, y1, w=0.65) = Plots.Shape(
        [x - w/2, x + w/2, x + w/2, x - w/2],
        [y0,      y0,      y1,      y1      ])

    ge_comp_colors = [:darkorange, :steelblue, :forestgreen, :crimson, :gray60]
    ge_comp_labels = ["Direct oil (pp)", "Network via prices (pp)",
                      "Non-oil imports (pp)", "Labor (pp)", "Residual (pp)"]

    function make_ge_decomp_fig(dir_h, net_h, noimp_h, lab_h, resid_h, mc_h,
                                 infl_h, infl_label, title_str, ylabel_str)
        comps   = hcat(dir_h, net_h, noimp_h, lab_h, resid_h)
        pos_mat = max.(comps, 0.0)
        neg_mat = min.(comps, 0.0)
        pos_tops = vec(sum(pos_mat, dims=2))
        neg_tops = vec(sum(neg_mat, dims=2))
        all_vals = vcat(pos_tops, neg_tops, infl_h, [0.0])
        ylo = min(minimum(all_vals), 0.0) * 1.40
        yhi = max(maximum(all_vals), 0.0) * 1.35
        p = Plots.plot(
            xticks=(1:nsec, bar_names), xrotation=55,
            ylabel=ylabel_str, title=title_str, titlefontsize=12,
            size=(1400, 660), legend=:topright, ylims=(ylo, yhi),
            bottom_margin=26Plots.mm, left_margin=14Plots.mm,
            right_margin=8Plots.mm, top_margin=3Plots.mm,
            xlims=(0.3, nsec + 0.7))
        Plots.hline!(p, [0.0], color=:black, lw=0.6, ls=:dash, label="")
        ncomp = 5; bw = 0.65
        for k in 1:ncomp
            pos_bot = k > 1 ? vec(sum(pos_mat[:, 1:k-1], dims=2)) : zeros(nsec)
            neg_bot = k > 1 ? vec(sum(neg_mat[:, 1:k-1], dims=2)) : zeros(nsec)
            lbl_used = false
            for i in 1:nsec
                lbl = lbl_used ? "" : ge_comp_labels[k]
                if pos_mat[i, k] > 1e-10
                    s = bar_rect(i, pos_bot[i], pos_bot[i] + pos_mat[i, k], bw)
                    Plots.plot!(p, s, color=ge_comp_colors[k], label=lbl,
                        alpha=0.85, linecolor=:white, linewidth=0.3)
                    lbl_used = true; lbl = ""
                end
                if abs(neg_mat[i, k]) > 1e-10
                    s = bar_rect(i, neg_bot[i] + neg_mat[i, k], neg_bot[i], bw)
                    Plots.plot!(p, s, color=ge_comp_colors[k], label=lbl,
                        alpha=0.85, linecolor=:white, linewidth=0.3)
                    lbl_used = true
                end
            end
            if !lbl_used
                Plots.plot!(p, [NaN], [NaN], color=ge_comp_colors[k],
                    label=ge_comp_labels[k], lw=4)
            end
        end
        Plots.scatter!(p, 1:nsec, infl_h,
            color=:black, markershape=:diamond,
            markersize=7, markerstrokewidth=1, label=infl_label)
        return p
    end

    # Generate decomposition figures (baseline only)
    epsY_str = "$(round(epsY_baseline, digits=2))"

    # Impact (h = 1)
    p1 = make_ge_decomp_fig(
        ge_h1.dir, ge_h1.net, ge_h1.noi, ge_h1.lab, ge_h1.res, ge_h1.mc,
        infl_irf_impact,
        "Home-price inflation (q-o-q pp, t=1)",
        "Oil Shock — GE MC Decomposition (Impact, εY=$(epsY_str))",
        "Percentage points (t = 1)")
    save_fig(p1, "decomposition_mc_inflation_oil_baseline.pdf")

    # 6-month (h = 2)
    p2 = make_ge_decomp_fig(
        ge_h2.dir, ge_h2.net, ge_h2.noi, ge_h2.lab, ge_h2.res, ge_h2.mc,
        infl_6m,
        "Home-price inflation (6-month cumulative, pp)",
        "Oil Shock — GE MC Decomposition (6-Month, εY=$(epsY_str))",
        "Percentage points")
    save_fig(p2, "decomposition_mc_inflation_6m_oil_baseline.pdf")

    # 12-month (h = 4)
    p4 = make_ge_decomp_fig(
        ge_h4.dir, ge_h4.net, ge_h4.noi, ge_h4.lab, ge_h4.res, ge_h4.mc,
        infl_12m,
        "Home-price inflation (12-month cumulative, pp)",
        "Oil Shock — GE MC Decomposition (12-Month, εY=$(epsY_str))",
        "Percentage points")
    save_fig(p4, "decomposition_mc_inflation_12m_oil_baseline.pdf")

    @printf "  Decomposition figures saved (baseline)\n"

    # ---- Output + MC decomposition (impact) ---- #
    p_out = groupedbar(
        [y_irf_impact mc_irf_impact],
        xticks=(1:12, bar_names), xrotation=55,
        label=["Output Y_i (% dev.)" "MC_i (% dev.)"],
        color=[:steelblue :firebrick],
        ylabel="% change from SS (t = 1)",
        title="Oil Shock Impact: Output & MC (εY=$(epsY_str))",
        titlefontsize=12, size=(1400, 600), legend=:topright,
        bottom_margin=24Plots.mm, left_margin=14Plots.mm, right_margin=5Plots.mm)
    save_fig(p_out, "decomposition_output_oil_baseline.pdf")


else
    @printf "\n  [Plots.jl not available — skipping figures]\n"
end


# =========================================================================== #
#  GENERATE LaTeX TABLES                                                       #
# =========================================================================== #

@printf "\n--- Generating LaTeX tables (baseline εY = %.2f) ---\n" epsY_baseline

# All variables are already in scope from the single baseline computation

# ---- Table 1: Sectoral Decomposition ---- #
# English abbreviated sector names (consistent with bar_names used in figures)
tab_names = [
    "Agric.\\ \\& Fishing",
    "Mining",
    "Manufacturing",
    "Utilities",
    "Construction",
    "Trade \\& Hotels",
    "Transport \\& ICT",
    "Finance",
    "Real Estate",
    "Business Serv.",
    "Personal Serv.",
    "Public Admin.",
]

open(joinpath(TABLES_DIR, "oil_shock_decomposition.tex"), "w") do f
    println(f, "\\begin{table}[htbp]")
    println(f, "\\centering")
    println(f, "\\caption{Sectoral Decomposition of a 10\\% Oil Price Shock}")
    println(f, "\\label{tab:oil_decomp}")
    println(f, "\\begin{threeparttable}")
    println(f, "\\small")
    println(f, "\\begin{tabular}{@{}l c c c c c c@{}}")
    println(f, "\\toprule")
    println(f, " & Oil & \\multicolumn{3}{c}{MC Increase (\\%)} & Amplif. & Output \\\\")
    println(f, "\\cmidrule(lr){3-5}")
    println(f, "Sector & Share (\\%) & Direct & Network & Total & Ratio & IRF (\\%) \\\\")
    println(f, "\\midrule")

    for i in 1:nsec
        gs = goods[i] ? "G" : "S"
        sign_str = y_irf_impact[i] < 0 ? "\$-\$$(abs(round(y_irf_impact[i], digits=3)))" :
                                          @sprintf("%.3f", y_irf_impact[i])
        @printf(f, "%s (%s) & %4.1f & %.3f & %.3f & %.3f & %5.2f & \$%+.3f\$ \\\\\n",
            tab_names[i], gs,
            modalphaOil[i]*100, direct_mc[i], network_mc[i],
            total_mc_approx[i], amp_ratio[i], y_irf_impact[i])
    end

    println(f, "\\midrule")
    @printf(f, "Aggregate (output-weighted) & --- & %.3f & %.3f & %.3f & %5.2f & \$%+.3f\$ \\\\\n",
        agg_direct, agg_network, agg_total, agg_amp, gdp_irf[1])
    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\begin{tablenotes}[flushleft]")
    println(f, "\\footnotesize")
    println(f, "\\item \\textit{Notes:} Oil Share is \$\\alpha^{\\text{Oil}}_i\$, the fraction of sector \$i\$'s intermediate imports that are oil/refined fuels (Chile 2021 IO tables). Direct MC is \$\\alpha_{Vi}\\,\\alpha^{\\text{Oil}}_i\\,\\Delta\\log P^O\$ (equation~\\eqref{eq:direct_mc}). Network is additional cost via IO linkages under the Leontief PE approximation (equation~\\eqref{eq:network_mc}). Total = Direct + Network (equation~\\eqref{eq:leontief_mc}). Amplification ratio = Total/Direct. Output IRF is the impact-period sectoral output response from the linearized model. G\\,=\\,Goods, S\\,=\\,Services.")
    println(f, "\\end{tablenotes}")
    println(f, "\\end{threeparttable}")
    println(f, "\\end{table}")
end
@printf "  Saved: oil_shock_decomposition.tex\n"
sync_table("oil_shock_decomposition.tex")


# ---- Table 2: Aggregate Responses ---- #
open(joinpath(TABLES_DIR, "oil_shock_aggregates.tex"), "w") do f
    println(f, "\\begin{table}[htbp]")
    println(f, "\\centering")
    println(f, "\\caption{Aggregate Responses to a 10\\% Oil Price Shock}")
    println(f, "\\label{tab:oil_agg}")
    println(f, "\\begin{threeparttable}")
    println(f, "\\begin{tabular}{@{}l c c c@{}}")
    println(f, "\\toprule")
    println(f, "Variable & Impact & Trough/Peak & Quarter \\\\")
    println(f, "\\midrule")

    vars_agg = [
        ("GDP",            "GDP",  gdp_irf,              "\\% dev."),
        ("CPI Inflation",  "pi",   pi_irf .* 4,          "ann.\\ pp"),
        ("Consumption",    "C",    c_irf,                 "\\% dev."),
        ("Employment",     "N",    n_irf_v,               "\\% dev."),
        ("Real Exch.\\ Rate","Q",  q_irf,                 "\\% dev."),
        ("Trade Balance",  "TB",   tb_irf .* (TB_ss/GDP_ss), "pp of GDP"),
        ("Policy Rate",    "r",    r_irf .* 4,            "ann.\\ pp"),
    ]
    for (label, _, irf_v, unit) in vars_agg
        impact = irf_v[1]
        ext_val = abs(minimum(irf_v)) > abs(maximum(irf_v)) ? minimum(irf_v) : maximum(irf_v)
        ext_q   = abs(minimum(irf_v)) > abs(maximum(irf_v)) ? argmin(irf_v) : argmax(irf_v)
        @printf(f, "%s & %+.3f & %+.3f & %d \\\\\n", label, impact, ext_val, ext_q)
    end

    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\begin{tablenotes}[flushleft]")
    println(f, "\\footnotesize")
    println(f, "\\item \\textit{Notes:} GDP, consumption, employment, and RER in \\% deviations from SS; inflation and policy rate in annualized pp; trade balance in pp of GDP. The shock is a one-time 10\\% increase in the world oil price \$P^{O*}_t\$ with persistence \$\\rho = $(rho_postar_val)\$. Impact = quarter 1 response. Trough/Peak = extreme response over 40 quarters.")
    println(f, "\\end{tablenotes}")
    println(f, "\\end{threeparttable}")
    println(f, "\\end{table}")
end
@printf "  Saved: oil_shock_aggregates.tex\n"
sync_table("oil_shock_aggregates.tex")


# ---- Sync section text to Overleaf ---- #
# oil_shock_section.tex lives next to the other tables and is \input{}-ed by Paper_soe.tex.
let sec_src = joinpath(TABLES_DIR, "oil_shock_section.tex")
    if isfile(sec_src)
        @printf "  Saved: oil_shock_section.tex (local)\n"
        sync_table("oil_shock_section.tex")
    end
end


# ---- Save IRF data as CSV for external use (baseline case) ---- #
# Reconstruct the IRF long-format DataFrame from the baseline case's stored get_irf
# Save IRF data as CSV (baseline case)
let
    vars_to_save = ["GDP", "pi", "Q", "TB", "C", "N", "w", "r"]
    for i in 1:nsec
        push!(vars_to_save, "Y_$(i)", "PH_$(i)", "MC_$(i)", "L_$(i)")
    end
    rows = []
    for vn in vars_to_save
        irf_v = get_irf(vn)
        for h in 1:nT
            push!(rows, (period=h, variable=vn, value=irf_v[h]))
        end
    end
    df_out = DataFrame(rows)
    CSV.write(joinpath(TABLES_DIR, "oil_shock_irfs_baseline.csv"), df_out)
    CSV.write(joinpath(TABLES_DIR, "oil_shock_irfs.csv"), df_out)
    @printf "  Saved: oil_shock_irfs_baseline.csv / oil_shock_irfs.csv\n"
end

# Save decomposition data
df_decomp = DataFrame(
    sector    = 1:nsec,
    name      = names_vec,
    oil_share = modalphaOil,
    alphaV    = modalphaV,
    direct_mc = direct_mc,
    network_mc= network_mc,
    total_mc  = total_mc_approx,
    amp_ratio = amp_ratio,
    mc_irf    = mc_irf_impact,
    y_irf     = y_irf_impact,
    ph_irf    = ph_irf_impact,
    yi_ss     = Yi_ss,
)
CSV.write(joinpath(TABLES_DIR, "oil_shock_decomposition.csv"), df_decomp)
@printf "  Saved: oil_shock_decomposition.csv\n"


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
