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
    decomposition_output_oil.pdf    — direct vs network effect on output
    decomposition_mc_oil.pdf        — direct vs network effect on MC
    oil_intensity_exposure.pdf      — sector oil exposure map
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
proc = run(pipeline(
    `$julia_exe --project=$project_dir $dynare_script $MOD_DIR`,
    stdout=dynare_out, stderr=stderr), wait=true)

dynare_stdout = String(take!(dynare_out))
print(dynare_stdout)

if !occursin("DYNARE_SUCCESS", dynare_stdout)
    error("Dynare subprocess failed. See output above.")
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
irf_pct = zeros(n_endo, n_irf)
for i in 1:n_endo
    ssv = abs(ss_vec[i])
    if ssv > 1e-10
        irf_pct[i, :] .= 100.0 .* irf_mat[i, :] ./ ssv
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
n_irf_v = get_irf_var("N")
r_irf   = get_irf_var("r")

@printf "--- Aggregate IRFs (10%% oil shock, impact period) ---\n"
@printf "  GDP    : %+.3f%%\n" gdp_irf[1]
@printf "  CPI (π): %+.3f%%\n" pi_irf[1]
@printf "  RER (Q): %+.3f%%\n" q_irf[1]
@printf "  TB     : %+.3f%%\n" tb_irf[1]
@printf "  C      : %+.3f%%\n" c_irf[1]
@printf "  N      : %+.3f%%\n" n_irf_v[1]
@printf "\n"

# Peak effects
gdp_trough = minimum(gdp_irf)
gdp_trough_q = argmin(gdp_irf)
pi_peak = maximum(pi_irf)
pi_peak_q = argmax(pi_irf)

@printf "  GDP trough : %+.3f%% at quarter %d\n" gdp_trough gdp_trough_q
@printf "  π peak     : %+.3f%% at quarter %d\n" pi_peak pi_peak_q


# =========================================================================== #
#  GENERATE FIGURES                                                            #
# =========================================================================== #

if _HAS_PLOTS[]
    @printf "\n--- Generating figures ---\n"
    periods = 1:n_irf

    # ---- Figure 1: Aggregate IRFs ---- #
    p_agg = Plots.plot(layout=(2,2), size=(1000,700),
        plot_title="10% Oil Price Shock — Aggregate Responses",
        titlefontsize=10, margin=5Plots.mm)

    Plots.plot!(p_agg, periods, gdp_irf, subplot=1,
        label="GDP", color=:steelblue, lw=2,
        ylabel="% dev.", title="GDP")
    Plots.hline!(p_agg, [0.0], subplot=1, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, pi_irf, subplot=2,
        label="Inflation (π)", color=:firebrick, lw=2,
        title="CPI Inflation")
    Plots.hline!(p_agg, [0.0], subplot=2, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, q_irf, subplot=3,
        label="RER (Q)", color=:forestgreen, lw=2,
        xlabel="Quarters", ylabel="% dev.", title="Real Exchange Rate")
    Plots.hline!(p_agg, [0.0], subplot=3, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, tb_irf, subplot=4,
        label="Trade Balance", color=:darkorange, lw=2,
        xlabel="Quarters", title="Trade Balance")
    Plots.hline!(p_agg, [0.0], subplot=4, color=:black, lw=0.6, ls=:dash, label="")

    Plots.savefig(p_agg, joinpath(FIGURES_DIR, "irf_aggregate_oil_shock.pdf"))
    @printf "  Saved: irf_aggregate_oil_shock.pdf\n"


    # ---- Figure 2: Sectoral Output IRFs (4×3 panel) ---- #
    short_names = [s[1:min(18,length(s))] for s in names_vec]
    sector_colors = vcat(
        fill(:steelblue, 5),   # goods sectors
        fill(:firebrick, 7)    # service sectors
    )

    p_sec_y = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Output Responses",
        titlefontsize=8)
    for i in 1:12
        yi = get_irf_var("Y_$(i)")
        Plots.plot!(p_sec_y, periods, yi, subplot=i,
            label="", color=sector_colors[i], lw=1.8,
            title=short_names[i], titlefontsize=7)
        Plots.hline!(p_sec_y, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    Plots.savefig(p_sec_y, joinpath(FIGURES_DIR, "irf_sectoral_Y_oil_shock.pdf"))
    @printf "  Saved: irf_sectoral_Y_oil_shock.pdf\n"


    # ---- Figure 3: Sectoral Price IRFs ---- #
    p_sec_ph = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Price Responses",
        titlefontsize=8)
    for i in 1:12
        phi = get_irf_var("PH_$(i)")
        Plots.plot!(p_sec_ph, periods, phi, subplot=i,
            label="", color=sector_colors[i], lw=1.8,
            title=short_names[i], titlefontsize=7)
        Plots.hline!(p_sec_ph, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    Plots.savefig(p_sec_ph, joinpath(FIGURES_DIR, "irf_sectoral_PH_oil_shock.pdf"))
    @printf "  Saved: irf_sectoral_PH_oil_shock.pdf\n"


    # ---- Figure 4: Decomposition Bar Chart — MC Direct vs Network ---- #
    p_decomp = groupedbar(
        [direct_mc network_mc],
        bar_position=:stack,
        xticks=(1:12, [s[1:min(12,length(s))] for s in names_vec]),
        xrotation=45,
        label=["Direct oil cost" "Network amplification"],
        color=[:darkorange :steelblue],
        ylabel="MC increase (%)",
        title="Oil Shock Decomposition: Direct vs. Network Effect on Marginal Cost",
        titlefontsize=10,
        size=(1000, 500),
        legend=:topright,
        bottom_margin=10Plots.mm
    )
    Plots.savefig(p_decomp, joinpath(FIGURES_DIR, "decomposition_mc_oil.pdf"))
    @printf "  Saved: decomposition_mc_oil.pdf\n"


    # ---- Figure 5: Oil Intensity Map ---- #
    p_oil = Plots.bar(1:12, modalphaOil .* 100,
        xticks=(1:12, [s[1:min(12,length(s))] for s in names_vec]),
        xrotation=45,
        label="Oil share of imports (%)",
        color=:darkorange,
        ylabel="Oil share (%)",
        title="Oil Intensity of Intermediate Imports by Sector",
        titlefontsize=10,
        size=(900, 450),
        bottom_margin=10Plots.mm
    )
    Plots.savefig(p_oil, joinpath(FIGURES_DIR, "oil_intensity_exposure.pdf"))
    @printf "  Saved: oil_intensity_exposure.pdf\n"


    # ---- Figure 6: Decomposition — Model IRFs for Output (impact) ---- #
    # Compare direct analytical effect with full model response
    p_out_decomp = groupedbar(
        [y_irf_impact mc_irf_impact],
        xticks=(1:12, [s[1:min(12,length(s))] for s in names_vec]),
        xrotation=45,
        label=["Output (Y_i, impact)" "Marginal Cost (MC_i, impact)"],
        color=[:steelblue :firebrick],
        ylabel="% deviation from SS",
        title="Oil Shock Impact: Sectoral Output and Marginal Cost",
        titlefontsize=10,
        size=(1000, 500),
        legend=:bottomright,
        bottom_margin=10Plots.mm
    )
    Plots.savefig(p_out_decomp, joinpath(FIGURES_DIR, "decomposition_output_oil.pdf"))
    @printf "  Saved: decomposition_output_oil.pdf\n"


    # ---- Figure 7: Sectoral MC IRFs (4×3 panel) ---- #
    p_sec_mc = Plots.plot(layout=(4,3), size=(1200,900),
        plot_title="10% Oil Shock — Sectoral Marginal Cost Responses",
        titlefontsize=8)
    for i in 1:12
        mci = get_irf_var("MC_$(i)")
        Plots.plot!(p_sec_mc, periods, mci, subplot=i,
            label="", color=sector_colors[i], lw=1.8,
            title=short_names[i], titlefontsize=7)
        Plots.hline!(p_sec_mc, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    Plots.savefig(p_sec_mc, joinpath(FIGURES_DIR, "irf_sectoral_MC_oil_shock.pdf"))
    @printf "  Saved: irf_sectoral_MC_oil_shock.pdf\n"

else
    @printf "\n  [Plots.jl not available — skipping figures]\n"
end


# =========================================================================== #
#  GENERATE LaTeX TABLES                                                       #
# =========================================================================== #

@printf "\n--- Generating LaTeX tables ---\n"

# ---- Table 1: Sectoral Decomposition ---- #
open(joinpath(TABLES_DIR, "oil_shock_decomposition.tex"), "w") do f
    println(f, "\\begin{table}[htbp]")
    println(f, "\\centering")
    println(f, "\\caption{Sectoral Decomposition of a 10\\% Oil Price Shock}")
    println(f, "\\label{tab:oil_decomp}")
    println(f, "\\begin{threeparttable}")
    println(f, "\\begin{tabular}{@{}l c c c c c c@{}}")
    println(f, "\\toprule")
    println(f, " & Oil & Direct & Network & Total & Amplif. & Output \\\\")
    println(f, "Sector & Share (\\%) & \\multicolumn{3}{c}{MC Increase (\\%)} & Ratio & IRF (\\%) \\\\")
    println(f, "\\midrule")

    for i in 1:nsec
        gs = goods[i] ? "G" : "S"
        @printf(f, "%s (%s) & %.1f & %.3f & %.3f & %.3f & %.2f & %+.3f \\\\\n",
            names_vec[i], gs,
            modalphaOil[i]*100, direct_mc[i], network_mc[i],
            total_mc_approx[i], amp_ratio[i], y_irf_impact[i])
    end

    println(f, "\\midrule")
    @printf(f, "Aggregate & --- & %.3f & %.3f & %.3f & %.2f & %+.3f \\\\\n",
        agg_direct, agg_network, agg_total, agg_amp, gdp_irf[1])
    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\begin{tablenotes}[flushleft]")
    println(f, "\\footnotesize")
    println(f, "\\item \\textit{Notes:} Oil Share is the fraction of sector \$i\$'s intermediate imports that are oil/refined fuels (from Chile's 2021 IO tables). Direct MC increase is \$\\alpha_{Vi} \\times \\alpha_{\\text{Oil},i} \\times \\Delta \\log P^O\$. Network effect is the additional cost transmitted through IO linkages via the Leontief inverse \$(I - \\text{diag}(\\alpha_m) \\cdot \\Gamma)^{-1}\$. Amplification ratio is Total/Direct. Output IRF is the impact-period response of sectoral gross output from the linearized model. G = Goods, S = Services.")
    println(f, "\\end{tablenotes}")
    println(f, "\\end{threeparttable}")
    println(f, "\\end{table}")
end
@printf "  Saved: oil_shock_decomposition.tex\n"


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
        ("GDP",          "GDP",   gdp_irf),
        ("CPI Inflation","pi",    pi_irf),
        ("Consumption",  "C",     c_irf),
        ("Employment",   "N",     n_irf_v),
        ("Real Exch. Rate","Q",   q_irf),
        ("Trade Balance","TB",    tb_irf),
        ("Policy Rate",  "r",     r_irf),
    ]
    for (label, _, irf_v) in vars_agg
        impact = irf_v[1]
        ext_val = abs(minimum(irf_v)) > abs(maximum(irf_v)) ? minimum(irf_v) : maximum(irf_v)
        ext_q   = abs(minimum(irf_v)) > abs(maximum(irf_v)) ? argmin(irf_v) : argmax(irf_v)
        @printf(f, "%s & %+.3f & %+.3f & %d \\\\\n", label, impact, ext_val, ext_q)
    end

    println(f, "\\bottomrule")
    println(f, "\\end{tabular}")
    println(f, "\\begin{tablenotes}[flushleft]")
    println(f, "\\footnotesize")
    println(f, "\\item \\textit{Notes:} All values in percentage deviations from steady state. The shock is a one-time 10\\% increase in the world oil price \$P^{O*}_t\$ with persistence \$\\rho = $(rho_postar_val)\$. Impact = quarter 1 response. Trough/Peak = extreme response over 40 quarters.")
    println(f, "\\end{tablenotes}")
    println(f, "\\end{threeparttable}")
    println(f, "\\end{table}")
end
@printf "  Saved: oil_shock_aggregates.tex\n"


# ---- Save IRF data as CSV for external use ---- #
df_irf_out = DataFrame(period = repeat(1:n_irf, outer=n_endo))
df_irf_out.variable = repeat(endo_names, inner=n_irf)
df_irf_out.value = vec(irf_pct')
CSV.write(joinpath(TABLES_DIR, "oil_shock_irfs.csv"), df_irf_out)
@printf "  Saved: oil_shock_irfs.csv\n"

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
@printf "  Figures → %s\n" FIGURES_DIR
@printf "  Tables  → %s\n" TABLES_DIR
@printf "%s\n\n" repeat("=", 70)

end  # function _main()

Base.invokelatest(_main)
