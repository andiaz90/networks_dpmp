"""
agr_shock_analysis.jl
=====================
NK-IOSOE 12-sector model for Chile — Agriculture TFP Shock Analysis (El Niño)

Computes the response of the Chilean economy to a negative 10% TFP shock in
the Agriculture sector (sector 1), motivated by El Niño events that
disrupt crop yields and agricultural production in Chile.

Generates publication-quality figures (PDF) saved directly to the
Overleaf Dropbox folder.

Workflow:
  1.  Load sector-level data (same as main_SOE_gap.jl / oil_shock_analysis.jl)
  2.  Set structural parameters (SMM estimates from smm_estimates.csv)
  3.  Configure Exercise 2: Agriculture TFP shock only (epsA_1)
  4.  Solve steady state and write params_jl.mod
  5.  Run Dynare (subprocess) to get decision rules
  6.  Compute IRFs to a negative 1% TFP shock in sector 1
  7.  Compute expenditure-weighted sectoral inflation decomposition
  8.  Generate PDF figures to figures/agr_shock/ AND Overleaf folder
  9.  (Optional) Generate comparison figures vs. oil shock

Usage:
  julia --project=. agr_shock_analysis.jl

REQUIREMENTS:
  Same as oil_shock_analysis.jl (Julia >= 1.9, Dynare.jl, CSV, DataFrames, etc.)

OUTPUT FILES:
  figures/agr_shock/irf_aggregate_agr_shock.pdf       — GDP, CPI, Q, r, TB
  figures/agr_shock/irf_sectoral_Y_agr_shock.pdf      — 12-sector output IRFs
  figures/agr_shock/irf_sectoral_inflation_agr_shock.pdf — 12-sector ΔPH_i inflation
  figures/agr_shock/irf_aggregate_inflation_agr_shock.pdf — goods/services/agg inflation
  figures/agr_shock/irf_comparison_oil_agr.pdf        — comparison vs oil (if oil data exists)
  figures/agr_shock/irf_comparison_sectoral_inflation.pdf — sectoral comparison (if oil data exists)

  All figures also copied to:
  ../../Network DPMP-DME/Figures/agr_shock/   (Overleaf Dropbox sync folder)
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
#  INCLUDE HELPERS                                                             #
# =========================================================================== #

SCRIPT_DIR = @__DIR__

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))

# Plotting (graceful failure if not installed)
const _HAS_PLOTS = Ref(false)
try
    @eval Main using Plots
    @eval Main gr(dpi=200)
    _HAS_PLOTS[] = true
catch
    @warn "Plots.jl not available — figures will be skipped. Install with: ] add Plots"
end

function _main()

@printf "\n%s\n" repeat("=", 70)
@printf "  NK-SOE 12-sector model — Agriculture TFP Shock Analysis (El Niño)\n"
@printf "%s\n\n" repeat("=", 70)


# =========================================================================== #
#  PATHS                                                                       #
# =========================================================================== #

MOD_DIR    = joinpath(SCRIPT_DIR, "mod")
REPO_ROOT  = abspath(joinpath(SCRIPT_DIR, "..", ".."))
DATA_DIR   = joinpath(REPO_ROOT, "Data")

FIGURES_DIR = joinpath(SCRIPT_DIR, "figures", "agr_shock")
TABLES_DIR  = joinpath(SCRIPT_DIR, "tables")
mkpath(FIGURES_DIR)
mkpath(TABLES_DIR)

# Overleaf Dropbox sync folder — figures are written here for direct inclusion
OVERLEAF_DIR = get(ENV, "OVERLEAF_AGR",
    abspath(joinpath(homedir(), "Library", "CloudStorage",
        "Dropbox", "Apps", "Overleaf", "Network DPMP-DME", "Figures", "agr_shock")))
overleaf_ok = try mkpath(OVERLEAF_DIR); true catch; false end
if overleaf_ok
    @printf "  Overleaf output: %s\n\n" OVERLEAF_DIR
else
    @printf "  [Overleaf folder not found — saving locally only]\n\n"
end

function save_fig(p, fname)
    local_path = joinpath(FIGURES_DIR, fname)
    @eval Main savefig($p, $local_path)
    @printf "  Saved: %s\n" local_path
    if overleaf_ok
        ol_path = joinpath(OVERLEAF_DIR, fname)
        @eval Main savefig($p, $ol_path)
        @printf "  Synced: %s\n" ol_path
    end
end


# =========================================================================== #
#  READ DATA (identical to oil_shock_analysis.jl)                              #
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

# ---- Oil parameters (needed for params_jl.mod even in Exercise 2) ----
modalphaOil = [0.1635, 0.2164, 0.1890, 0.0871, 0.0417,
               0.0793, 0.3734, 0.0043, 0.0449, 0.0633, 0.0391, 0.0447]
epsilonV_oil_val = 0.5
rho_postar_val   = 0.9
sigma_postar_val = 0.02
POstar_ss_val    = 1.0

# ---- Exercise 2: AGRICULTURE TFP SHOCK ONLY (El Niño) ----
agr_sector = 1   # sector 1 = Agriculture, Forestry & Fishing

sigma_i_val     = 0.0
rho_om1_val     = 0.1
sigma_om_vec    = zeros(nsec)
rho_tfp1_val    = 0.95;  rho_tfp2_val = 0.0   # Exercise 2 default persistence
isigma_tfp_val  = zeros(nsec)
isigma_tfp_val[agr_sector] = 0.01              # only Agriculture TFP active
rho_val         = 0.1
sigma_L_agg_val = 0.0
rho_pvstar_val  = 0.9; sigma_pvstar_val = 0.0
rho_xi_val      = 0.80; sigma_xi_val = 0.0

# Shock activation flags
shock_eps_om_vec     = zeros(nsec)
shock_eps_i_val      = 0.0
shock_eps_pvstar_val = 0.0
shock_eps_xi_val     = 0.0
shock_eps_postar_val = 0.0         # oil OFF
shock_epsA_val       = zeros(nsec)
shock_epsA_val[agr_sector] = 1.0  # << Agriculture TFP ON

# ---- Load SMM estimates (override defaults) ----
smm_est_file = joinpath(DATA_DIR, "smm_estimates.csv")
if isfile(smm_est_file)
    est_df = CSV.read(smm_est_file, DataFrame)
    est    = Dict(String(r.param) => Float64(r.value) for r in eachrow(est_df))
    ilabcosts_val   = est["ilabcosts"]
    modepsY         = fill(est["epsY"], nsec)
    modepsM         = fill(est["epsM"], nsec)
    kappaV_val      = exp(est["log_kappaV"])
    rho_om1_val     = est["rho_om"]
    rho_tfp1_val    = est["rho_A"]          # use SMM-estimated AR(1) persistence
    isigma_tfp_val  = [est["isigma_tfp_$(i)"] for i in 1:nsec]
    # Override: all TFP sigmas zero except Agriculture
    sigma_A1        = isigma_tfp_val[agr_sector]   # keep for scale factor
    fill!(isigma_tfp_val, 0.0)
    isigma_tfp_val[agr_sector] = sigma_A1
    haskey(est, "etastar") && (etastar_val = est["etastar"])
    @printf "  Loaded SMM estimates from %s\n" smm_est_file
    @printf "  σ_A1 (Agriculture TFP std) = %.5f\n" sigma_A1
else
    sigma_A1 = isigma_tfp_val[agr_sector]
    @printf "  SMM estimates not found — using defaults. σ_A1 = %.4f\n" sigma_A1
end


# =========================================================================== #
#  STEADY STATE (identical to oil_shock_analysis.jl)                           #
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

PO_ss      = Q_ss * POstar_ss_val
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

# Validate decision rules
_nz_ghx = sum(abs.(ghx) .> 1e-12)
_nz_ghu = sum(abs.(ghu) .> 1e-12)
@printf "  Decision rule nonzeros:  ghx=%d  ghu=%d\n" _nz_ghx _nz_ghu
if _nz_ghx == 0 || _nz_ghu == 0
    error("Decision rules are all zeros — Dynare LRE solver failed")
end

# varexo declaration order (from NK_SOE_lev_gap2.mod):
#   eps_om(1) eps_i(2) epschi(3) eps_pvstar(4) eps_postar(5)
#   epsA_1(6) epsA_2(7) epsA_3(8) ... epsA_12(17) eps_xi(18)
agr_col = 5 + agr_sector   # = 6 for sector 1
@printf "  Agriculture TFP shock (epsA_%d) column: %d\n\n" agr_sector agr_col


# =========================================================================== #
#  COMPUTE IRFs: NEGATIVE 1% AGRICULTURE TFP SHOCK (El Niño)                   #
#                                                                              #
#  In the .mod: A_1 = rho_tfp1*A_1(-1) + isigma_tfp_1 * epsA_1              #
#  A negative 1% TFP shock means d(log A_1) = -0.01 on impact.              #
#  Scale: eps_A1 = -0.01 / sigma_A1   (negative sign for El Niño)            #
# =========================================================================== #

shock_pct    = -0.10   # NEGATIVE 10% TFP shock (El Niño reduces productivity)
scale_factor = shock_pct / sigma_A1

@printf "--- Computing IRFs (negative 10%% Agriculture TFP shock — El Niño) ---\n"
@printf "  σ_A1 = %.5f → scale factor = %.4f s.d.\n\n" sigma_A1 scale_factor

n_irf = 40   # quarters

# State-space matrices
A_mat = ghx[state_rows, :]  # n_states × n_states
B_mat = ghu[state_rows, :]  # n_states × n_shocks

# Compute scaled IRFs
irf_mat = zeros(n_endo, n_irf)
x = zeros(n_states)

for h in 1:n_irf
    if h == 1
        y_h = ghu[:, agr_col] .* scale_factor
        x   = B_mat[:, agr_col] .* scale_factor
    else
        y_h = ghx * x
        x   = A_mat * x
    end
    irf_mat[:, h] .= y_h
end

# Convert to percentage deviations from SS
irf_pct = zeros(n_endo, n_irf)
for i in 1:n_endo
    ssv = abs(ss_vec[i])
    if ssv > 1e-10
        irf_pct[i, :] .= 100.0 .* irf_mat[i, :] ./ ssv
    else
        irf_pct[i, :] .= 100.0 .* irf_mat[i, :]
    end
end

function get_irf_var(varname::String)
    idx = get(endo_idx, varname, 0)
    idx == 0 && return zeros(n_irf)
    return irf_pct[idx, :]
end

# Extract aggregate CPI inflation
pi_cpi_irf = get_irf_var("pi")

# Sectoral nominal home-price inflation
function get_inflation_irf(i::Int)
    ph = get_irf_var("PH_$(i)")
    dph = similar(ph)
    dph[1] = ph[1] + pi_cpi_irf[1]
    for t in 2:n_irf
        dph[t] = ph[t] - ph[t-1] + pi_cpi_irf[t]
    end
    return dph .* 4   # annualize
end

# Build infl_irf as a proper (nsec × n_irf) Matrix
infl_irf = zeros(nsec, n_irf)
for i in 1:nsec
    infl_irf[i, :] .= get_inflation_irf(i)
end

# Expenditure-weighted aggregate inflation
omega_G = sum(spend_good)   # ≈ 0.57
omega_S = sum(spend_serv)   # ≈ 0.43

agg_infl_g   = (spend_good' * infl_irf)[:]  ./ omega_G
agg_infl_s   = (spend_serv' * infl_irf)[:]  ./ omega_S
agg_infl_tot = omega_G .* agg_infl_g .+ omega_S .* agg_infl_s


# =========================================================================== #
#  PRINT KEY RESULTS                                                           #
# =========================================================================== #

gdp_irf = get_irf_var("GDP")
pi_irf  = pi_cpi_irf
q_irf   = get_irf_var("Q")
tb_irf  = get_irf_var("TB")
c_irf   = get_irf_var("C")
cg_irf  = get_irf_var("C_g")
cs_irf  = get_irf_var("C_s")
n_irf_v = get_irf_var("N")
r_irf   = get_irf_var("r")
w_irf   = get_irf_var("w")

@printf "--- Aggregate IRFs (negative 10%% Agriculture TFP shock — El Niño, impact period) ---\n"
@printf "  GDP      : %+.3f%%\n"   gdp_irf[1]
@printf "  CPI (π)  : %+.3f ann pp\n" pi_irf[1] * 4
@printf "  RER (Q)  : %+.3f%%\n"   q_irf[1]
@printf "  r        : %+.3f ann pp\n" r_irf[1] * 4
@printf "  TB/GDP   : %+.3f pp\n"  100 * (tb_irf[1] / 100 * TB_ss / GDP_ss)
@printf "  PH_1     : %+.3f%%\n"   get_irf_var("PH_$(agr_sector)")[1]
@printf "\n"
@printf "  Home-price inflation (annualized pp, impact):\n"
@printf "    Goods    : %+.3f\n"  agg_infl_g[1]
@printf "    Services : %+.3f\n"  agg_infl_s[1]
@printf "    Aggregate: %+.3f\n"  agg_infl_tot[1]
@printf "\n"

# Sectoral inflation impact
@printf "  %-45s  %8s\n" "Sector" "ΔPH (pp)"
@printf "  %s\n" repeat("-", 58)
for i in 1:nsec
    gs = goods[i] ? "G" : "S"
    @printf "  %-43s (%s)  %+7.3f\n" names_vec[i] gs infl_irf[i, 1]
end
@printf "\n"


# =========================================================================== #
#  DECOMPOSITION: DIRECT TFP COST vs. NETWORK AMPLIFICATION                   #
#                                                                              #
#  For the agriculture TFP shock, the "direct" effect hits only sector 1's    #
#  marginal cost (through lower TFP → higher unit cost).  The "network"       #
#  effect captures how that cost increase propagates to downstream sectors    #
#  via the IO matrix Γ and the Leontief inverse.                              #
#                                                                              #
#  GE decomposition of MC at each horizon h:                                  #
#    direct_h(i)  = -ΔA^1_h  if i = agriculture, 0 otherwise                 #
#    network_h(i) = α_mi × (Γ × PH_h)_i   (IO cost pass-through)            #
#    oil_h(i)     = αV_i × αO_i × PO_h     (oil import cost)                 #
#    nonoil_h(i)  = αV_i × (1-αO_i) × PV_h (non-oil import cost)            #
#    labor_h(i)   = αL_i × w_h              (wage channel)                    #
#    residual     = MC_h - sum of above                                       #
# =========================================================================== #

@printf "--- Decomposing direct vs. network effects (Agriculture shock) ---\n"

# Build MC and PH IRF matrices (nsec × n_irf)
mc_irf_mat = zeros(nsec, n_irf)
ph_irf_mat = zeros(nsec, n_irf)
for i in 1:nsec
    mc_irf_mat[i, :] .= get_irf_var("MC_$(i)")
    ph_irf_mat[i, :] .= get_irf_var("PH_$(i)")
end

po_irf = get_irf_var("PO")
pv_irf = get_irf_var("PV")

alpha_L_vec = 1.0 .- modalpha .- modalphaV

# GE decomposition at a given horizon
function _mc_decomp_agr(h, _ph, _mc, _po, _pv, _w, _aV, _aO, _am, _aL, _beta)
    # Direct oil channel (small for agr shock but nonzero through GE)
    dir_oil = _aV .* _aO .* _po[h]
    # Network channel: IO cost propagation
    net_h   = _am .* (_beta * _ph[:, h])
    # Non-oil import channel
    noi_h   = _aV .* (1.0 .- _aO) .* _pv[h]
    # Labor channel
    lab_h   = _aL .* _w[h]
    # Residual (includes direct TFP effect on sector 1 + higher-order terms)
    res_h   = _mc[:, h] .- dir_oil .- net_h .- noi_h .- lab_h
    return (dir=copy(dir_oil), net=copy(net_h), noi=copy(noi_h),
            lab=copy(lab_h), res=copy(res_h), mc=copy(_mc[:, h]))
end

ge_agr_h1 = _mc_decomp_agr(1, ph_irf_mat, mc_irf_mat, po_irf, pv_irf, w_irf,
                            modalphaV, modalphaOil, modalpha, alpha_L_vec, modbeta)
ge_agr_h2 = _mc_decomp_agr(2, ph_irf_mat, mc_irf_mat, po_irf, pv_irf, w_irf,
                            modalphaV, modalphaOil, modalpha, alpha_L_vec, modbeta)
ge_agr_h4 = _mc_decomp_agr(4, ph_irf_mat, mc_irf_mat, po_irf, pv_irf, w_irf,
                            modalphaV, modalphaOil, modalpha, alpha_L_vec, modbeta)

# Print decomposition at impact
@printf "\n  %-30s  %8s  %8s  %8s  %8s  %8s  %8s\n" "Sector" "Oil" "Network" "NonOilImp" "Labor" "Resid" "MC(model)"
@printf "  %s\n" repeat("-", 100)
for i in 1:nsec
    @printf "  %-30s  %+7.3f  %+7.3f  %+7.3f  %+7.3f  %+7.3f  %+7.3f\n" names_vec[i] ge_agr_h1.dir[i] ge_agr_h1.net[i] ge_agr_h1.noi[i] ge_agr_h1.lab[i] ge_agr_h1.res[i] ge_agr_h1.mc[i]
end
@printf "\n"

# Inflation at various horizons (for diamond markers)
infl_irf_impact_agr = ph_irf_mat[:, 1] .+ pi_irf[1]
infl_6m_agr  = ph_irf_mat[:, 2] .+ sum(pi_irf[1:2])
infl_12m_agr = ph_irf_mat[:, 4] .+ sum(pi_irf[1:4])


# =========================================================================== #
#  GENERATE FIGURES                                                            #
# =========================================================================== #

if !_HAS_PLOTS[]
    @printf "\n  [Plots.jl not available — skipping figures]\n"
else
    @printf "\n--- Generating figures ---\n"
    periods = 1:n_irf
    short_names = [s[1:min(18, length(s))] for s in names_vec]
    sector_colors = vcat(fill(:steelblue, 5), fill(:firebrick, 7))

    # ------------------------------------------------------------------ #
    # Figure 1: Aggregate IRFs (5 panels)                                 #
    # ------------------------------------------------------------------ #
    p_agg = Plots.plot(layout=(2, 5), size=(2000, 700),
        plot_title="10% Agriculture TFP Decline (El Niño) — Aggregate Responses",
        titlefontsize=10, margin=5Plots.mm)

    # Row 1: GDP, CPI, RER, TB, Goods Consumption
    Plots.plot!(p_agg, periods, gdp_irf, subplot=1,
        label="GDP", color=:steelblue, lw=2,
        ylabel="% dev. from SS", title="GDP")
    Plots.hline!(p_agg, [0.0], subplot=1, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, pi_irf .* 4, subplot=2,
        label="CPI inflation", color=:firebrick, lw=2,
        title="CPI Inflation (ann. pp)")
    Plots.hline!(p_agg, [0.0], subplot=2, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, q_irf, subplot=3,
        label="RER (Q)", color=:forestgreen, lw=2,
        title="Real Exchange Rate")
    Plots.hline!(p_agg, [0.0], subplot=3, color=:black, lw=0.6, ls=:dash, label="")

    tb_gdp_irf = tb_irf .* (TB_ss / GDP_ss)
    Plots.plot!(p_agg, periods, tb_gdp_irf, subplot=4,
        label="TB/GDP", color=:darkorange, lw=2,
        title="Trade Balance (% GDP)")
    Plots.hline!(p_agg, [0.0], subplot=4, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, cg_irf, subplot=5,
        label="Goods C", color=:steelblue, lw=2,
        title="Goods Consumption")
    Plots.hline!(p_agg, [0.0], subplot=5, color=:black, lw=0.6, ls=:dash, label="")

    # Row 2: Policy Rate, Real Wage, Goods Inflation, Services Inflation, Services Consumption
    Plots.plot!(p_agg, periods, r_irf .* 4, subplot=6,
        label="Policy rate", color=:darkorchid, lw=2,
        xlabel="Quarters", ylabel="ann. pp dev.", title="Policy Rate")
    Plots.hline!(p_agg, [0.0], subplot=6, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, w_irf, subplot=7,
        label="Real wage", color=:teal, lw=2,
        xlabel="Quarters", title="Real Wage")
    Plots.hline!(p_agg, [0.0], subplot=7, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, agg_infl_g, subplot=8,
        label="Goods", color=:steelblue, lw=2,
        xlabel="Quarters", title="Goods Inflation (ann. pp)")
    Plots.hline!(p_agg, [0.0], subplot=8, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, agg_infl_s, subplot=9,
        label="Services", color=:firebrick, lw=2,
        xlabel="Quarters", title="Services Inflation (ann. pp)")
    Plots.hline!(p_agg, [0.0], subplot=9, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_agg, periods, cs_irf, subplot=10,
        label="Services C", color=:firebrick, lw=2,
        xlabel="Quarters", title="Services Consumption")
    Plots.hline!(p_agg, [0.0], subplot=10, color=:black, lw=0.6, ls=:dash, label="")

    save_fig(p_agg, "irf_aggregate_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 2: Sectoral Output IRFs (4×3 panel)                          #
    # ------------------------------------------------------------------ #
    p_sec_y = Plots.plot(layout=(4, 3), size=(1200, 900),
        plot_title="10% Agriculture TFP Decline (El Niño) — Sectoral Output",
        titlefontsize=8)
    for i in 1:12
        yi = get_irf_var("Y_$(i)")
        clr = (i == agr_sector) ? :firebrick : sector_colors[i]
        lw_i = (i == agr_sector) ? 2.5 : 1.8
        Plots.plot!(p_sec_y, periods, yi, subplot=i,
            label="", color=clr, lw=lw_i,
            title=short_names[i], titlefontsize=7,
            ylabel=(i % 3 == 1 ? "% dev." : ""))
        Plots.hline!(p_sec_y, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_y, "irf_sectoral_Y_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 3: Sectoral Home-Price Inflation IRFs (4×3 panel)            #
    # ------------------------------------------------------------------ #
    p_sec_pi = Plots.plot(layout=(4, 3), size=(1200, 900),
        plot_title="10% Agriculture TFP Decline (El Niño) — Sectoral Inflation (ann. pp)",
        titlefontsize=8)
    for i in 1:12
        pi_i = get_inflation_irf(i)
        Plots.plot!(p_sec_pi, periods, pi_i, subplot=i,
            label="", color=sector_colors[i], lw=1.8,
            title=short_names[i], titlefontsize=7,
            ylabel=(i % 3 == 1 ? "ann. pp" : ""))
        Plots.hline!(p_sec_pi, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_pi, "irf_sectoral_inflation_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 4: Aggregate Inflation Decomposition (Goods/Services/Total)  #
    # ------------------------------------------------------------------ #
    p_agg_pi = Plots.plot(size=(800, 500),
        title="10% Agriculture TFP Decline (El Niño) — Expenditure-Weighted Inflation",
        titlefontsize=10, xlabel="Quarters", ylabel="Annualized pp dev. from SS",
        legend=:topright, margin=5Plots.mm)
    Plots.plot!(p_agg_pi, periods, agg_infl_tot, label="Aggregate",
        color=:black, lw=2.5, ls=:solid)
    Plots.plot!(p_agg_pi, periods, agg_infl_g, label="Goods (sectors 1–5)",
        color=:steelblue, lw=2, ls=:dash)
    Plots.plot!(p_agg_pi, periods, agg_infl_s, label="Services (sectors 6–12)",
        color=:firebrick, lw=2, ls=:dot)
    Plots.hline!(p_agg_pi, [0.0], color=:black, lw=0.6, ls=:dash, label="")
    save_fig(p_agg_pi, "irf_aggregate_inflation_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 4b: Aggregate Labor Market IRFs (N and w)                    #
    # ------------------------------------------------------------------ #
    w_irf_full = get_irf_var("w")
    n_irf_lab  = get_irf_var("N")

    p_lab_agg = Plots.plot(layout=(1, 2), size=(1000, 400),
        plot_title="10% Agriculture TFP Decline (El Niño) — Aggregate Labor Market",
        titlefontsize=10, margin=5Plots.mm)

    Plots.plot!(p_lab_agg, periods, n_irf_lab, subplot=1,
        label="Employment (N)", color=:steelblue, lw=2,
        xlabel="Quarters", ylabel="% dev. from SS", title="Aggregate Employment")
    Plots.hline!(p_lab_agg, [0.0], subplot=1, color=:black, lw=0.6, ls=:dash, label="")

    Plots.plot!(p_lab_agg, periods, w_irf_full, subplot=2,
        label="Wage (w)", color=:firebrick, lw=2,
        xlabel="Quarters", ylabel="% dev. from SS", title="Real Wage")
    Plots.hline!(p_lab_agg, [0.0], subplot=2, color=:black, lw=0.6, ls=:dash, label="")

    save_fig(p_lab_agg, "irf_labor_aggregate_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 4c: Sectoral Employment IRFs (4×3 panel)                     #
    # ------------------------------------------------------------------ #
    p_sec_l = Plots.plot(layout=(4, 3), size=(1200, 900),
        plot_title="10% Agriculture TFP Decline (El Niño) — Sectoral Employment",
        titlefontsize=8)
    for i in 1:12
        li = get_irf_var("L_$(i)")
        clr = (i == agr_sector) ? :firebrick : sector_colors[i]
        lw_i = (i == agr_sector) ? 2.5 : 1.8
        Plots.plot!(p_sec_l, periods, li, subplot=i,
            label="", color=clr, lw=lw_i,
            title=short_names[i], titlefontsize=7,
            ylabel=(i % 3 == 1 ? "% dev." : ""))
        Plots.hline!(p_sec_l, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_l, "irf_sectoral_L_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 4d: GDP Gap IRF                                              #
    # ------------------------------------------------------------------ #
    gdpgap_irf = get_irf_var("GDPgap")

    p_gap_agg = Plots.plot(size=(700, 400),
        title="10% Agriculture TFP Decline (El Niño) — GDP Gap",
        titlefontsize=10, margin=5Plots.mm,
        xlabel="Quarters", ylabel="% dev. from SS")

    Plots.plot!(p_gap_agg, periods, gdpgap_irf,
        label="GDP gap", color=:firebrick, lw=2)
    Plots.hline!(p_gap_agg, [0.0], color=:black, lw=0.6, ls=:dash, label="")

    save_fig(p_gap_agg, "irf_gdpgap_aggregate_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 4e: Sectoral Output Gap IRFs (4×3 panel)                     #
    # ------------------------------------------------------------------ #
    p_sec_ygap = Plots.plot(layout=(4, 3), size=(1200, 900),
        plot_title="10% Agriculture TFP Decline (El Niño) — Sectoral Output Gap",
        titlefontsize=8)
    for i in 1:12
        ygi = get_irf_var("Ygap_$(i)")
        clr = (i == agr_sector) ? :firebrick : sector_colors[i]
        lw_i = (i == agr_sector) ? 2.5 : 1.8
        Plots.plot!(p_sec_ygap, periods, ygi, subplot=i,
            label="", color=clr, lw=lw_i,
            title=short_names[i], titlefontsize=7,
            ylabel=(i % 3 == 1 ? "% dev." : ""))
        Plots.hline!(p_sec_ygap, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
    end
    save_fig(p_sec_ygap, "irf_sectoral_Ygap_agr_shock.pdf")


    # ------------------------------------------------------------------ #
    # Figure 5: MC Decomposition — direct vs. network (stacked bars)      #
    # ------------------------------------------------------------------ #
    bar_names_agr = [
        "Agric. & Fishing", "Mining", "Manufacturing", "Utilities",
        "Construction", "Trade & Hotels", "Transport & ICT", "Finance",
        "Real Estate", "Business Serv.", "Personal Serv.", "Public Admin.",
    ]

    bar_rect_agr(x, y0, y1, w=0.65) = Plots.Shape(
        [x - w/2, x + w/2, x + w/2, x - w/2],
        [y0, y0, y1, y1])

    # For the agriculture shock, relabel: "Direct oil" becomes "Direct TFP + oil"
    # since the residual captures the TFP effect on sector 1
    ge_comp_colors_agr = [:darkorange, :steelblue, :forestgreen, :crimson, :gray60]
    ge_comp_labels_agr = ["Direct oil (pp)", "Network via prices (pp)",
                          "Non-oil imports (pp)", "Labor (pp)", "TFP + residual (pp)"]

    function make_ge_decomp_fig_agr(dir_h, net_h, noimp_h, lab_h, resid_h, mc_h,
                                     infl_h, infl_label, title_str, ylabel_str)
        comp_mat = hcat(dir_h, net_h, noimp_h, lab_h, resid_h)
        pos_mat = max.(comp_mat, 0.0)
        neg_mat = min.(comp_mat, 0.0)
        pos_tops = vec(sum(pos_mat, dims=2))
        neg_tops = vec(sum(neg_mat, dims=2))
        all_vals = vcat(pos_tops, neg_tops, infl_h, [0.0])
        ylo = min(minimum(all_vals), 0.0) * 1.40
        yhi = max(maximum(all_vals), 0.0) * 1.35
        p = Plots.plot(
            xticks=(1:nsec, bar_names_agr), xrotation=55,
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
                lbl = lbl_used ? "" : ge_comp_labels_agr[k]
                if pos_mat[i, k] > 1e-10
                    s = bar_rect_agr(i, pos_bot[i], pos_bot[i] + pos_mat[i, k], bw)
                    Plots.plot!(p, s, color=ge_comp_colors_agr[k], label=lbl,
                        alpha=0.85, linecolor=:white, linewidth=0.3)
                    lbl_used = true; lbl = ""
                end
                if abs(neg_mat[i, k]) > 1e-10
                    s = bar_rect_agr(i, neg_bot[i] + neg_mat[i, k], neg_bot[i], bw)
                    Plots.plot!(p, s, color=ge_comp_colors_agr[k], label=lbl,
                        alpha=0.85, linecolor=:white, linewidth=0.3)
                    lbl_used = true
                end
            end
            if !lbl_used
                Plots.plot!(p, [NaN], [NaN], color=ge_comp_colors_agr[k],
                    label=ge_comp_labels_agr[k], lw=4)
            end
        end
        Plots.scatter!(p, 1:nsec, infl_h,
            color=:black, markershape=:diamond,
            markersize=7, markerstrokewidth=1, label=infl_label)
        return p
    end

    # Impact (h = 1)
    p_decomp1 = make_ge_decomp_fig_agr(
        ge_agr_h1.dir, ge_agr_h1.net, ge_agr_h1.noi, ge_agr_h1.lab, ge_agr_h1.res, ge_agr_h1.mc,
        infl_irf_impact_agr,
        "Home-price inflation (q-o-q pp, t=1)",
        "Agriculture TFP Shock — GE MC Decomposition (Impact)",
        "Percentage points (t = 1)")
    save_fig(p_decomp1, "decomposition_mc_inflation_agr_baseline.pdf")

    # 6-month (h = 2)
    p_decomp2 = make_ge_decomp_fig_agr(
        ge_agr_h2.dir, ge_agr_h2.net, ge_agr_h2.noi, ge_agr_h2.lab, ge_agr_h2.res, ge_agr_h2.mc,
        infl_6m_agr,
        "Home-price inflation (6-month cumulative, pp)",
        "Agriculture TFP Shock — GE MC Decomposition (6-Month)",
        "Percentage points")
    save_fig(p_decomp2, "decomposition_mc_inflation_6m_agr_baseline.pdf")

    # 12-month (h = 4)
    p_decomp4 = make_ge_decomp_fig_agr(
        ge_agr_h4.dir, ge_agr_h4.net, ge_agr_h4.noi, ge_agr_h4.lab, ge_agr_h4.res, ge_agr_h4.mc,
        infl_12m_agr,
        "Home-price inflation (12-month cumulative, pp)",
        "Agriculture TFP Shock — GE MC Decomposition (12-Month)",
        "Percentage points")
    save_fig(p_decomp4, "decomposition_mc_inflation_12m_agr_baseline.pdf")

    @printf "  Decomposition figures saved (agriculture shock)\n"


    # ------------------------------------------------------------------ #
    # Figure 6 & 7: Comparison with Oil Shock (if oil IRF data exists)    #
    # ------------------------------------------------------------------ #
    oil_irfs_path = joinpath(TABLES_DIR, "oil_shock_irfs.csv")
    oil_decomp_path = joinpath(TABLES_DIR, "oil_shock_decomposition.csv")

    if isfile(oil_irfs_path)
        @printf "  Loading oil shock IRFs from %s\n" oil_irfs_path
        df_oil = CSV.read(oil_irfs_path, DataFrame)

        function get_oil_irf(varname::String)
            sub = filter(r -> String(r.variable) == varname, df_oil)
            isempty(sub) && return zeros(n_irf)
            s = sort(sub, :period)
            out = zeros(n_irf)
            n = min(n_irf, nrow(s))
            out[1:n] .= s.value[1:n]
            return out
        end

        # Oil sectoral nominal inflation
        pi_oil_cpi = get_oil_irf("pi")
        function get_oil_inflation_irf(i::Int)
            ph = get_oil_irf("PH_$(i)")
            dph = similar(ph)
            dph[1] = ph[1] + pi_oil_cpi[1]
            for t in 2:n_irf
                dph[t] = ph[t] - ph[t-1] + pi_oil_cpi[t]
            end
            return dph .* 4
        end

        # --- Figure 5: Aggregate comparison (4 panels) ---
        p_comp = Plots.plot(layout=(2, 2), size=(1000, 700),
            plot_title="Oil Price Shock vs. Agriculture TFP Shock (El Niño)",
            titlefontsize=10, margin=5Plots.mm)

        Plots.plot!(p_comp, periods, get_oil_irf("GDP"), subplot=1,
            label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
        Plots.plot!(p_comp, periods, gdp_irf, subplot=1,
            label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
            ylabel="% dev.", title="GDP")
        Plots.hline!(p_comp, [0.0], subplot=1, color=:black, lw=0.5, ls=:dot, label="")

        Plots.plot!(p_comp, periods, get_oil_irf("pi") .* 4, subplot=2,
            label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
        Plots.plot!(p_comp, periods, pi_irf .* 4, subplot=2,
            label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
            title="CPI Inflation (ann. pp)")
        Plots.hline!(p_comp, [0.0], subplot=2, color=:black, lw=0.5, ls=:dot, label="")

        Plots.plot!(p_comp, periods, get_oil_irf("Q"), subplot=3,
            label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
        Plots.plot!(p_comp, periods, q_irf, subplot=3,
            label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
            xlabel="Quarters", ylabel="% dev.", title="Real Exchange Rate (Q)")
        Plots.hline!(p_comp, [0.0], subplot=3, color=:black, lw=0.5, ls=:dot, label="")

        Plots.plot!(p_comp, periods, get_oil_irf("r") .* 4, subplot=4,
            label="Oil shock (10%)", color=:darkorange, lw=2, ls=:solid)
        Plots.plot!(p_comp, periods, r_irf .* 4, subplot=4,
            label="Agr TFP (−10%)", color=:steelblue, lw=2, ls=:dash,
            xlabel="Quarters", title="Policy Rate (ann. pp)")
        Plots.hline!(p_comp, [0.0], subplot=4, color=:black, lw=0.5, ls=:dot, label="")

        save_fig(p_comp, "irf_comparison_oil_agr.pdf")


        # --- Figure 6: Sectoral inflation comparison (4×3 panel) ---
        p_comp_sec = Plots.plot(layout=(4, 3), size=(1200, 900),
            plot_title="Sectoral Inflation: Oil vs. Agriculture TFP Shock (El Niño)",
            titlefontsize=8)
        for i in 1:12
            pi_oil_i = get_oil_inflation_irf(i)
            pi_agr_i = get_inflation_irf(i)
            Plots.plot!(p_comp_sec, periods, pi_oil_i, subplot=i,
                label=(i == 1 ? "Oil (10%)" : ""), color=:darkorange, lw=1.8, ls=:solid)
            Plots.plot!(p_comp_sec, periods, pi_agr_i, subplot=i,
                label=(i == 1 ? "Agr TFP (−10%)" : ""), color=:steelblue, lw=1.8, ls=:dash,
                title=short_names[i], titlefontsize=7,
                ylabel=(i % 3 == 1 ? "ann. pp" : ""))
            Plots.hline!(p_comp_sec, [0.0], subplot=i, color=:black, lw=0.5, ls=:dot, label="")
        end
        save_fig(p_comp_sec, "irf_comparison_sectoral_inflation.pdf")

    else
        @printf "  [Oil IRF data not found at %s — skipping comparison figures]\n" oil_irfs_path
        @printf "  Run oil_shock_analysis.jl first to generate comparison figures.\n"
    end

end  # _HAS_PLOTS


# =========================================================================== #
#  SAVE IRF DATA AS CSV                                                        #
# =========================================================================== #

df_irf_out = DataFrame(period=repeat(1:n_irf, outer=n_endo))
df_irf_out.variable = repeat(endo_names, inner=n_irf)
df_irf_out.value    = vec(irf_pct')
CSV.write(joinpath(TABLES_DIR, "agr_shock_irfs.csv"), df_irf_out)
@printf "\n  Saved: agr_shock_irfs.csv\n"


# =========================================================================== #
#  DONE                                                                        #
# =========================================================================== #

@printf "\n%s\n" repeat("=", 70)
@printf "  Agriculture TFP shock analysis (El Niño) complete.\n"
@printf "  Figures → %s\n" FIGURES_DIR
if overleaf_ok; @printf "  Overleaf → %s\n" OVERLEAF_DIR; end
@printf "  Tables  → %s\n" TABLES_DIR
@printf "%s\n\n" repeat("=", 70)

end  # function _main()

Base.invokelatest(_main)
