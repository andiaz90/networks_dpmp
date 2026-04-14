"""
main_SOE_gap.jl
===============
NK-IOSOE 12-sector model for Chile — Julia translation of main_SOE_gap.m

Workflow:
  1.  Select exercise (0–3)
  2.  Load sector-level data (Excel + CSV)
  3.  Set structural parameters
  4.  Configure shocks for the selected exercise
  5.  (Optional) load SMM estimates from smm_estimates.mat
  6.  Solve outer steady-state system (NLsolve)
  7.  Evaluate full steady-state quantities
  8.  Save params_val_ul.mat for Dynare to read
  9.  Run Dynare (via Octave — free MATLAB alternative)
  10. Load Dynare results; check BK conditions
  11. Compute rank correlations (Lyapunov-based std devs)
  12. Save results to <exercise_label>_results.mat
  13. Print moment-fit table

Usage:
  julia --project=. main_SOE_gap.jl

or from a Julia REPL:
  include("main_SOE_gap.jl")

REQUIREMENTS (all free):
  - Julia ≥ 1.9  (https://julialang.org/downloads/)
  - Julia packages in Project.toml  (] instantiate)
  - Octave ≥ 8   (https://octave.org/download)
  - Dynare ≥ 6   (https://dynare.org/download) — installed for Octave use
  - Data files: see DATA PATHS section below
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
using XLSX
using MAT          # for reading smm_estimates.mat / smm_best_so_far.mat
using StatsBase
using Dynare        # native Julia reimplementation — no MATLAB, no Octave needed

# =========================================================================== #
#  INCLUDE HELPERS                                                              #
# =========================================================================== #

SCRIPT_DIR = @__DIR__    # directory of this file

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))


# =========================================================================== #
#  EXERCISE SELECTOR                                                           #
#  0 = Baseline (all shocks)                                                  #
#  1 = Consumption preference shock only                                      #
#  2 = Manufacturing TFP shock only                                           #
#  3 = Monetary policy shock only                                             #
# =========================================================================== #

EXERCISE = 2    # <<< CHANGE THIS (0=Baseline matches SMM calibration)

exercise_labels = [
    "Baseline (all shocks)",
    "Exercise 1: Preference shock",
    "Exercise 2: TFP shock to Manufacturing",
    "Exercise 3: Monetary policy shock",
]

# Allow override via environment variable (used by SMM estimation scripts)
smm_ex_env = get(ENV, "SMM_EXERCISE", "")
smm_called = !isempty(smm_ex_env)
if smm_called
    EXERCISE = parse(Int, smm_ex_env)
end

@printf "\n%s\n" repeat("=", 60)
@printf "  NK-SOE 12-sector model for Chile\n"
if smm_called
    @printf "  [SMM mode] Exercise %d: %s\n" EXERCISE exercise_labels[EXERCISE+1]
else
    @printf "  %s\n" exercise_labels[EXERCISE+1]
end
@printf "%s\n\n" repeat("=", 60)

@assert EXERCISE in 0:3 "Invalid EXERCISE value. Must be 0, 1, 2, or 3."


# =========================================================================== #
#  DATA PATHS                                                                  #
#  Adjust DYNARE_MATLAB_PATH and DATA_DIR to match your installation.          #
# =========================================================================== #

# Directory with .mod files (relative to this script)
MOD_DIR = joinpath(SCRIPT_DIR, "mod")

# ---- Data files ----
# NOTE: 'Stata_to_excel_few_industries_chile.xls' must be saved as .xlsx
#       (Excel→Save As→.xlsx) or the variable name below adjusted.
#       CSV files can be used directly.

# Try to find data files in common locations (mirrors MATLAB path search)
DATA_CANDIDATES = [
    SCRIPT_DIR,
    joinpath(SCRIPT_DIR, "..", "modelo_chile"),   # sibling folder
    joinpath(SCRIPT_DIR, "..", "..", "Data"),
    get(ENV, "NKIOSOE_DATA_DIR", ""),
]

function find_file(candidates, fname)
    for d in candidates
        isempty(d) && continue
        p = joinpath(d, fname)
        isfile(p) && return p
    end
    return ""
end

fname_industries = "Stata_to_excel_few_industries_chile.xlsx"   # saved as xlsx
fname_io         = "IO_2021_chile.csv"
fname_fpa        = "fpa_vector_few_industries_chile.csv"

path_industries = find_file(DATA_CANDIDATES, fname_industries)
path_io         = find_file(DATA_CANDIDATES, fname_io)
path_fpa        = find_file(DATA_CANDIDATES, fname_fpa)

for (p, n) in [(path_industries, fname_industries),
               (path_io, fname_io),
               (path_fpa, fname_fpa)]
    isempty(p) && error("""
    Required data file not found: $n
    Searched in: $(join(filter(!isempty, DATA_CANDIDATES), ", "))
    Set the environment variable NKIOSOE_DATA_DIR to the folder containing the data files,
    or copy them next to this script.
    """)
end


# =========================================================================== #
#  READ DATA (12 sectors)                                                      #
# =========================================================================== #

nsec = 12

@printf "--- Loading data ---\n"

# ---- Sector-level Excel file ----
# Columns (1-indexed):
#   1=kk  2=BEAName  3=Nameshort  4=dp  5=dy  6=dl  7=dri
#   8=share  9=industrytype  10=spend_good  11=spend_serv
#   12=alpha  13=alpha_V  14=var_rho
wb = XLSX.readxlsx(path_industries)
ws = wb[1]   # first sheet; change index if needed

# Read rows 2:nsec+1 (skip header row 1)
names_vec  = String[ws[i, 3] for i in 2:nsec+1]
p_d        = Float64[ws[i, 4]  for i in 2:nsec+1]
y_d        = Float64[ws[i, 5]  for i in 2:nsec+1]
l_d        = Float64[ws[i, 6]  for i in 2:nsec+1]
ri_d       = Float64[ws[i, 7]  for i in 2:nsec+1]
spend_good = Float64[ws[i, 10] for i in 2:nsec+1]
spend_serv = Float64[ws[i, 11] for i in 2:nsec+1]
alpha      = Float64[ws[i, 12] for i in 2:nsec+1]
alpha_V    = Float64[ws[i, 13] for i in 2:nsec+1]
var_rho    = Float64[ws[i, 14] for i in 2:nsec+1]

# ---- IO matrix ----
betaio_df = CSV.read(path_io, DataFrame, header=false)
betaio    = Matrix{Float64}(betaio_df[1:nsec, 1:nsec])

# Normalize columns (each column j sums to 1 over supplying sectors)
col_sums  = sum(betaio, dims=1)
betax     = betaio ./ col_sums
modbeta   = betax'          # modbeta[i,j] = share of inputs sector i gets from sector j

# ---- Price adjustment frequency (Rotemberg κ) ----
fpa_df    = CSV.read(path_fpa, DataFrame, header=false)
theta_vec = vec(Matrix{Float64}(fpa_df))   # frequency of price adjustment

# Material / import shares
modalpha   = alpha       # intermediate material share
modalphaV  = alpha_V     # imported inputs share

@printf "  Loaded %d sectors from data files.\n\n" nsec


# =========================================================================== #
#  SET STRUCTURAL PARAMETERS FROM DATA                                         #
# =========================================================================== #

beta_val = 0.986        # quarterly discount factor (~5.7% annual discount rate)
epsilon  = 10.0         # variety elasticity (Ferrante et al. 2023)
gamma    = 2.0          # risk aversion / inverse IES
psi      = 1.0          # inverse Frisch elasticity
chi      = 1.0          # labor disutility weight

# Rotemberg adjustment costs (κ_i = θ(ε-1) / [(1-θ)(1-θβ)])
modkappa = theta_vec .* (epsilon - 1) ./ ((1 .- theta_vec) .* (1 .- theta_vec .* beta_val))

# Goods vs. services classification
goods    = spend_good .> spend_serv    # Bool vector
services = spend_serv .> spend_good    # Bool vector

# Consumption basket weights
modgammag = spend_good ./ sum(spend_good)  # within-goods sector shares γ^g_i
modgammas = spend_serv ./ sum(spend_serv)  # within-services sector shares γ^s_i

# Labor adjustment cost (set to zero in baseline — controlled by ilabcosts_val below)
modcl    = zeros(nsec)
modclneg = zeros(nsec)
modcm    = zeros(nsec)

# Elasticities (homogeneous across sectors)
modepsM  = fill(0.1, nsec)    # materials substitution elasticity ε^m
modepsY  = fill(0.8, nsec)    # production input substitution elasticity ε^Y

# Labor supply parameters
modpsil = fill(-1000.0, nsec)   # large negative = no sectoral adjustment cost
modpsim = fill(-1000.0, nsec)

# Monetary policy / Taylor rule
phi_val      = 2.5      # Taylor inflation coefficient φ_π
rhoi_val     = 0.6      # interest rate persistence
rhoirule_val = 0.74     # smoothing parameter

# Structural
gammaind_val = 0.0       # price indexation
ilabcosts_val = 0.1      # inverse of aggregate labor adjustment cost (1/Ψ^L)
ombar_val     = 0.57     # steady-state goods share in consumption basket

# SOE / foreign sector
Pistar_ss    = 1.00      # steady-state foreign inflation (zero net inflation)
Rworld_ss    = Pistar_ss / beta_val   # from Euler equation in SS

kappaV_val   = 1e13      # import price adjustment cost (very large ≈ fixed imports)
epsilonV_val = 1e13      # import demand elasticity (very large)
epsilonX_val = 1.0       # export demand elasticity
omegaX_val   = 1.0       # export demand scale

chii_b_val   = 0.001     # debt-risk premium elasticity χ_b
etastar_val  = 3.5       # foreign demand elasticity η*

xi_rstar_val    = 0.2    # world rate shock persistence
ystar_ss_val    = 1.0    # steady-state foreign output
PVstar_ss       = 1.0    # steady-state foreign price (normalized)
sigmaH_val      = 0.999  # Armington elasticity (home vs. foreign)

# Export and home-bias
modchiX   = fill(1/nsec, nsec)  # uniform export shares (placeholder)
modvarrho = var_rho              # home-bias parameters ϱ_i from data
modA      = ones(nsec)           # TFP normalized to 1 in SS

# Default shock parameters (overridden by exercise below)
rho_val          = 0.1
sigma_L_agg_val  = 0.0
sigma_L_het      = zeros(nsec)

# Import price (PVstar) shock
rho_pvstar_val   = 0.9
sigma_pvstar_val = 0.03

# Preference (xi) shock
rho_xi_val   = 0.80
sigma_xi_val = 0.005


# =========================================================================== #
#  SHOCK PARAMETERS BY EXERCISE                                                #
# =========================================================================== #

if EXERCISE == 0
    # --- Baseline: all shocks active ---
    sigma_i_val    = 0.001
    rho_om1_val    = 0.1;   rho_om2_val = 0.0
    sigma_om_val   = 0.0001
    rho_tfp1_val   = 0.5;   rho_tfp2_val = 0.0
    isigma_tfp_val = fill(0.05, nsec)

elseif EXERCISE == 1
    # --- Exercise 1: preference shock only ---
    sigma_i_val    = 0.0
    rho_om1_val    = 0.95;  rho_om2_val = 0.0
    sigma_om_val   = 0.01
    rho_tfp1_val   = 0.5;   rho_tfp2_val = 0.0
    isigma_tfp_val = zeros(nsec)

elseif EXERCISE == 2
    # --- Exercise 2: manufacturing TFP shock only ---
    mfg = 3   # sector 3 = Manufacturing
    sigma_i_val    = 0.0
    rho_om1_val    = 0.1;   rho_om2_val = 0.0
    sigma_om_val   = 0.0
    rho_tfp1_val   = 0.95;  rho_tfp2_val = 0.0
    isigma_tfp_val = zeros(nsec)
    isigma_tfp_val[mfg] = 0.01

elseif EXERCISE == 3
    # --- Exercise 3: monetary policy shock only ---
    sigma_i_val    = 0.01
    rho_om1_val    = 0.1;   rho_om2_val = 0.0
    sigma_om_val   = 0.0
    rho_tfp1_val   = 0.5;   rho_tfp2_val = 0.0
    isigma_tfp_val = zeros(nsec)
end

# Shock variance indicators for Dynare shocks block
shock_eps_om_val     = Float64(sigma_om_val > 0)
shock_eps_i_val      = Float64(sigma_i_val  > 0)
shock_eps_pvstar_val = Float64(EXERCISE == 0 && sigma_pvstar_val > 0)
shock_eps_xi_val     = Float64(EXERCISE == 0 && sigma_xi_val > 0)
shock_epsA_val       = ones(nsec)


# =========================================================================== #
#  LOAD SMM ESTIMATES (override defaults when available)                       #
# =========================================================================== #

smm_est_file  = joinpath(SCRIPT_DIR, "..", "modelo_chile", "smm_estimates.mat")
smm_ckpt_file = joinpath(SCRIPT_DIR, "..", "modelo_chile", "smm_best_so_far.mat")

param_names = ["ilabcosts", "epsY", "epsM", "kappaV", "rho_om1", "sigma_om",
               "rho_tfp1", "isigma_tfp(1)", "rho_pvstar", "sigma_pvstar"]

smm_param_source = "hard-coded defaults"

if isfile(smm_est_file)
    est = matread(smm_est_file)
    ilabcosts_val    = est["ilabcosts_val"]
    modepsY          = fill(est["modepsY"] isa AbstractArray ? est["modepsY"][1] : est["modepsY"], nsec)
    modepsM          = fill(est["modepsM"] isa AbstractArray ? est["modepsM"][1] : est["modepsM"], nsec)
    kappaV_val       = est["kappaV_val"]
    rho_om1_val      = est["rho_om1_val"]
    sigma_om_val     = est["sigma_om_val"]
    rho_tfp1_val     = est["rho_tfp1_val"]
    isigma_tfp_val   = vec(est["isigma_tfp_val"])
    rho_pvstar_val   = est["rho_pvstar_val"]
    sigma_pvstar_val = est["sigma_pvstar_val"]
    haskey(est, "rho_xi_val")   && (rho_xi_val   = est["rho_xi_val"])
    haskey(est, "sigma_xi_val") && (sigma_xi_val = est["sigma_xi_val"])
    smm_param_source = "smm_estimates.mat"

elseif isfile(smm_ckpt_file)
    ckpt  = matread(smm_ckpt_file)
    tb    = ckpt["smm_best_so_far"]
    θ_best = vec(tb["theta_best"])
    ilabcosts_val    = θ_best[1]
    modepsY          = fill(θ_best[2], nsec)
    modepsM          = fill(θ_best[3], nsec)
    kappaV_val       = exp(θ_best[4])
    rho_om1_val      = θ_best[5]
    sigma_om_val     = θ_best[6]
    rho_tfp1_val     = θ_best[7]
    isigma_tfp_val   = θ_best[8:19]
    rho_pvstar_val   = θ_best[20]
    sigma_pvstar_val = θ_best[21]
    length(θ_best) >= 22 && (rho_xi_val   = θ_best[22])
    length(θ_best) >= 23 && (sigma_xi_val = θ_best[23])
    smm_param_source = @sprintf "smm_best_so_far.mat (obj=%.6f)" tb["obj_best"]
end

@printf "--- Parameters (%s) ---\n" smm_param_source
param_vals = [ilabcosts_val, modepsY[1], modepsM[1], kappaV_val,
              rho_om1_val, sigma_om_val, rho_tfp1_val, isigma_tfp_val[1],
              rho_pvstar_val, sigma_pvstar_val]
for (nm, vl) in zip(param_names, param_vals)
    @printf "  %-18s %g\n" nm vl
end
@printf "\n"


# =========================================================================== #
#  STEADY STATE                                                                #
# =========================================================================== #

@printf "--- Steady State ---\n"

# Initial guesses (in levels → must be positive)
pHvec_guess = ones(nsec)
w_guess     = 1.0
Q_guess     = 1.0
C_guess     = 1.0
x_guess     = [pHvec_guess; w_guess; Q_guess; C_guess]

# Pack calibration parameters for the SS function
sigmaH    = sigmaH_val
varrho_val = modvarrho
gammag_vec = modgammag
gammas_vec = modgammas
om_g       = ombar_val
om_s       = 1 - ombar_val
chiX_vec   = modchiX
omegaX     = omegaX_val
etastar    = etastar_val
Ystar      = ystar_ss_val
alpha_vec  = modalpha
alphaV_vec = modalphaV
beta_mat   = modbeta
epsY_vec   = modepsY
epsM_vec   = modepsM
A_vec      = modA

# Load trade-balance target from data_moments_chile.mat if available
tb_target = 0.02
dm_file_tb = joinpath(SCRIPT_DIR, "..", "modelo_chile", "data_moments_chile.mat")
if isfile(dm_file_tb)
    try
        tmp_dm = matread(dm_file_tb)
        if haskey(tmp_dm, "dm_chile")
            dm_ch = tmp_dm["dm_chile"]
            if haskey(dm_ch, "d_TBGDP") && isfinite(dm_ch["d_TBGDP"])
                tb_target = dm_ch["d_TBGDP"]
            end
        end
    catch
    end
end

ss_result = nlsolve(
    (F, x) -> F .= steady_ntwsoe(
        x, PVstar_ss, epsilon, varrho_val, sigmaH,
        gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar,
        alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec,
        gamma, chi, psi, A_vec, tb_target
    ),
    x_guess;
    ftol = 1e-14,
    show_trace = !smm_called,
    method = :trust_region,
)

if !converged(ss_result)
    @warn "Outer steady-state solver did not converge (residual=$(norm(ss_result.residual, Inf))). Proceeding anyway — results may be inaccurate."
end

pH_ss = ss_result.zero[1:nsec]
w_ss  = ss_result.zero[nsec+1]
Q_ss  = ss_result.zero[nsec+2]
C_ss  = ss_result.zero[nsec+3]


# =========================================================================== #
#  EVALUATE FULL STEADY STATE                                                  #
# =========================================================================== #

PL_ss   = fill(w_ss, nsec)
PV_ss   = Q_ss * PVstar_ss

MCi_ss  = (epsilon - 1) / epsilon .* pH_ss
PMi_ss  = (beta_mat * (pH_ss .^ (1 .- epsM_vec))) .^ (1 ./ (1 .- epsM_vec))

P_ss    = (varrho_val .^ sigmaH .* pH_ss .^ (1 - sigmaH)
         .+ (1 .- varrho_val) .^ sigmaH .* PV_ss .^ (1 - sigmaH)) .^ (1/(1-sigmaH))

p_g_ss  = prod(P_ss .^ gammag_vec)
p_s_ss  = prod(P_ss .^ gammas_vec)

C_g_ss  = om_g * C_ss / p_g_ss
C_s_ss  = om_s * C_ss / p_s_ss

C_gi_ss = gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss
C_si_ss = gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss

CHg_ss  = varrho_val .^ sigmaH .* (pH_ss ./ P_ss) .^ (-sigmaH) .* C_gi_ss
CHs_ss  = varrho_val .^ sigmaH .* (pH_ss ./ P_ss) .^ (-sigmaH) .* C_si_ss
CFg_ss  = (1 .- varrho_val) .^ sigmaH .* (PV_ss ./ P_ss) .^ (-sigmaH) .* C_gi_ss
CFs_ss  = (1 .- varrho_val) .^ sigmaH .* (PV_ss ./ P_ss) .^ (-sigmaH) .* C_si_ss

CHi_ss  = CHg_ss .+ CHs_ss
CFi_ss  = CFg_ss .+ CFs_ss

PX_ss   = prod(pH_ss .^ chiX_vec)
X_ss    = omegaX * (PX_ss / Q_ss)^(-etastar) * Ystar
Xi_ss   = chiX_vec .* X_ss .* PX_ss ./ pH_ss

# Inner solve for production inputs
ig_ss = zeros(nsec)
for i in 1:nsec, j in 1:nsec
    ig_ss[i] += beta_mat[j, i]
end
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
    (F, x) -> steady_ntwsoe_system!(
        F, x, alpha_vec, alphaV_vec, beta_mat,
        MCi_ss, PMi_ss, PL_ss, PV_ss, CHi_ss, Xi_ss,
        epsY_vec, epsM_vec, A_vec, pH_ss
    ),
    [M_init; L_init; Vi_init; Yi_init];
    ftol = 1e-10, show_trace = false, method = :trust_region,
)

M_ss   = inner_sol.zero[1:nsec]
L_ss   = inner_sol.zero[nsec+1:2*nsec]
Vi_ss  = inner_sol.zero[2*nsec+1:3*nsec]
Yi_ss  = inner_sol.zero[3*nsec+1:4*nsec]

V_ss          = sum(Vi_ss)
CF_ss         = sum(CFi_ss)
mkupV         = 1.0
IMP_tot_ss    = mkupV * (V_ss + CF_ss)
TB_ss         = PX_ss * X_ss - PV_ss * IMP_tot_ss
GDP_ss        = C_ss + TB_ss
N_ss          = sum(L_ss)
Y_ss          = sum(Yi_ss)
r_star_ss     = Rworld_ss
pi_ss         = 1.0
r_ss          = 1 / beta_val
tbgdp         = TB_ss / GDP_ss
Bstar_ss      = -TB_ss / (Q_ss * (1 - r_star_ss / Pistar_ss))
bbar_val      = Q_ss * Bstar_ss / GDP_ss

ygdp_target   = 2.0
fdebt_target  = 0.70
ygdp_ss       = Y_ss / GDP_ss
fdebt_ss      = Q_ss * abs(Bstar_ss) / GDP_ss

# Aggregates for Dynare
M_tot_ss      = sum(M_ss)
Y_tot_ss      = sum(Yi_ss)
VA_ss_val     = sum(Yi_ss .- M_ss)
Ctotg_ss_val  = sum(gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss)
Ctots_ss_val  = sum(gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss)
Ctot_ss_val   = Ctotg_ss_val + Ctots_ss_val
IMP_ss_val    = IMP_tot_ss

Pistar_ss_val  = Pistar_ss
Rworld_ss_val  = Rworld_ss
PVstar_ss_val  = PVstar_ss

# Foreign demand shock params (kept for smm_estimates.mat compatibility)
rho_psi_val   = 0.5
sigma_psi_val = 0.001

@printf "  GDP=%.4f  TB/GDP=%.4f  X=%.4f  IMP=%.4f  bbar=%.4f\n" GDP_ss tbgdp X_ss IMP_tot_ss bbar_val
@printf "  TB/GDP target: %.4f  |  world rate (ann.): %.2f%%\n" tb_target 400*(Rworld_ss-1)
@printf "  Bstar=%.4f  Q*Bstar=%.4f  Q*Bstar/GDP=%.3f (%.1f%%)\n" Bstar_ss (Q_ss*Bstar_ss) (Q_ss*Bstar_ss/GDP_ss) (100*Q_ss*Bstar_ss/GDP_ss)
@printf "  Gross output/GDP : %.3f  (IMF target: %.1f)%s\n" ygdp_ss ygdp_target (ternary_str(abs(ygdp_ss-ygdp_target)>0.20, "  << off target", ""))
@printf "  Foreign debt/GDP : %.3f  (IMF target: %.2f)%s\n\n" fdebt_ss fdebt_target (ternary_str(abs(fdebt_ss-fdebt_target)>0.10, "  << off target", ""))


# =========================================================================== #
#  WRITE params_jl.mod  (read by NK_SOE_lev_gap2.mod via @#include)          #
#  Dynare.jl processes this as a Dynare parameter file — no MATLAB/Octave.   #
# =========================================================================== #

params_nt = (
    nsec           = nsec,
    # Shock scalars
    sigma_i_val    = sigma_i_val,
    sigma_L_agg_val= sigma_L_agg_val,
    sigma_om_val   = sigma_om_val,
    ilabcosts_val  = ilabcosts_val,
    gamma_val      = gamma,
    beta_val       = beta_val,
    phi_val        = phi_val,
    rho_val        = rho_val,
    rho_om1_val    = rho_om1_val,
    rho_om2_val    = rho_om2_val,
    rho_tfp1_val   = rho_tfp1_val,
    rho_tfp2_val   = rho_tfp2_val,
    rhoi_val       = rhoi_val,
    rhoirule_val   = rhoirule_val,
    ombar_val      = ombar_val,
    chii_b_val     = chii_b_val,
    bbar_val       = bbar_val,
    epsilonX_val   = epsilonX_val,
    omegaX_val     = omegaX_val,
    ystar_ss_val   = ystar_ss_val,
    etastar_val    = etastar_val,
    epsilonV_val   = epsilonV_val,
    kappaV_val     = kappaV_val,
    sigmaH_val     = sigmaH_val,
    Pistar_ss_val  = Pistar_ss_val,
    PVstar_ss_val  = PVstar_ss_val,
    Rworld_ss_val  = Rworld_ss_val,
    rho_pvstar_val = rho_pvstar_val,
    sigma_pvstar_val=sigma_pvstar_val,
    rho_xi_val     = rho_xi_val,
    sigma_xi_val   = sigma_xi_val,
    # Steady state scalars
    w_ss           = w_ss,
    C_ss           = C_ss,
    GDP_ss         = GDP_ss,
    N_ss           = N_ss,
    p_s_ss         = p_s_ss,
    p_g_ss         = p_g_ss,
    C_s_ss         = C_s_ss,
    C_g_ss         = C_g_ss,
    Bstar_ss       = Bstar_ss,
    IMP_ss_val     = IMP_ss_val,
    Ctot_ss_val    = Ctot_ss_val,
    Ctotg_ss_val   = Ctotg_ss_val,
    Ctots_ss_val   = Ctots_ss_val,
    VA_ss_val      = VA_ss_val,
    M_tot_ss       = M_tot_ss,
    Y_tot_ss       = Y_tot_ss,
    # Shock flags
    shock_eps_om_val     = shock_eps_om_val,
    shock_eps_i_val      = shock_eps_i_val,
    shock_eps_pvstar_val = shock_eps_pvstar_val,
    shock_eps_xi_val     = shock_eps_xi_val,
    shock_epsA_val       = shock_epsA_val,
    # Sectoral vectors
    modgammag    = modgammag,
    modgammas    = modgammas,
    modalpha     = modalpha,
    modalphaV    = modalphaV,
    modepsY      = modepsY,
    modepsM      = modepsM,
    modkappa     = modkappa,
    goods        = goods,
    services     = services,
    modcl        = modcl,
    modclneg     = modclneg,
    modcm        = modcm,
    modchiX      = modchiX,
    modvarrho    = modvarrho,
    isigma_tfp_val = isigma_tfp_val,
    PL_ss        = PL_ss,
    CFg_ss       = CFg_ss,
    CFs_ss       = CFs_ss,
    CHg_ss       = CHg_ss,
    CHs_ss       = CHs_ss,
    Vi_ss        = Vi_ss,
    pH_ss        = pH_ss,
    MCi_ss       = MCi_ss,
    Yi_ss        = Yi_ss,
    L_ss         = L_ss,
    C_gi_ss      = C_gi_ss,
    C_si_ss      = C_si_ss,
    P_ss         = P_ss,
    PMi_ss       = PMi_ss,
    M_ss         = M_ss,
    modbeta      = modbeta,
)

params_mod_path = write_params_mod(MOD_DIR, params_nt)
@printf "--- Wrote params_jl.mod to %s ---\n\n" params_mod_path


# =========================================================================== #
#  RUN DYNARE.jl  — pure Julia, no MATLAB, no Octave                         #
# =========================================================================== #

@printf "--- Running Dynare.jl ---\n"
@printf "  Model dir : %s\n\n" MOD_DIR

# @dynare must be called with the path to the .mod file (without extension).
# It reads params_jl.mod via @#include before solving the model.
context = cd(MOD_DIR) do
    @dynare "NK_SOE_lev_gap2"
end

@printf "\n--- Dynare.jl completed ---\n\n"


# =========================================================================== #
#  EXTRACT RESULTS FROM DYNARE.jl CONTEXT                                     #
#                                                                              #
#  Dynare.jl stores everything in the `context` object:                       #
#    Steady state    : context.results.model_results[1].trends                #
#                      .endogenous_steady_state                                #
#    Decision rule   : context.results.model_results[1].linearrationalexpectations #
#                      .g1_1  (state feedback ≈ ghx)                          #
#                      .g1_2  (shock impact  ≈ ghu)                           #
#    Simulations     : context.results.model_results[1].simulations[1].data   #
#    Variable names  : get_endogenous(context.symboltable)                     #
#    Shock cov matrix: context.models[1].Sigma_e                               #
# =========================================================================== #

mr        = context.results.model_results[1]
endo_names = Dynare.get_endogenous(context.symboltable)

# Steady state vector (one value per endogenous variable, declaration order)
ss_vec    = mr.trends.endogenous_steady_state

# Decision rule (first-order approximation)
# g1_1 : n_endo × n_states  (≈ oo_.dr.ghx in MATLAB)
# g1_2 : n_endo × n_shocks  (≈ oo_.dr.ghu in MATLAB)
lre = mr.linearrationalexpectations
ghx_jl = lre.g1_1
ghu_jl = lre.g1_2
Σe_jl  = context.models[1].Sigma_e

ss_ok  = !any(isnan, ss_vec) && !any(isinf, ss_vec)
sim_ok = !isempty(mr.simulations)

@printf "--- Dynare.jl results ---\n"
@printf "  Steady state : %s\n" (ss_ok  ? "OK" : "FAILED")
@printf "  Simulation   : %s\n\n" (sim_ok ? "OK" : "FAILED")

(!ss_ok || !sim_ok) && (@printf "Model failed – aborting.\n"; exit(1))

# Build variable-name → index map
endo_idx = Dict(nm => i for (i, nm) in enumerate(endo_names))

# Simulation time series (AxisArrayTable → named-column table)
sim_data  = mr.simulations[1].data   # AxisArrayTable, rows=periods, cols=variables
# Access as a matrix: columns correspond to endo_names order
sim_matrix = Matrix(sim_data)        # periods × n_endo


# =========================================================================== #
#  RANK CORRELATIONS (Lyapunov-based, mirrors smm_model_moments.m)            #
# =========================================================================== #

@printf "--- Rank correlations (model vs data, Spearman) ---\n"

std_Y_m  = fill(NaN, nsec)
std_PH_m = fill(NaN, nsec)
std_L_m  = fill(NaN, nsec)

rc_lyap_ok = false
try
    # State-variable indices in the DR ordering
    n_state = size(ghx_jl, 2)
    n_endo  = length(endo_names)

    # Lyapunov: unconditional variance of all endogenous variables
    # P_state = ghx * P_state * ghx' + ghu * Σe * ghu'
    ghx_s = ghx_jl[1:n_state, :]   # state rows of ghx
    ghu_s = ghu_jl[1:n_state, :]   # state rows of ghu
    P_st  = local_dlyap(ghx_s, ghu_s * Σe_jl * ghu_s')
    Γ_rc  = ghx_jl * P_st * ghx_jl' + ghu_jl * Σe_jl * ghu_jl'
    Γ_rc  = (Γ_rc + Γ_rc') / 2

    for i in 1:nsec
        for (kv, vn) in enumerate(["Y_$(i)", "PH_$(i)", "L_$(i)"])
            dr_idx = get(endo_idx, vn, nothing)
            dr_idx === nothing && continue
            ss_val = abs(ss_vec[dr_idx])
            ss_val < 1e-12 && (ss_val = 1.0)
            pstd = sqrt(max(Γ_rc[dr_idx, dr_idx], 0.0)) / ss_val
            kv == 1 && (std_Y_m[i]  = pstd)
            kv == 2 && (std_PH_m[i] = pstd)
            kv == 3 && (std_L_m[i]  = pstd)
        end
    end
    rc_lyap_ok = true
catch e
    @printf "  [rank corr] Lyapunov failed (%s); using simulation std devs.\n" string(e)
end

# Fallback: simulation std devs
if !rc_lyap_ok
    for i in 1:nsec
        for (kv, vn) in enumerate(["Y_$(i)", "PH_$(i)", "L_$(i)"])
            col = get(endo_idx, vn, nothing)
            col === nothing && continue
            ts     = sim_matrix[:, col]
            ss_val = abs(ss_vec[col])
            ss_val < 1e-12 && (ss_val = 1.0)
            pstd   = std(ts) / ss_val
            kv == 1 && (std_Y_m[i]  = pstd)
            kv == 2 && (std_PH_m[i] = pstd)
            kv == 3 && (std_L_m[i]  = pstd)
        end
    end
end

# Compute Spearman rank correlations
valid_y = isfinite.(std_Y_m)  .& isfinite.(y_d)
valid_p = isfinite.(std_PH_m) .& isfinite.(p_d)
valid_l = isfinite.(std_L_m)  .& isfinite.(l_d)

rho_y = sum(valid_y) >= 3 ? safe_spearman(std_Y_m[valid_y],  y_d[valid_y])  : NaN
rho_p = sum(valid_p) >= 3 ? safe_spearman(std_PH_m[valid_p], p_d[valid_p])  : NaN
rho_l = sum(valid_l) >= 3 ? safe_spearman(std_L_m[valid_l],  l_d[valid_l])  : NaN

@printf "  %-30s  %8s\n" "Variable" "Spearman r"
@printf "  %s\n" repeat("-", 42)
@printf "  %-30s  %8.4f\n" "Sectoral output (Y)" rho_y
@printf "  %-30s  %8.4f\n" "Prices (PH)"         rho_p
@printf "  %-30s  %8.4f\n" "Labor (L)"           rho_l
(rho_y == 0 || rho_p == 0 || rho_l == 0) &&
    @printf "  (0 = model std devs are uniform across sectors for this exercise)\n"
@printf "\n"

rank_corr = (rho_output=rho_y, rho_price=rho_p, rho_labor=rho_l)


# =========================================================================== #
#  SAVE RESULTS                                                                #
# =========================================================================== #

output_files = ["model_output_IOSOE", "model_output_IOSOE_ex1",
                "model_output_IOSOE_ex2", "model_output_IOSOE_ex3"]
output_path  = joinpath(MOD_DIR, "$(output_files[EXERCISE+1]).mat")

results_out = Dict{String, Any}(
    "EXERCISE"     => EXERCISE,
    "nsec"         => nsec,
    "pH_ss"        => pH_ss,
    "w_ss"         => w_ss,
    "Q_ss"         => Q_ss,
    "C_ss"         => C_ss,
    "Yi_ss"        => Yi_ss,
    "L_ss"         => L_ss,
    "M_ss"         => M_ss,
    "Vi_ss"        => Vi_ss,
    "GDP_ss"       => GDP_ss,
    "TB_ss"        => TB_ss,
    "Bstar_ss"     => Bstar_ss,
    "rank_corr_Y"  => rho_y,
    "rank_corr_P"  => rho_p,
    "rank_corr_L"  => rho_l,
    "std_Y_m"      => std_Y_m,
    "std_PH_m"     => std_PH_m,
    "std_L_m"      => std_L_m,
    "y_d"          => y_d,
    "p_d"          => p_d,
    "l_d"          => l_d,
)
matwrite(output_path, results_out)
@printf "--- Results saved to: %s ---\n\n" output_path


# =========================================================================== #
#  DONE                                                                        #
# =========================================================================== #

@printf "%s\n" repeat("=", 60)
@printf "  Done: %s\n" exercise_labels[EXERCISE+1]
@printf "%s\n\n" repeat("=", 60)
