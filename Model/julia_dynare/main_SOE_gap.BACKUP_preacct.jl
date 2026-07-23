"""
main_SOE_gap.jl
===============
NK-IOSOE 12-sector model for Chile — Julia/Dynare.jl translation of main_SOE_gap.m
No MATLAB. No Octave. No licenses.

Workflow:
  1.  Select exercise (0–3)
  2.  Load sector-level data (Excel + CSV from Data/)
  3.  Set structural parameters
  4.  Configure shocks for the selected exercise
  5.  (Optional) load SMM estimates from modelo_chile/smm_estimates.mat
  6.  Solve outer steady-state system (NLsolve)
  7.  Evaluate full steady-state quantities
  8.  Write mod/params_jl.mod (read by NK_SOE_lev_gap2.mod via @#include)
  9.  Run Dynare.jl: context = @dynare "NK_SOE_lev_gap2"
  10. Extract results from context object
  11. Compute rank correlations (Lyapunov-based std devs)
  12. Save results to julia_dynare/<exercise>_results.mat

Usage:
  julia --project=. main_SOE_gap.jl

REQUIREMENTS (all free):
  - Julia ≥ 1.9  (https://julialang.org/downloads/)
  - Julia packages: run  julia --project=. -e "import Pkg; Pkg.instantiate()"
  - Data files in Data/ folder (see DATA PATHS section below)
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
using Dynare        # native Julia reimplementation — no MATLAB, no Octave needed

# =========================================================================== #
#  INCLUDE HELPERS                                                              #
# =========================================================================== #

SCRIPT_DIR = @__DIR__    # directory of this file

include(joinpath(SCRIPT_DIR, "steady_ntwsoe_system.jl"))
include(joinpath(SCRIPT_DIR, "steady_ntwsoe.jl"))
include(joinpath(SCRIPT_DIR, "utils.jl"))
include(joinpath(SCRIPT_DIR, "figs_SOE_gap.jl"))    # generate_figures()
include(joinpath(SCRIPT_DIR, "plot_scripts.jl"))    # run_all_plots() — all MATLAB plot scripts


# =========================================================================== #
#  EXERCISE SELECTOR                                                           #
#  0 = Baseline (all shocks)                                                  #
#  1 = Consumption preference shock only                                      #
#  2 = Manufacturing TFP shock only                                           #
#  3 = Monetary policy shock only                                             #
# =========================================================================== #

# Wrap in _main() so that:
#   1. Julia 1.12 world-age issues are avoided inside a function scope.
#   2. Soft-scope ambiguities go away inside a function.
#   3. We can use `return` for early exit.
function _main()

# <<<  CHANGE THIS to select exercise  >>>
# 0 = Baseline (all shocks)          1 = Preference shock only
# 2 = Manufacturing TFP shock only   3 = Monetary policy shock only
EXERCISE = 2

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

# Short tag for filenames — defined early so all sections can use it
tag = ["baseline", "ex1_pref", "ex2_mfg", "ex3_mp"][EXERCISE + 1]


# =========================================================================== #
#  PATHS                                                                       #
# =========================================================================== #

# Directory with .mod files — always next to this script
MOD_DIR = joinpath(SCRIPT_DIR, "mod")

# Repo root is two levels up:  julia_dynare/ → Model/ → Networks-DPMP/
REPO_ROOT = abspath(joinpath(SCRIPT_DIR, "..", ".."))

# SMM estimates live in the sibling modelo_chile/ folder
MODELO_DIR = abspath(joinpath(SCRIPT_DIR, "..", "modelo_chile"))

# Output directories inside julia_dynare/ (created automatically)
FIGURES_DIR      = joinpath(SCRIPT_DIR, "figures")          # parent
FIGURES_IRF_DIR  = joinpath(FIGURES_DIR, "irfs")           # Dynare-style IRF plots
FIGURES_EX_DIR   = joinpath(FIGURES_DIR, "exercises")      # exercise analysis plots
TABLES_DIR       = joinpath(SCRIPT_DIR,  "tables")
mkpath(FIGURES_IRF_DIR)
mkpath(FIGURES_EX_DIR)
mkpath(TABLES_DIR)

# ---- Data file search ----
DATA_DIR = joinpath(REPO_ROOT, "Data")

DATA_CANDIDATES = filter!(!isempty, [
    DATA_DIR,              # Networks-DPMP/Data/     ← primary
    MODELO_DIR,            # Model/modelo_chile/      (some files live here too)
    SCRIPT_DIR,            # same folder as this script
    get(ENV, "NKIOSOE_DATA_DIR", ""),   # user override
])

# find_file: search candidates for any of the given filenames
function find_file(candidates, fnames...)
    for d in candidates, fname in fnames
        p = joinpath(d, fname)
        isfile(p) && return p
    end
    return ""
end

# ---- All data files are CSV (no XLS / no MAT) ----
# Run bootstrap_csv.jl once if sector_calibration.csv or sectoral_moments.csv
# are missing (they are generated from the MATLAB-produced .mat files).
path_cal     = find_file(DATA_CANDIDATES, "sector_calibration.csv")
path_io      = find_file(DATA_CANDIDATES, "IO_2021_chile.csv")
path_fpa     = find_file(DATA_CANDIDATES, "fpa_vector_few_industries_chile.csv")
path_sec_mom = find_file(DATA_CANDIDATES, "sectoral_moments.csv")
path_agg_mom = find_file(DATA_CANDIDATES, "aggregate_moments.csv")

@printf "  Data root : %s\n" DATA_DIR
for (p, label) in [
    (path_cal,     "sector_calibration.csv"),
    (path_io,      "IO_2021_chile.csv"),
    (path_fpa,     "fpa_vector CSV"),
    (path_sec_mom, "sectoral_moments.csv"),
    (path_agg_mom, "aggregate_moments.csv"),
]
    @printf "  %-30s : %s\n" label (isempty(p) ? "NOT FOUND" : "OK")
end
println()

missing_files = String[]
isempty(path_cal)     && push!(missing_files, "sector_calibration.csv")
isempty(path_io)      && push!(missing_files, "IO_2021_chile.csv")
isempty(path_fpa)     && push!(missing_files, "fpa_vector_few_industries_chile.csv")

if !isempty(missing_files)
    error("""
    Required CSV file(s) not found: $(join(missing_files, ", "))
    Searched in: $(join(DATA_CANDIDATES, ", "))

    If sector_calibration.csv is missing, run once from julia_dynare/:
        julia --project=. bootstrap_csv.jl
    This generates the CSV files from the existing MATLAB .mat outputs.
    """)
end


# =========================================================================== #
#  READ DATA (12 sectors, all from CSV)                                        #
# =========================================================================== #

nsec = 12

@printf "--- Loading data ---\n"

# ---- Sector calibration CSV ----
# Columns: sector, name, alpha, alpha_V, var_rho, spend_good, spend_serv,
#          kappa, is_goods
cal_df     = CSV.read(path_cal, DataFrame)
names_vec  = cal_df.name
alpha      = Float64.(cal_df.alpha)
alpha_V    = Float64.(cal_df.alpha_V)
var_rho    = Float64.(cal_df.var_rho)
spend_good = Float64.(cal_df.spend_good)
spend_serv = Float64.(cal_df.spend_serv)

# ---- Data moments from CSV (y_d, p_d, l_d for rank correlations) ----
if !isempty(path_sec_mom)
    sec_mom = CSV.read(path_sec_mom, DataFrame)
    y_d = Float64.(sec_mom.std_Y)
    p_d = Float64.(sec_mom.std_PH)
    l_d = Float64.(sec_mom.std_L)
    c_d = hasproperty(sec_mom, :corr_YPH) ? Float64.(sec_mom.corr_YPH) : fill(NaN, nsec)
else
    @printf "  WARNING: sectoral_moments.csv not found — rank correlations will use zeros.\n"
    y_d = zeros(nsec); p_d = zeros(nsec); l_d = zeros(nsec); c_d = fill(NaN, nsec)
end

# ---- IO matrix ----
betaio_df = CSV.read(path_io, DataFrame, header=false)
betaio    = Matrix{Float64}(betaio_df[1:nsec, 1:nsec])

col_sums  = sum(betaio, dims=1)
betax     = betaio ./ col_sums
modbeta   = Matrix(betax')   # modbeta[i,j] = share of inputs sector i gets from sector j

# ---- Price adjustment frequency (Rotemberg κ) ----
fpa_df    = CSV.read(path_fpa, DataFrame, header=false)
theta_vec = vec(Matrix{Float64}(fpa_df))

# Material / import shares
modalpha  = alpha
modalphaV = alpha_V

@printf "  Loaded %d sectors from data files.\n\n" nsec


# =========================================================================== #
#  SET STRUCTURAL PARAMETERS FROM DATA                                         #
# =========================================================================== #

beta_val = 0.986        # quarterly discount factor (~5.7% annual discount rate)
epsilon  = 10.0         # variety elasticity (Ferrante et al. 2023)
gamma    = 2.0          # risk aversion / inverse IES
psi      = parse(Float64, get(ENV, "SMM_PSI", "1.5"))  # inverse Frisch (XMAS posterior ≈1.5, García et al. WP833); env SMM_PSI
chi      = 1.0          # labor disutility weight

# Rotemberg adjustment costs (κ_i = θ(ε-1) / [(1-θ)(1-θβ)])
# theta_vec is the FREQUENCY of price adjustment (fraction of firms that reset each
# quarter). Calvo stickiness = probability of NOT adjusting = 1 - theta_vec.
# theta_vec = MONTHLY frequency of price adjustment. Chilean vector built by
# Data/build_fpa_vector_chile.py: Pastén et al. (2020) US PPI monthly freqs
# rescaled US→Chile, anchored to Albagli et al. (2026) regular-PPI manufacturing
# frequency (0.252). The model is QUARTERLY, so time-aggregate to a quarterly
# frequency before the Calvo→Rotemberg map: θ_q = 1 - (1-θ_m)^3 (constant hazard).
# Default ON; set SMM_FPA_QUARTERLY=0 to restore the legacy (monthly-as-quarterly).
theta_q  = get(ENV, "SMM_FPA_QUARTERLY", "1") == "1" ? (1 .- (1 .- theta_vec).^3) : theta_vec
stick    = 1 .- theta_q
modkappa = stick .* (epsilon - 1) ./ ((1 .- stick) .* (1 .- stick .* beta_val))

# Goods vs. services classification
goods    = spend_good .> spend_serv    # Bool vector
services = spend_serv .> spend_good    # Bool vector

# Consumption basket weights
modgammag = spend_good ./ sum(spend_good)  # within-goods sector shares γ^g_i
modgammas = spend_serv ./ sum(spend_serv)  # within-services sector shares γ^s_i

# Labor adjustment cost (set to zero in baseline — controlled by ilabcosts_val below)
modcl    = fill(parse(Float64, get(ENV, "SMM_CL", "0.0")), nsec)  # sectoral labor reallocation cost (env SMM_CL; 0 = free reallocation, the current baseline)
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

chii_b_val   = 0.0024    # debt-risk premium elasticity χ_b — XMAS posterior mean
                         # (García et al., BCCh; 100ψ = 0.24 [0.18, 0.30], Bayesian
                         # with Chile EMBIG observable; NFA/quarterly-GDP ratio,
                         # same units as here). Was 0.001 (hand-set) pre-2026-07-08.
etastar_val  = 3.5       # foreign demand elasticity η*
kappaw_val   = 115.0     # Rotemberg wage adj. cost (≈4q Calvo at epsw=10);
                         # overridden by smm_estimates.csv when available;
                         # 0 = flexible wages (nests pre-2026-07 model)

xi_rstar_val    = 0.2    # world rate shock persistence
ystar_ss_val    = 1.0    # steady-state foreign output
PVstar_ss       = 1.0    # steady-state foreign price (normalized)
sigmaH_val      = 0.999  # Armington elasticity (home vs. foreign)

# Export and home-bias
modchiX   = let   # sectoral export shares chi_i^X from Chilean 2021 supply-use table (Data/computed)
    _f = joinpath(DATA_DIR, "computed", "export_shares_chile.csv")
    isfile(_f) ? (v = Float64.(CSV.read(_f, DataFrame).chi_x); v ./ sum(v)) : fill(1/nsec, nsec)
end
modvarrho = var_rho              # home-bias parameters ϱ_i from data
modA      = ones(nsec)           # TFP normalized to 1 in SS

# Default shock parameters (overridden by exercise below)
rho_val          = 0.1
sigma_L_agg_val  = 0.0
sigma_L_het      = zeros(nsec)

# Import price (PVstar) shock
rho_pvstar_val   = 0.9
sigma_pvstar_val = 0.03

# Oil price (POstar) shock
# alphaOilShare: fraction of each sector's composite imports V_i that is oil/refined fuels.
# Source: Cuadro 21 (intermediate imported use at basic prices), MIP Chile 2021,
#         products 28,77-81 (crude oil + diesel, gasoline, kerosene, fuel oils, LPG).
#         Verified against 2021_Cuadros_12x12.xlsx.
# Sector order: 1=Agro/Fishing, 2=Mining, 3=Manufacturing, 4=Utilities,
#               5=Construction, 6=Trade/Hotels, 7=Transport/Comms, 8=Finance,
#               9=Real Estate, 10=Business Svcs, 11=Personal Svcs, 12=Public Admin
modalphaOil = [0.1635,   # 1  Agropecuario-silvícola y Pesca
               0.2164,   # 2  Minería
               0.1890,   # 3  Industria manufacturera
               0.0871,   # 4  Electricidad, gas, agua y gestión de desechos
               0.0417,   # 5  Construcción
               0.0793,   # 6  Comercio, hoteles y restaurantes
               0.3734,   # 7  Transporte, comunicaciones y servicios de información
               0.0043,   # 8  Intermediación financiera
               0.0449,   # 9  Servicios inmobiliarios y de vivienda
               0.0633,   # 10 Servicios empresariales
               0.0391,   # 11 Servicios personales
               0.0447]   # 12 Administración pública
epsilonV_oil_val  = 0.5           # CES elasticity between oil and non-oil imports
rho_postar_val    = 0.9           # AR(1) persistence of world oil price
sigma_postar_val  = 0.02          # std dev of oil price innovation
POstar_ss_val     = 1.0           # SS world oil price (normalized, same as PVstar_ss)

# Preference (xi) shock
rho_xi_val   = 0.80
sigma_xi_val = 0.005


# =========================================================================== #
#  SHOCK PARAMETERS BY EXERCISE                                                #
# =========================================================================== #

if EXERCISE == 0
    # --- Baseline: all shocks active (Option-A: 12 sectoral demand shocks) ---
    sigma_i_val    = 0.001
    rho_om1_val    = 0.1;   rho_tfp2_val = 0.0
    sigma_om_vec   = fill(0.03, nsec)   # uniform initial amplitude
    rho_tfp1_val   = 0.5;   rho_tfp2_val = 0.0
    isigma_tfp_val = fill(0.05, nsec)

elseif EXERCISE == 1
    # --- Exercise 1: sectoral demand shocks only ---
    sigma_i_val    = 0.0
    rho_om1_val    = 0.95
    sigma_om_vec   = fill(0.01, nsec)
    rho_tfp1_val   = 0.5;   rho_tfp2_val = 0.0
    isigma_tfp_val = zeros(nsec)

elseif EXERCISE == 2
    # --- Exercise 2: manufacturing TFP shock only ---
    mfg = 3   # sector 3 = Manufacturing
    sigma_i_val    = 0.0
    rho_om1_val    = 0.1
    sigma_om_vec   = zeros(nsec)
    rho_tfp1_val   = 0.95;  rho_tfp2_val = 0.0
    isigma_tfp_val = zeros(nsec)
    isigma_tfp_val[mfg] = 0.01

elseif EXERCISE == 3
    # --- Exercise 3: monetary policy shock only ---
    sigma_i_val    = 0.01
    rho_om1_val    = 0.1
    sigma_om_vec   = zeros(nsec)
    rho_tfp1_val   = 0.5;   rho_tfp2_val = 0.0
    isigma_tfp_val = zeros(nsec)
end

# Shock variance indicators for Dynare shocks block
shock_eps_om_vec     = Float64.(sigma_om_vec .> 0)
shock_eps_i_val      = Float64(sigma_i_val  > 0)
shock_eps_pvstar_val = Float64(EXERCISE == 0 && sigma_pvstar_val > 0)
shock_eps_xi_val     = Float64(EXERCISE == 0 && sigma_xi_val > 0)
shock_eps_postar_val = Float64(EXERCISE == 0 && sigma_postar_val > 0 && any(modalphaOil .> 0))
shock_epsA_val       = ones(nsec)


# =========================================================================== #
#  LOAD SMM ESTIMATES (override defaults when available)                       #
# =========================================================================== #

# SMM estimates: smm_estimates.csv is written by a COMPLETED estimation;
# smm_checkpoint.csv holds the live best of a running/interrupted one.
# Use whichever is NEWER (2026-07-10) — no more manual cp needed.
smm_est_file  = joinpath(ESTIMATION_DIR, "smm_estimates.csv")
smm_ckpt_file = joinpath(ESTIMATION_DIR, "smm_checkpoint.csv")

param_names = ["ilabcosts", "epsY", "epsM", "kappaV", "rho_om1",
               "rho_tfp1", "isigma_tfp(1)", "sigma_om(avg)", "rho_pvstar", "sigma_pvstar"]

smm_param_source = "hard-coded defaults"

_theta_candidates = [(f, mtime(f)) for f in (smm_est_file, smm_ckpt_file) if isfile(f)]
sort!(_theta_candidates, by = x -> x[2], rev = true)   # newest first
_theta_file = isempty(_theta_candidates) ? "" : _theta_candidates[1][1]

if !isempty(_theta_file)
    est_df = CSV.read(_theta_file, DataFrame)
    if nrow(est_df) != length(CSV_PARAM_NAMES)
        @printf "  WARNING: %s has %d params (current layout: %d) — stale estimation layout.\n" basename(_theta_file) nrow(est_df) length(CSV_PARAM_NAMES)
    end
    _obj_prev = try Float64(est_df.obj_hat[1]) catch; NaN end
    @printf "  θ source: %s  (obj=%.4f, written %s)%s\n" basename(_theta_file) _obj_prev Libc.strftime("%Y-%m-%d %H:%M", mtime(_theta_file)) (endswith(_theta_file, "checkpoint.csv") ? "  [live best of an unfinished run]" : "")
    est    = Dict(String(r.param) => Float64(r.value) for r in eachrow(est_df))

    ilabcosts_val    = est["ilabcosts"]
    modepsY          = fill(est["epsY"], nsec)
    modepsM          = fill(est["epsM"], nsec)
    kappaV_val       = exp(est["log_kappaV"])
    rho_om1_val      = est["rho_om"]
    rho_tfp1_val     = est["rho_A"]
    isigma_tfp_val   = [est["isigma_tfp_$(i)"] for i in 1:nsec]
    # Option-A: 12 sectoral demand shock amplitudes
    sigma_om_vec     = [haskey(est, "sigma_om_$(i)") ? est["sigma_om_$(i)"] : 0.03 for i in 1:nsec]
    rho_pvstar_val   = est["rho_pvstar"]
    sigma_pvstar_val = est["sigma_pvstar"]
    haskey(est, "rho_xi")   && (rho_xi_val   = est["rho_xi"])
    haskey(est, "sigma_xi") && (sigma_xi_val = est["sigma_xi"])
    haskey(est, "etastar")  && (etastar_val  = est["etastar"])
    haskey(est, "kappaw")   && (kappaw_val   = est["kappaw"])   # θ[36], added 2026-07-08
    shock_eps_om_vec = Float64.(sigma_om_vec .> 0)   # recompute flags from loaded values
    smm_param_source = basename(_theta_file)
end

# Enforce the AR-persistence box [0, 0.99] on the parameters that will drive the
# simulation, regardless of source. A stored estimate from BEFORE the [0,1)
# non-negativity floor (e.g. rho_om = -0.82) would otherwise produce oscillatory
# (period-2, zig-zag) IRFs. Clamp and warn — for a proper fit, re-estimate the
# SMM under the new bounds rather than relying on this projection.
for (nm, v) in (("rho_om", rho_om1_val), ("rho_A", rho_tfp1_val),
                ("rho_pvstar", rho_pvstar_val), ("rho_xi", rho_xi_val))
    (v < 0 || v >= 1) && @printf "  WARNING: loaded %s = %.4f outside [0,1) — clamped to [0,0.99]; re-estimate for a proper fit\n" nm v
end
rho_om1_val    = clamp(rho_om1_val,    0.0, 0.99)
rho_tfp1_val   = clamp(rho_tfp1_val,   0.0, 0.99)
rho_pvstar_val = clamp(rho_pvstar_val, 0.0, 0.99)
rho_xi_val     = clamp(rho_xi_val,     0.0, 0.99)

# --- Redesign: fix production elasticities at Atalay (2017) complementarity
# values (epsY ~ 0.80 = eps_Q, epsM ~ 0.10 = eps_m), overriding any stored
# estimate. Low elasticities keep sector-specific shocks from washing out
# through the network (Atalay: eps_m=0.1 -> 83% sector-specific vs 21% at eps=1).
_epsY = parse(Float64, get(ENV, "SMM_EPSY", "0.80"))
_epsM = parse(Float64, get(ENV, "SMM_EPSM", "0.20"))
modepsY = fill(_epsY, nsec)
modepsM = fill(_epsM, nsec)
@printf "  [redesign] elasticities fixed: epsY=%.2f epsM=%.2f (Atalay eps_Q~0.84 eps_m~0.10; env SMM_EPSY/SMM_EPSM)\n" _epsY _epsM

# --- Redesign: override sectoral + external shock SIZES with the DATA-MEASURED
# calibrations (compute_sectoral_shocks.jl / compute_external_shock.jl), so the
# shocks are data-disciplined instead of SMM-inflated. No-op if the CSVs are absent.
if isfile(joinpath(DATA_DIR, "sectoral_shock_calibration.csv"))
    _sc = CSV.read(joinpath(DATA_DIR, "sectoral_shock_calibration.csv"), DataFrame)
    if nrow(_sc) == nsec
        isigma_tfp_val   = Float64.(_sc.isigma_tfp_init)
        sigma_om_vec     = Float64.(_sc.sigma_om_init)
        shock_eps_om_vec = Float64.(sigma_om_vec .> 0)
        @printf "  [redesign] sectoral shock sizes <- sectoral_shock_calibration.csv\n"
    else
        @printf "  [redesign] WARNING sectoral_shock_calibration.csv has %d rows != nsec=%d; ignored\n" nrow(_sc) nsec
    end
end
if isfile(joinpath(DATA_DIR, "external_shock_calibration.csv"))
    _ec = CSV.read(joinpath(DATA_DIR, "external_shock_calibration.csv"), DataFrame)
    _em = Dict(String.(_ec.param) .=> Float64.(_ec.value))
    rho_pvstar_val   = get(_em, "rho_pvstar",   rho_pvstar_val)
    sigma_pvstar_val = get(_em, "sigma_pvstar", sigma_pvstar_val)
    @printf "  [redesign] external shock <- external_shock_calibration.csv (rho=%.3f sigma=%.3f)\n" rho_pvstar_val sigma_pvstar_val
end

@printf "--- Parameters (%s) ---\n" smm_param_source
param_vals = [ilabcosts_val, modepsY[1], modepsM[1], kappaV_val,
              rho_om1_val, rho_tfp1_val, isigma_tfp_val[1], mean(sigma_om_vec),
              rho_pvstar_val, sigma_pvstar_val,
              rho_xi_val, sigma_xi_val, etastar_val, kappaw_val]
for (nm, vl) in zip(vcat(param_names, ["rho_xi", "sigma_xi", "etastar", "kappaw"]), param_vals)
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

# Load trade-balance target from aggregate_moments.csv if available
# Use assignment-expression form to avoid soft-scope ambiguity (Julia 1.12)
tb_target = let
    _tb = 0.02
    agg_mom_tb = joinpath(DATA_DIR, "aggregate_moments.csv")
    if isfile(agg_mom_tb)
        try
            agg_tb = CSV.read(agg_mom_tb, DataFrame)
            row    = filter(r -> String(r.moment) == "TBGDP", agg_tb)
            !isempty(row) && isfinite(row[1, :value]) && (_tb = Float64(row[1, :value]))
        catch
        end
    end
    _tb
end

ss_result = nlsolve(
    (F, x) -> F .= steady_ntwsoe(
        x, PVstar_ss, epsilon, varrho_val, sigmaH,
        gammag_vec, gammas_vec, om_g, om_s, chiX_vec, omegaX, etastar, Ystar,
        alpha_vec, alphaV_vec, beta_mat, epsY_vec, epsM_vec,
        gamma, chi, psi, A_vec, tb_target
    ),
    x_guess;
    ftol = 1e-10,
    show_trace = get(ENV, "SS_TRACE", "0") == "1",   # 1000-iter trace off by default (SS_TRACE=1 to debug)
    method = :trust_region,
)

if !converged(ss_result)
    @warn "Outer steady-state solver did not converge (residual_norm=$(ss_result.residual_norm)). Proceeding anyway — results may be inaccurate."
end

pH_ss = ss_result.zero[1:nsec]
w_ss  = ss_result.zero[nsec+1]
Q_ss  = ss_result.zero[nsec+2]
C_ss  = ss_result.zero[nsec+3]

# GHH preferences: disutility scale chi0 = chi_weight * C_ss^gamma. This holds
# the steady state EXACTLY at the KPR-solved values (the SS labor-supply residual
# in steady_ntwsoe.jl is left in KPR form N = (C^-gamma w / chi)^(1/psi), so the
# SS allocation is unchanged; chi0 only enters the DYNAMIC MRS = chi0*chi*N^psi
# and the composite Zc = C - chi0*chi*N^(1+psi)/(1+psi)). Recomputed from C_ss so
# it tracks the SS if C_ss moves during estimation.
chi0_val = chi * C_ss^gamma


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
Y_tot_ss      = sum(pH_ss .* Yi_ss)   # gross output at constant SS relative prices (matches Y in .mod)
VA_ss_val     = sum(pH_ss .* Yi_ss .- PMi_ss .* M_ss .- PV_ss .* Vi_ss)   # double-deflated VA at SS prices (matches VA in .mod; PIV_ss_i = PV_ss at SS)
Ctotg_ss_val  = sum(gammag_vec .* (p_g_ss ./ P_ss) .* C_g_ss)
Ctots_ss_val  = sum(gammas_vec .* (p_s_ss ./ P_ss) .* C_s_ss)
Ctot_ss_val   = Ctotg_ss_val + Ctots_ss_val
IMP_ss_val    = IMP_tot_ss

# Aggregates used in initval block (must be declared as Dynare parameters)
CFg_total_ss  = sum(CFg_ss)   # total goods consumer imports
CFs_total_ss  = sum(CFs_ss)   # total services consumer imports
# PX_ss, V_ss, CF_ss already computed above

Pistar_ss_val  = Pistar_ss
Rworld_ss_val  = Rworld_ss
PVstar_ss_val  = PVstar_ss

# Oil sector SS: PO_ss = Q_ss * POstar_ss; with POstar_ss = PVstar_ss = 1 -> PO_ss = PV_ss
PO_ss     = Q_ss * POstar_ss_val
# Sector-specific composite import price index at SS:
# PIV_i^(1-eps) = alphaOil_i * PO^(1-eps) + (1-alphaOil_i) * PV^(1-eps)
# When POstar_ss = PVstar_ss = 1 -> PO_ss = PV_ss -> PIV_i_ss = PV_ss for all i
PIV_ss_vec = [
    (modalphaOil[i] * PO_ss^(1 - epsilonV_oil_val)
     + (1 - modalphaOil[i]) * PV_ss^(1 - epsilonV_oil_val))^(1 / (1 - epsilonV_oil_val))
    for i in 1:nsec
]

# Foreign demand shock params (kept for smm_estimates.mat compatibility)
rho_psi_val   = 0.5
sigma_psi_val = 0.001

@printf "  GDP=%.4f  TB/GDP=%.4f  X=%.4f  IMP=%.4f  bbar=%.4f\n" GDP_ss tbgdp X_ss IMP_tot_ss bbar_val
@printf "  TB/GDP target: %.4f  |  world rate (ann.): %.2f%%\n" tb_target 400*(Rworld_ss-1)
@printf "  Bstar=%.4f  Q*Bstar=%.4f  Q*Bstar/GDP=%.3f (%.1f%%)\n" Bstar_ss (Q_ss*Bstar_ss) (Q_ss*Bstar_ss/GDP_ss) (100*Q_ss*Bstar_ss/GDP_ss)
@printf "  Gross output/GDP : %.3f  (IMF target: %.1f)%s\n" ygdp_ss ygdp_target (ternary_str(abs(ygdp_ss-ygdp_target)>0.20, "  << off target", ""))
@printf "  Foreign debt/GDP : %.3f  (IMF target: %.2f)%s\n\n" fdebt_ss fdebt_target (ternary_str(abs(fdebt_ss-fdebt_target)>0.10, "  << off target", ""))


# =========================================================================== #
#  SECTORAL GDP / VALUE-ADDED SHARES                                          #
#  VA_net = Yi - M - Vi  (net of domestic materials AND imported inputs)      #
#  VA_mdl = Yi - M       (model's VA_ss_val convention: domestic mats only)   #
#  %GDP shares normalize by GDP_ss = C_ss + TB_ss (expenditure side).         #
# =========================================================================== #

VAnet_vec = Yi_ss .- M_ss .- Vi_ss
VAmdl_vec = Yi_ss .- M_ss
VAnet_tot = sum(VAnet_vec)
VAmdl_tot = sum(VAmdl_vec)

# Fixed-width sector name: %-Ns pads but does NOT truncate, so long names
# ("Electricidad gas agua y gestion de desechos") destroyed column alignment.
_secname(s, w=32) = length(s) <= w ? s : first(s, w - 1) * "…"

@printf "--- Sectoral GDP / value-added shares ---\n"
@printf "  GDP_ss = %.4f   ΣVA(Y-M-V) = %.4f   ΣVA(Y-M) = %.4f\n\n" GDP_ss VAnet_tot VAmdl_tot
@printf "%-32s %11s %8s %8s %11s %8s %8s\n" "Sector" "VA(Y-M-V)" "%GDP" "%ΣVA" "VA(Y-M)" "%GDP" "%ΣVA"
@printf "%s\n" repeat("-", 92)
for i in 1:nsec
    @printf("%-32s %11.4f %7.2f%% %7.2f%% %11.4f %7.2f%% %7.2f%%\n",
            _secname(names_vec[i]),
            VAnet_vec[i], 100*VAnet_vec[i]/GDP_ss, 100*VAnet_vec[i]/VAnet_tot,
            VAmdl_vec[i], 100*VAmdl_vec[i]/GDP_ss, 100*VAmdl_vec[i]/VAmdl_tot)
end
@printf "%s\n" repeat("-", 92)
@printf("%-32s %11.4f %7.2f%% %7.2f%% %11.4f %7.2f%% %7.2f%%\n",
        "TOTAL",
        VAnet_tot, 100*VAnet_tot/GDP_ss, 100.0,
        VAmdl_tot, 100*VAmdl_tot/GDP_ss, 100.0)
@printf "\n"

# Save to tables/gdp_va_shares_<tag>.csv
df_va_shares = DataFrame(
    sector        = 1:nsec,
    name          = names_vec,
    VA_net        = VAnet_vec,                 # Yi - M - Vi
    VA_mdl        = VAmdl_vec,                 # Yi - M
    share_GDP_net = VAnet_vec ./ GDP_ss,
    share_VA_net  = VAnet_vec ./ VAnet_tot,
    share_GDP_mdl = VAmdl_vec ./ GDP_ss,
    share_VA_mdl  = VAmdl_vec ./ VAmdl_tot,
)
CSV.write(joinpath(TABLES_DIR, "gdp_va_shares_$(tag).csv"), df_va_shares)
@printf "  → saved to: %s\n\n" joinpath(TABLES_DIR, "gdp_va_shares_$(tag).csv")


# =========================================================================== #
#  EXPANDED VALUE ADDED  (VAE — "valor agregado expandido")                   #
#  Method: Foster & Valdés (2012); BCCh Estudio Económico Estadístico N°148   #
#  (Chovar & Leiva, 2026).                                                    #
#                                                                             #
#  Idea: a sector's relevance is not only its own value added but also the    #
#  slices of OTHER sectors' value added that exist because of its position in #
#  the production network — both as an input SUPPLIER (forward links) and as  #
#  an input BUYER (backward links). Intra-sector links are excluded (they are #
#  a zero-sum transfer already counted in the sector's own VA).               #
#                                                                             #
#  All quantities are NOMINAL (price × steady-state quantity).                #
#    pY_j        = pH_j · Y_j                       (gross output, "TVT")      #
#    VA_j        = pY_j − PMi_j·M_j − PV·Vi_j       (traditional value added)  #
#    Z[i,j]      = nominal intermediates sector i buys from sector j          #
#                = Γ_{i,j}·PMi_i^{εm}·pH_j^{1-εm}·M_i   (CES cost-min, ΣⱼZ=PMi·M)#
#    totinput_j  = PMi_j·M_j + PV·Vi_j              (all intermediates: nat+imp)#
#                                                                             #
#  Forward (eq. 5): sector j buys S's output ⇒ slice of VA_j credited to S    #
#    ΔVAᶠ_j(S) = [Z[j,S]/totinput_j] · VA_j                                    #
#  Backward (eq. 6): S buys from j ⇒ slice of VA_j credited to S              #
#    ΔVAᵇ_j(S) = [Z[S,j]/pY_j]      · VA_j                                     #
#  Expanded: VAE_S = VA_S + Σ_{j≠S} ΔVAᶠ_j(S) + Σ_{j≠S} ΔVAᵇ_j(S)             #
#                                                                             #
#  NOTE / open issue (flag for Agustín): in this model ALL bilateral inter-   #
#  mediate flows are domestic (modbeta is the domestic IO matrix); imported   #
#  inputs Vi are an aggregate composite, not a bilateral flow. Hence the      #
#  national/total ratios Xᴺ/Xᵀ in eqs (5)-(6) collapse to 1, so model VAE is  #
#  an upper bound vs the data version that discounts imported-input content.  #
# =========================================================================== #

pY_ss_nom   = pH_ss .* Yi_ss                       # nominal gross output (TVT_j)
matcost_nom = PMi_ss .* M_ss                       # domestic materials cost
impcost_nom = PV_ss .* Vi_ss                       # imported-input cost
VA_nom      = pY_ss_nom .- matcost_nom .- impcost_nom   # traditional nominal VA
totinput    = matcost_nom .+ impcost_nom           # total intermediate cost (nat+imp)
VA_nom_tot  = sum(VA_nom)

# Bilateral nominal domestic intermediate-flow matrix: Z[i,j] = i buys from j
Zflow = [ modbeta[i, j] * PMi_ss[i]^epsM_vec[i] * pH_ss[j]^(1 - epsM_vec[i]) * M_ss[i]
          for i in 1:nsec, j in 1:nsec ]

# Import discounts (proxy for the paper's national/total ratio Xᴺ/Xᵀ, which the
# model lacks at the bilateral level). matcost/totinput = domestic share of a
# sector's intermediate bill (buyer side); (1-α_V) = domestic share (supplier side).
dom_input_share  = [totinput[j] > 0 ? matcost_nom[j] / totinput[j] : 1.0 for j in 1:nsec]
dom_supply_share = 1 .- modalphaV

VAE_fwd = zeros(nsec)   # forward-link contributions received by each sector S
VAE_bwd = zeros(nsec)   # backward-link contributions received by each sector S
VAE_fwd_adj = zeros(nsec)   # import-adjusted versions
VAE_bwd_adj = zeros(nsec)
for S in 1:nsec, j in 1:nsec
    j == S && continue
    fwd = (totinput[j]   > 0 ? Zflow[j, S] / totinput[j]   : 0.0) * VA_nom[j]
    bwd = (pY_ss_nom[j]  > 0 ? Zflow[S, j] / pY_ss_nom[j]  : 0.0) * VA_nom[j]
    VAE_fwd[S] += fwd
    VAE_bwd[S] += bwd
    VAE_fwd_adj[S] += fwd * dom_input_share[j]   # discount by buyer j's domestic-input share
    VAE_bwd_adj[S] += bwd * dom_supply_share[j]  # discount by supplier j's domestic share (1-α_V)
end
VAE_vec     = VA_nom .+ VAE_fwd     .+ VAE_bwd
VAE_vec_adj = VA_nom .+ VAE_fwd_adj .+ VAE_bwd_adj

@printf "--- Expanded value added (VAE — Foster & Valdés 2012 / EEE-148) ---\n"
@printf "  Shares are %% of total nominal value added (ΣVA = %.4f).\n" VA_nom_tot
@printf "  VAE shares do NOT sum to 100%% — VAE is a network-influence measure, not a partition.\n\n"
@printf "%-32s %8s %8s %8s %8s %9s\n" "Sector" "VA(dir)" "Fwd" "Bwd" "VAE" "VAE/VA"
@printf "%s\n" repeat("-", 80)
for i in 1:nsec
    @printf("%-32s %7.2f%% %7.2f%% %7.2f%% %7.2f%% %8.2fx\n",
            _secname(names_vec[i]),
            100*VA_nom[i]/VA_nom_tot,
            100*VAE_fwd[i]/VA_nom_tot,
            100*VAE_bwd[i]/VA_nom_tot,
            100*VAE_vec[i]/VA_nom_tot,
            VA_nom[i] > 0 ? VAE_vec[i]/VA_nom[i] : NaN)
end
@printf "%s\n" repeat("-", 80)
@printf("%-32s %7.2f%% %7.2f%% %7.2f%% %7.2f%%\n",
        "TOTAL (direct sums to 100%)",
        100*sum(VA_nom)/VA_nom_tot, 100*sum(VAE_fwd)/VA_nom_tot,
        100*sum(VAE_bwd)/VA_nom_tot, 100*sum(VAE_vec)/VA_nom_tot)
@printf "\n"

# Import-adjusted VAE (discounts imported-input content; comparable to EEE-148)
@printf "--- Import-adjusted VAE (fwd × buyer domestic-input share; bwd × supplier (1-α_V)) ---\n"
@printf "%-32s %8s %8s %8s %8s %9s\n" "Sector" "VA(dir)" "Fwd*" "Bwd*" "VAE*" "VAE*/VA"
@printf "%s\n" repeat("-", 80)
for i in 1:nsec
    @printf("%-32s %7.2f%% %7.2f%% %7.2f%% %7.2f%% %8.2fx\n",
            _secname(names_vec[i]),
            100*VA_nom[i]/VA_nom_tot,
            100*VAE_fwd_adj[i]/VA_nom_tot,
            100*VAE_bwd_adj[i]/VA_nom_tot,
            100*VAE_vec_adj[i]/VA_nom_tot,
            VA_nom[i] > 0 ? VAE_vec_adj[i]/VA_nom[i] : NaN)
end
@printf "%s\n" repeat("-", 80)
@printf("%-32s %7.2f%% %7.2f%% %7.2f%% %7.2f%%\n",
        "TOTAL",
        100*sum(VA_nom)/VA_nom_tot, 100*sum(VAE_fwd_adj)/VA_nom_tot,
        100*sum(VAE_bwd_adj)/VA_nom_tot, 100*sum(VAE_vec_adj)/VA_nom_tot)
@printf "\n"

df_vae = DataFrame(
    sector            = 1:nsec,
    name              = names_vec,
    VA_direct         = VA_nom,                        # own nominal value added
    fwd_links         = VAE_fwd,                       # forward-link VA received (raw)
    bwd_links         = VAE_bwd,                        # backward-link VA received (raw)
    VAE               = VAE_vec,                        # expanded value added (raw)
    fwd_links_adj     = VAE_fwd_adj,                    # forward, import-adjusted
    bwd_links_adj     = VAE_bwd_adj,                    # backward, import-adjusted
    VAE_adj           = VAE_vec_adj,                    # expanded VA, import-adjusted
    share_VA_direct   = VA_nom      ./ VA_nom_tot,      # direct share of total VA
    share_VAE         = VAE_vec     ./ VA_nom_tot,      # expanded share (raw)
    share_VAE_adj     = VAE_vec_adj ./ VA_nom_tot,      # expanded share (import-adjusted)
    vae_multiplier    = VAE_vec     ./ VA_nom,          # raw VAE / direct VA
    vae_multiplier_adj= VAE_vec_adj ./ VA_nom,          # adjusted VAE / direct VA
)
CSV.write(joinpath(TABLES_DIR, "vae_$(tag).csv"), df_vae)
@printf "  → saved to: %s\n\n" joinpath(TABLES_DIR, "vae_$(tag).csv")


# =========================================================================== #
#  WRITE params_jl.mod  (read by NK_SOE_lev_gap2.mod via @#include)          #
#  Dynare.jl processes this as a Dynare parameter file — no MATLAB/Octave.   #
# =========================================================================== #

params_nt = (
    nsec           = nsec,
    # Shock scalars
    sigma_i_val    = sigma_i_val,
    sigma_L_agg_val= sigma_L_agg_val,
    ilabcosts_val  = ilabcosts_val,
    gamma_val      = gamma,
    psi_val        = psi,       # inverse Frisch (env SMM_PSI); written to params_jl.mod
    chi0_val       = chi0_val,   # GHH disutility scale = chi * C_ss^gamma
    beta_val       = beta_val,
    phi_val        = phi_val,
    rho_val        = rho_val,
    rho_om1_val    = rho_om1_val,
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
    kappaw_val     = kappaw_val,   # write_params_mod defaults to 115.0 if absent
    # Oil sector
    modalphaOil         = modalphaOil,
    epsilonV_oil_val    = epsilonV_oil_val,
    rho_postar_val      = rho_postar_val,
    sigma_postar_val    = sigma_postar_val,
    POstar_ss_val       = POstar_ss_val,
    shock_eps_postar_val = shock_eps_postar_val,
    PIV_ss_vec          = PIV_ss_vec,
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
    Q_ss           = Q_ss,
    TB_ss          = TB_ss,
    PX_ss          = PX_ss,
    V_ss           = V_ss,
    CF_ss          = CF_ss,
    CFg_total_ss   = CFg_total_ss,
    CFs_total_ss   = CFs_total_ss,
    IMP_ss_val     = IMP_ss_val,
    Ctot_ss_val    = Ctot_ss_val,
    Ctotg_ss_val   = Ctotg_ss_val,
    Ctots_ss_val   = Ctots_ss_val,
    VA_ss_val      = VA_ss_val,
    M_tot_ss       = M_tot_ss,
    Y_tot_ss       = Y_tot_ss,
    # Shock flags
    shock_eps_i_val      = shock_eps_i_val,
    shock_eps_pvstar_val = shock_eps_pvstar_val,
    shock_eps_xi_val     = shock_eps_xi_val,
    shock_epsA_val       = shock_epsA_val,
    # Option-A: 12 sectoral demand shock parameters
    sigma_om_vec         = sigma_om_vec,
    shock_eps_om_vec     = shock_eps_om_vec,
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
# @dynare macro requires a string LITERAL (not a variable) and the
# Dynare C++ preprocessor runs with the CWD at the time @dynare expands.
# Strategy: include() a small helper file from inside this function.
# include() called at runtime (inside a function) re-compiles and runs
# the included file at call time — so the preprocessor sees CWD = MOD_DIR.
# =========================================================================== #
#  RUN DYNARE IN A SUBPROCESS                                                  #
#  @dynare inside a function triggers Julia 1.12 world-age crashes.           #
#  Solution: call a standalone script in a fresh Julia process.               #
#  The subprocess writes CSV files; we read them back here.                   #
# =========================================================================== #

dynare_script = joinpath(SCRIPT_DIR, "run_dynare_subprocess.jl")
julia_exe     = joinpath(Sys.BINDIR, "julia")

# Use the PARENT process's active project (where Dynare is already installed).
# The julia_dynare/Project.toml lists Dynare but requires `Pkg.instantiate()`
# before the subprocess can use it — the parent's env already has it.
project_dir = dirname(Base.active_project())

@printf "--- Running Dynare.jl (subprocess) ---\n"
@printf "  Script : %s\n" dynare_script
@printf "  Mod dir: %s\n\n" MOD_DIR

# Capture stdout to check for DYNARE_SUCCESS sentinel.
# ignorestatus: the subprocess exits 1 on a failed solve (DYNARE_FAILED);
# we want to report that cleanly instead of a raw ProcessFailedException.
dynare_out = IOBuffer()
proc = run(pipeline(
    ignorestatus(`$julia_exe --project=$project_dir $dynare_script $MOD_DIR`),
    stdout=dynare_out, stderr=stderr), wait=true)

dynare_stdout = String(take!(dynare_out))
print(dynare_stdout)   # echo subprocess output

if occursin("DYNARE_FAILED", dynare_stdout)
    error("""
    Dynare solved the steady state but the FIRST-ORDER (LRE) SOLVE FAILED:
    decision rule is all zeros (see DECISION_RULE_NONZERO=0 above, and any
    MethodError printed by Dynare.jl during parse_statements!).
    No context was saved — the SMM stage will not run.
    Likely causes: parameterization outside determinacy region, or a
    Dynare.jl dependency version mismatch (check LinearRationalExpectations /
    FastLapackInterface / PolynomialMatrixEquations against the working env).
    """)
elseif !occursin("DYNARE_SUCCESS", dynare_stdout) || !success(proc)
    error("""
    Dynare subprocess did not complete successfully (exit code $(proc.exitcode)).
    Check the output above for errors.
    Common fixes:
      - Open pib_sectorial_bc.xlsx in Excel and File→Save As .xlsx
      - Make sure params_jl.mod was written to: $MOD_DIR
    """)
end
@printf "\n--- Dynare subprocess completed ---\n\n"

# =========================================================================== #
#  READ DYNARE RESULTS FROM CSV                                                #
# =========================================================================== #

@printf "--- Loading Dynare results ---\n"

df_names   = CSV.read(joinpath(MOD_DIR, "dynare_endo_names.csv"),  DataFrame)
df_ss      = CSV.read(joinpath(MOD_DIR, "dynare_ss.csv"),          DataFrame)
df_g1_1    = CSV.read(joinpath(MOD_DIR, "dynare_g1_1.csv"),        DataFrame)
df_g1_2    = CSV.read(joinpath(MOD_DIR, "dynare_g1_2.csv"),        DataFrame)
df_sigma_e = CSV.read(joinpath(MOD_DIR, "dynare_sigma_e.csv"),     DataFrame)
df_states  = CSV.read(joinpath(MOD_DIR, "dynare_state_rows.csv"),  DataFrame)

endo_names = String.(df_names.variable)
ss_vec     = Float64.(df_ss.ss_value)
ghx_jl     = Matrix{Float64}(df_g1_1)    # n_endo × n_states
ghu_jl     = Matrix{Float64}(df_g1_2)    # n_endo × n_shocks
Σe_jl      = Matrix{Float64}(df_sigma_e)
state_rows = Int.(df_states.state_row)

# Simulated paths (periods × variables)
sim_file   = joinpath(MOD_DIR, "dynare_sim.csv")
sim_matrix = isfile(sim_file) ? Matrix{Float64}(CSV.read(sim_file, DataFrame)) :
                                 fill(NaN, 0, length(endo_names))

ss_ok  = !any(isnan, ss_vec) && !any(isinf, ss_vec)
sim_ok = !isempty(sim_matrix)

@printf "  Steady state : %s  (%d variables)\n" (ss_ok ? "OK" : "FAILED") length(ss_vec)
@printf "  Simulation   : %s  (%d periods)\n\n" (sim_ok ? "OK" : "no data") size(sim_matrix,1)

# Build variable-name → index map
endo_idx = Dict(nm => i for (i, nm) in enumerate(endo_names))


# =========================================================================== #
#  HP-FILTERED MOMENTS  (replicates smm_estimation.jl / smm_model_moments.jl)#
# =========================================================================== #

@printf "--- Computing HP-filtered moments (matching SMM estimation) ---\n"

std_Y_m  = fill(NaN, nsec)
std_PH_m = fill(NaN, nsec)
std_L_m  = fill(NaN, nsec)

rc_lyap_ok = false
m_std_GDP    = NaN
m_std_pi     = NaN
m_std_Q      = NaN
m_corr_GDPpi = NaN
m_corr_GDPQ  = NaN
m_omG        = ombar_val   # fixed by calibration
m_autocorr_Q = NaN
m_std_TBGDP  = NaN
rho_y = 0.0; rho_p = 0.0; rho_l = 0.0
corr_YPH_m   = fill(NaN, nsec)   # corr(Y_i,PH_i): supply vs demand identifier
m_corr_NGDP  = NaN               # corr(N,GDP)
m_corr_NAPL  = NaN               # corr(N,GDP/N)

try
    n_exo   = size(ghu_jl, 2)
    n_state = length(state_rows)

    # Shock covariance: activate BY NAME, matching active_shock_indices() in
    # smm_model_moments.jl exactly (eps_i, eps_pvstar, eps_xi, epsA_1:12,
    # eps_om_1:12 — epschi and eps_postar stay off). The previous hardcoded
    # index map was from the old 18-shock layout: it missed eps_pvstar and all
    # eps_om_i and wrongly activated epschi/eps_postar (fixed 2026-07-10).
    # Shock std devs are in-equation parameters, so unit innovation variances
    # are correct here.
    exo_names_file = joinpath(MOD_DIR, "dynare_exo_names.csv")
    Σe_smm = zeros(n_exo, n_exo)
    n_active = 0
    if isfile(exo_names_file)
        exo_names_jl = String.(CSV.read(exo_names_file, DataFrame).shock)
        active_set = Set(vcat(["eps_i", "eps_pvstar", "eps_xi"],
                              ["epsA_$(i)"   for i in 1:nsec],
                              ["eps_om_$(i)" for i in 1:nsec]))
        for (k, nm) in enumerate(exo_names_jl)
            k > n_exo && break
            if nm in active_set
                Σe_smm[k,k] = 1.0; n_active += 1
            end
        end
        @printf "  Active shocks (by name): %d of %d exogenous\n" n_active n_exo
    end
    # GUARD (2026-07-10): an empty/missing/unmatched name file must NOT
    # silently produce a zero-shock economy (all model moments = 0, useless
    # fit table). Fall back to all-shocks-active with a loud warning —
    # slightly overstates volatility vs the SMM's 27-shock set (adds epschi,
    # eps_postar) but yields an informative table.
    if n_active == 0
        for k in 1:n_exo; Σe_smm[k,k] = 1.0; end
        @printf "  WARNING: no shocks matched by name (%s) — ALL %d shocks active; rerun the Dynare subprocess for the exact SMM shock set\n" (isfile(exo_names_file) ? "empty/unmatched dynare_exo_names.csv" : "file missing") n_exo
    end

    # State-space matrices
    Tsr = ghx_jl[state_rows, :]   # n_state × n_state
    Rsr = ghu_jl[state_rows, :]   # n_state × n_exo
    B_Σ_Bt = Rsr * Σe_smm * Rsr'
    B_Σ_Bt = (B_Σ_Bt + B_Σ_Bt') / 2

    # Solve state covariance via Lyapunov
    P_st = local_dlyap(Tsr, B_Σ_Bt)
    P_st = (P_st + P_st') / 2

    # Sub-matrix for the ~39 needed variables (25× faster HP filter)
    needed_names = vcat(
        ["Y_$(i)"  for i in 1:nsec],
        ["PH_$(i)" for i in 1:nsec],
        ["L_$(i)"  for i in 1:nsec],
        ["GDP", "pi", "Q", "TB", "N"]   # N added for corr(N,GDP), corr(N,GDP/N)
    )
    needed_idx  = [get(endo_idx, nm, 0) for nm in needed_names]
    valid_mask  = needed_idx .> 0
    nidx_valid  = needed_idx[valid_mask]

    T_sub = Matrix{Float64}(ghx_jl[nidx_valid, :])   # n_needed × n_state
    R_sub = Matrix{Float64}(ghu_jl[nidx_valid, :])   # n_needed × n_exo
    RsubΣRsub = R_sub * Σe_smm * R_sub'
    RsubΣRsub = (RsubΣRsub + RsubΣRsub') / 2

    # HP-filtered covariance (spectral, λ=1600, 256 frequencies)
    w_var, w_lag1 = build_hp_weights(1600.0, 256)
    Γ_sub, Γ1_sub = hp_filtered_cov_fast(
        Matrix{Float64}(Tsr), B_Σ_Bt, T_sub, RsubΣRsub, w_var, w_lag1)

    # Map sub-matrix results back to needed_names ordering
    n_needed = length(needed_names)
    Γ_val  = zeros(n_needed, n_needed)
    Γ1_val = zeros(n_needed, n_needed)
    sub_positions = findall(valid_mask)
    for (si, pi) in enumerate(sub_positions), (sj, pj) in enumerate(sub_positions)
        Γ_val[pi, pj]  = Γ_sub[si, sj]
        Γ1_val[pi, pj] = Γ1_sub[si, sj]
    end

    # Local index lookup for needed_names
    ei_sub = Dict(nm => i for (i, nm) in enumerate(needed_names))

    # Helpers matching smm_estimation.jl exactly
    _pstd_hp(vn) = let k = get(ei_sub, vn, 0)
        k == 0 ? 0.0 : sqrt(max(Γ_val[k, k], 0.0)) / max(abs(ss_vec[needed_idx[k]]), 1e-12)
    end
    _xcorr_hp(v1, v2) = let i1 = get(ei_sub, v1, 0), i2 = get(ei_sub, v2, 0)
        (i1 == 0 || i2 == 0) ? NaN :
        let d = sqrt(max(Γ_val[i1,i1], 0.0) * max(Γ_val[i2,i2], 0.0))
            d < 1e-15 ? 0.0 : clamp(Γ_val[i1,i2] / d, -1.0, 1.0)
        end
    end

    # Sectoral moments
    for i in 1:nsec
        std_Y_m[i]  = _pstd_hp("Y_$(i)")
        std_PH_m[i] = _pstd_hp("PH_$(i)")
        std_L_m[i]  = _pstd_hp("L_$(i)")
        corr_YPH_m[i] = _xcorr_hp("Y_$(i)", "PH_$(i)")
    end

    # Employment comovement (same formulas as smm_model_moments.jl):
    # corr(N,GDP) directly; corr(N, GDP/N) via log-deviation (co)variances
    # (Γ is in LEVEL deviations → divide by steady states).
    m_corr_NGDP = _xcorr_hp("N", "GDP")
    m_corr_NAPL = let iN = get(ei_sub, "N", 0), iG = get(ei_sub, "GDP", 0)
        if iN == 0 || iG == 0
            NaN
        else
            Nbar = max(abs(ss_vec[needed_idx[iN]]), 1e-12)
            Gbar = max(abs(ss_vec[needed_idx[iG]]), 1e-12)
            v_n  = max(Γ_val[iN,iN], 0.0)/Nbar^2
            v_g  = max(Γ_val[iG,iG], 0.0)/Gbar^2
            c_ng = Γ_val[iN,iG]/(Nbar*Gbar)
            den  = sqrt(max(v_n, 0.0)*max(v_g + v_n - 2c_ng, 0.0))
            den < 1e-15 ? 0.0 : clamp((c_ng - v_n)/den, -1.0, 1.0)
        end
    end

    # Aggregate moments
    m_std_GDP    = _pstd_hp("GDP")
    m_std_pi     = _pstd_hp("pi")
    m_std_Q      = _pstd_hp("Q")
    m_corr_GDPpi = _xcorr_hp("GDP", "pi")
    m_corr_GDPQ  = _xcorr_hp("GDP", "Q")

    # std(TB/GDP): normalize by steady-state GDP (TB is a level variable)
    i_TB_sub = get(ei_sub, "TB", 0)
    GDP_ss_val = abs(ss_vec[get(endo_idx, "GDP", 1)])
    GDP_ss_val = GDP_ss_val > 1e-12 ? GDP_ss_val : 1.0
    m_std_TBGDP = (i_TB_sub > 0) ?
        sqrt(max(Γ_val[i_TB_sub, i_TB_sub], 0.0)) / GDP_ss_val : NaN

    # Autocorrelation of Q
    i_Q_sub = get(ei_sub, "Q", 0)
    m_autocorr_Q = (i_Q_sub > 0 && Γ_val[i_Q_sub, i_Q_sub] > 1e-15) ?
        Γ1_val[i_Q_sub, i_Q_sub] / Γ_val[i_Q_sub, i_Q_sub] : NaN

    # Rank correlations (default 0.0, matching estimation)
    vy = isfinite.(std_Y_m)  .& isfinite.(y_d)
    vp = isfinite.(std_PH_m) .& isfinite.(p_d)
    vl = isfinite.(std_L_m)  .& isfinite.(l_d)
    rho_y = sum(vy) >= 3 ? safe_spearman(std_Y_m[vy],  y_d[vy])  : 0.0
    rho_p = sum(vp) >= 3 ? safe_spearman(std_PH_m[vp], p_d[vp])  : 0.0
    rho_l = sum(vl) >= 3 ? safe_spearman(std_L_m[vl],  l_d[vl])  : 0.0

    rc_lyap_ok = true
catch e
    @printf "  [HP-filtered moments] failed (%s)\n" string(e)
end

@printf "\n--- Rank correlations (Spearman) ---\n"
@printf "  Output (Y) :   %7.4f\n" rho_y
@printf "  Prices (PH):   %7.4f\n" rho_p
@printf "  Labor (L)  :   %7.4f\n" rho_l

# =========================================================================== #
#  MOMENT FIT TABLE  (mirrors MATLAB run_log format)                          #
# =========================================================================== #

# Load data moments from aggregate_moments.csv + sectoral_moments.csv
agg_mom_path = joinpath(DATA_DIR, "aggregate_moments.csv")
sec_mom_path = joinpath(DATA_DIR, "sectoral_moments.csv")

d_std_GDP = NaN; d_std_pi = NaN; d_corr_GDPpi = NaN; d_omG = 0.57
d_std_Q = NaN; d_autocorr_Q = NaN; d_corr_GDPQ = NaN; d_TBGDP = NaN
d_std_TBGDP = NaN   # std of HP-filtered TB/GDP — add "std_TBGDP" to aggregate_moments.csv
d_corr_NGDP = NaN; d_corr_NAPL = NaN
y_d_tab = y_d; p_d_tab = p_d; l_d_tab = l_d

if isfile(agg_mom_path)
    agg = CSV.read(agg_mom_path, DataFrame)
    agg_d = Dict(String(r.moment) => Float64(r.value) for r in eachrow(agg))
    d_std_GDP    = get(agg_d, "std_GDP",     NaN)
    d_std_pi     = get(agg_d, "std_pi",      NaN)
    d_corr_GDPpi = get(agg_d, "corr_GDPpi", NaN)
    d_omG        = get(agg_d, "omG",         0.57)
    d_std_Q      = get(agg_d, "std_Q",       NaN)
    d_autocorr_Q = get(agg_d, "autocorr_Q", NaN)
    d_corr_GDPQ  = get(agg_d, "corr_GDPQ",  NaN)
    d_TBGDP      = get(agg_d, "TBGDP",      NaN)
    d_std_TBGDP  = get(agg_d, "std_TBGDP",  NaN)
    d_corr_NGDP  = get(agg_d, "corr_NGDP",  NaN)
    d_corr_NAPL  = get(agg_d, "corr_NAPL",  NaN)
end

# Build the FULL 60-element data and model vectors — same layout, same
# weighting matrix, same printer as the SMM estimator (all shared via
# utils.jl since 2026-07-10), so the loss printed here IS the SMM objective.
data_vec = [y_d_tab; p_d_tab; l_d_tab;
            d_std_GDP; d_std_pi; d_corr_GDPpi; d_std_TBGDP;
            d_std_Q; d_autocorr_Q; d_corr_GDPQ;
            1.0; 1.0; 1.0;             # rank corr targets = 1
            c_d;                        # corr(Y_i,PH_i), 12 sectors
            d_corr_NGDP; d_corr_NAPL]  # labor comovement

model_vec = [std_Y_m; std_PH_m; std_L_m;
             m_std_GDP; m_std_pi; m_corr_GDPpi; m_std_TBGDP;
             m_std_Q; m_autocorr_Q; m_corr_GDPQ;
             rho_y; rho_p; rho_l;
             corr_YPH_m;
             m_corr_NGDP; m_corr_NAPL]

@assert length(data_vec)  == N_MOMENTS "data_vec has $(length(data_vec)) ≠ N_MOMENTS=$N_MOMENTS"
@assert length(model_vec) == N_MOMENTS "model_vec has $(length(model_vec)) ≠ N_MOMENTS=$N_MOMENTS"

# NaN guard: replace NaN entries by 0 with a warning (a NaN data moment means
# aggregate/sectoral CSVs are stale — rerun compute_data_moments.jl)
for (v, nm) in ((data_vec, "data"), (model_vec, "model"))
    bad = findall(isnan, v)
    if !isempty(bad)
        @printf "  WARNING: NaN in %s moments at %s — set to 0 for the fit table\n" nm join(MOMENT_NAMES[bad], ", ")
        v[bad] .= 0.0
    end
end

W_smm  = build_weighting_matrix(data_vec)   # shared with smm_estimation.jl
W_diag = [W_smm[i,i] for i in 1:N_MOMENTS]
moment_labels = MOMENT_NAMES

# =========================================================================== #
#  LOAD SMM-REPORTED MOMENTS (authoritative, from smm_results.csv)            #
#  The estimation uses Klein's method for decision rules; Dynare's QZ can     #
#  produce numerically different T,R → different moments.  When available,    #
#  use the estimation's own moments so the loss matches exactly.              #
# =========================================================================== #

smm_res_file = joinpath(ESTIMATION_DIR, "smm_results.csv")
smm_model_vec = nothing   # will hold the 60-element model moments from estimation

# Only valid when θ came from smm_estimates.csv — smm_results.csv holds the
# moments of THAT completed run. A newer checkpoint θ has different moments.
if isfile(smm_res_file) && smm_param_source == "smm_estimates.csv"
    smm_res = CSV.read(smm_res_file, DataFrame)
    if hasproperty(smm_res, :model) && nrow(smm_res) == N_MOMENTS
        smm_model_vec = Float64.(smm_res.model)
        @printf "\n--- Loaded SMM-reported moments from smm_results.csv ---\n"
    elseif hasproperty(smm_res, :model)
        @printf "\n--- smm_results.csv has %d rows ≠ %d — stale layout, ignoring ---\n" nrow(smm_res) N_MOMENTS
    end
end

# Use estimation moments when available, otherwise fall back to Dynare-based
model_vec_final = something(smm_model_vec, model_vec)
model_source    = smm_model_vec !== nothing ? "SMM (Klein)" : "Dynare (QZ)"

# ---- Parameter report (θ actually used in this run, with bounds) ---------- #
θ_report = [ilabcosts_val; modepsY[1]; modepsM[1]; log(kappaV_val);
            rho_om1_val; rho_tfp1_val;
            Float64.(isigma_tfp_val);
            Float64.(sigma_om_vec);
            rho_pvstar_val; sigma_pvstar_val; rho_xi_val; sigma_xi_val;
            etastar_val; kappaw_val]
@printf "\n--- Parameters in effect (%s) ---\n" smm_param_source
print_param_table(θ_report)

# ---- Moment fit (same 60 moments, weights, and layout as the estimator) --- #
@printf "\n--- Moment fit (%s) ---\n" model_source
print_fit_table(data_vec, model_vec_final, W_smm)

total_loss = dot(data_vec .- model_vec_final, W_smm*(data_vec .- model_vec_final))

# Persist the same report to tables/moment_fit_<tag>.txt
mom_table_path = joinpath(TABLES_DIR, "moment_fit_$(tag).txt")
open(mom_table_path, "w") do f_mom
    write(f_mom, "NK-SOE Chile — Moment fit (Exercise: $(exercise_labels[EXERCISE+1]), source: $model_source)\n")
    write(f_mom, "Same 60 moments and weighting matrix as the SMM objective.\n")
    print_param_table(θ_report; io=f_mom)
    print_fit_table(data_vec, model_vec_final, W_smm; io=f_mom)
end
@printf "  → saved to: %s\n" mom_table_path

# If using SMM moments, also show Dynare-based moments for diagnostic comparison
if smm_model_vec !== nothing
    dynare_loss = dot(data_vec .- model_vec, W_smm*(data_vec .- model_vec))
    @printf "\n--- Dynare (QZ) diagnostic comparison ---\n"
    @printf "  Dynare QZ loss:  %11.6f\n" dynare_loss
    @printf "  SMM Klein loss:  %11.6f\n" total_loss
    @printf "  Difference:      %11.6f  (numerical: Klein vs QZ decomposition)\n" abs(dynare_loss - total_loss)
end

# Also save moment fit as CSV for easy analysis
df_mom = DataFrame(
    moment   = moment_labels,
    data     = data_vec,
    model    = model_vec_final,
    diff     = data_vec .- model_vec_final,
    W_diag   = W_diag,
    wdiff2   = W_diag .* (data_vec .- model_vec_final) .^ 2,
)
CSV.write(joinpath(TABLES_DIR, "moment_fit_$(tag).csv"), df_mom)

@printf "\n--- Dynare results ---\n"
@printf "  Steady state  : %s\n" (ss_ok ? "OK" : "FAILED")
@printf "  Blanchard-Kahn: %s\n" (rc_lyap_ok ? "OK" : "NOT CHECKED")
@printf "  Simulation    : %s\n" (sim_ok ? "OK" : "no data (ARM gees workaround)")
(rho_y == 0 || rho_p == 0 || rho_l == 0) &&
    @printf "  (0 = model std devs are uniform across sectors for this exercise)\n"
@printf "\n"

rank_corr = (rho_output=rho_y, rho_price=rho_p, rho_labor=rho_l)


# =========================================================================== #
#  SAVE RESULTS                                                                #
# =========================================================================== #

# All tables go to julia_dynare/tables/
output_path = joinpath(TABLES_DIR, "model_output_$(tag).csv")
sec_path    = joinpath(TABLES_DIR, "model_output_$(tag)_sectoral.csv")

# Steady-state scalars
df_ss = DataFrame(
    variable = ["GDP_ss","TB_ss","Bstar_ss","w_ss","Q_ss","C_ss",
                "rank_corr_Y","rank_corr_P","rank_corr_L"],
    value    = [GDP_ss, TB_ss, Bstar_ss, w_ss, Q_ss, C_ss, rho_y, rho_p, rho_l],
)
CSV.write(output_path, df_ss)

# Sectoral results
df_sec_out = DataFrame(
    sector  = 1:nsec,
    pH_ss   = pH_ss,  Yi_ss = Yi_ss, L_ss = L_ss, M_ss = M_ss, Vi_ss = Vi_ss,
    std_Y_m = std_Y_m, std_PH_m = std_PH_m, std_L_m = std_L_m,
    y_d = y_d, p_d = p_d, l_d = l_d,
)
CSV.write(sec_path, df_sec_out)
@printf "--- Tables saved to: %s\n\n" TABLES_DIR


# =========================================================================== #
#  FIGURES AND TABLES (equivalent to figs_SOE_gap.m + plot_*.m)              #
# =========================================================================== #

sec_results_for_figs = DataFrame(
    sector  = 1:nsec,
    pH_ss   = pH_ss,  Yi_ss = Yi_ss, L_ss = L_ss,
    std_Y   = std_Y_m, std_PH = std_PH_m, std_L = std_L_m,
)

# ---- figs_SOE_gap.jl: Dynare-style IRF plots → figures/irfs/ -----------
# These mirror what MATLAB Dynare would auto-generate:
#   irf_aggregate_<tag>_<shock>.png  — GDP, π, Q, TB per shock
#   irf_sectoral_Y/PH/L_<tag>.png   — 12-sector IRF panels
generate_figures(
    MOD_DIR        = MOD_DIR,
    FIGURES_DIR    = FIGURES_IRF_DIR,   # → figures/irfs/
    TABLES_DIR     = TABLES_DIR,
    SCRIPT_DIR     = SCRIPT_DIR,
    EXERCISE       = EXERCISE,
    nsec           = nsec,
    names_vec      = names_vec,
    ss_results     = (GDP_ss=GDP_ss, TB_ss=TB_ss, Q_ss=Q_ss,
                      C_ss=C_ss, N_ss=N_ss, Bstar_ss=Bstar_ss, w_ss=w_ss),
    sec_results    = sec_results_for_figs,
    exercise_label = exercise_labels[EXERCISE+1],
)

# ---- plot_scripts.jl: exercise analysis plots → figures/exercises/ -----
# Translates: figs_SOE_gap.m, plot_manufacturing_shock.m,
#             plot_figure7_manufacturing_shock.m, plot_shock_effects.m,
#             steady_state_table.m
irf_path = joinpath(MOD_DIR, "dynare_irfs.csv")
if isfile(irf_path) && filesize(irf_path) > 10
    df_irf_plots = CSV.read(irf_path, DataFrame)
    if nrow(df_irf_plots) > 0
        run_all_plots(
            df_irf      = df_irf_plots,
            EXERCISE    = EXERCISE,
            names_vec   = names_vec,
            ss_results  = (GDP_ss=GDP_ss, TB_ss=TB_ss, Q_ss=Q_ss,
                           C_ss=C_ss, N_ss=N_ss, Bstar_ss=Bstar_ss, w_ss=w_ss),
            sec_results = sec_results_for_figs,
            FIGURES_DIR = FIGURES_EX_DIR,   # → figures/exercises/
            TABLES_DIR  = TABLES_DIR,
            tag         = tag,
            ombar       = ombar_val,
        )
    end
end


# =========================================================================== #
#  DONE                                                                        #
# =========================================================================== #

@printf "%s\n" repeat("=", 60)
@printf "  Done: %s\n" exercise_labels[EXERCISE+1]
@printf "%s\n\n" repeat("=", 60)

end  # function _main()

Base.invokelatest(_main)  # Julia 1.12: world-age fix
