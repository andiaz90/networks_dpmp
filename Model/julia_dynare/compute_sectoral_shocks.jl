# compute_sectoral_shocks.jl
# ==========================
# Data-driven calibration of the 12 sectoral SUPPLY (TFP) and DEMAND (taste) shock
# sizes, to REPLACE their SMM estimation (redesign per FGI 2023 / FSW 2011 /
# Atalay 2017 — see docs/redesign_FSW_Atalay_memo.md).
#
# IDEA. A sector-specific SUPPLY (TFP) shock raises output and LOWERS its price;
# a sector-specific DEMAND (taste) shock raises BOTH. So the sign of the sectoral
# output–price correlation identifies the supply/demand mix:
#     corr(Y_i, P_i) = -1  ->  pure supply        corr(Y_i, P_i) = +1  ->  pure demand
# We split each sector's output variance accordingly (transparent first pass;
# refine with the common-factor removal + model inversion noted at the bottom):
#     frac_demand_i = (1 + r_i)/2,   frac_supply_i = (1 - r_i)/2
#     sigma_supply_i = sqrt(frac_supply_i) * std(Y_i)
#     sigma_demand_i = sqrt(frac_demand_i) * std(Y_i)
#
# This targets the corr(Y_i,PH_i) SIGN failure directly: the data's negative
# sectoral output-price correlations imply supply-dominated sectors, which the
# current demand/external-heavy estimate gets backwards.
#
# INPUT : Data/sectoral_moments.csv        (std_Y, std_PH, corr_YPH per sector)
# OUTPUT: Data/sectoral_shock_calibration.csv
# USAGE : julia --project=. compute_sectoral_shocks.jl
using CSV, DataFrames, Printf, Statistics

const SCRIPT_DIR = @__DIR__
const DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

infile = joinpath(DATA_DIR, "sectoral_moments.csv")
isfile(infile) || error("sectoral_moments.csv not found — run compute_data_moments.jl first.")
sec = CSV.read(infile, DataFrame)

n = nrow(sec)
frac_supply = fill(NaN, n); frac_demand = fill(NaN, n)
sig_supply  = fill(NaN, n); sig_demand  = fill(NaN, n)

for i in 1:n
    sY = Float64(sec.std_Y[i]); r = Float64(sec.corr_YPH[i])
    (isnan(sY) || isnan(r)) && continue
    r  = clamp(r, -0.99, 0.99)
    fs = (1 - r) / 2; fd = (1 + r) / 2
    frac_supply[i] = fs;            frac_demand[i] = fd
    sig_supply[i]  = sqrt(fs) * sY; sig_demand[i]  = sqrt(fd) * sY
end

out = DataFrame(
    sector          = sec.sector,
    name            = sec.name,
    std_Y_data      = sec.std_Y,
    corr_YPH_data   = sec.corr_YPH,
    frac_supply     = frac_supply,
    frac_demand     = frac_demand,
    sigma_supply    = sig_supply,   # target std of the sectoral SUPPLY (TFP) component of output
    sigma_demand    = sig_demand,   # target std of the sectoral DEMAND (taste) component of output
    # First-pass shock stds for the .mod (assumes unit shock->output loading).
    # Rescale with ONE model solve so model std(Y_i) matches std_Y_data (note 2).
    isigma_tfp_init = sig_supply,
    sigma_om_init   = sig_demand,
)

outfile = joinpath(DATA_DIR, "sectoral_shock_calibration.csv")
CSV.write(outfile, out)

@printf "\n  Sectoral supply/demand shock calibration (first pass)\n"
@printf "  %-16s %8s %9s %10s %10s\n" "sector" "std_Y" "corr_YP" "sig_supply" "sig_demand"
@printf "  %s\n" repeat("-", 58)
for i in 1:n
    @printf "  %-16s %8.3f %9.3f %10.4f %10.4f\n" sec.name[i] sec.std_Y[i] sec.corr_YPH[i] sig_supply[i] sig_demand[i]
end
@printf "\n  -> wrote %s\n\n" outfile

# ---------------------------------------------------------------------------
# REFINEMENTS (apply before treating these as the final calibration):
#  1. COMMON-FACTOR REMOVAL. std_Y here is TOTAL sectoral vol (idiosyncratic +
#     the common external/aggregate factor). To get the truly sector-specific
#     piece, regenerate the HP-filtered sectoral cycles Y~_i, P~_i, strip the
#     first principal component (the common factor), and recompute std_Y/corr_YPH
#     on the residual. Matters less once the external shock is reined in (see the
#     external-shock note), because then TOTAL ~= idiosyncratic.
#  2. MODEL INVERSION. isigma_tfp_init/sigma_om_init assume a unit shock->output
#     loading. Solve the model once at these values, read model std(Y_i), and
#     rescale:  isigma_tfp_i *= std_Y_data_i / std_Y_model_i  (same for sigma_om).
#     One or two iterations converges. Then FIX these in params_jl.mod and drop
#     the 24 sectoral shock variances from the estimated theta.
# ---------------------------------------------------------------------------
