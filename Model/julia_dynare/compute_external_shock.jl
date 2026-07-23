# compute_external_shock.jl
# =========================
# Data-driven calibration of the external (import-price / terms-of-trade) shock
# process — rho_pvstar and sigma_pvstar — to REPLACE its SMM estimation, so the
# external shock is disciplined by an OBSERVED external-price series instead of by
# domestic RER/TB moments (Schmitt-Grohé & Uribe 2018; García-Cicco-Pancrazi-Uribe
# 2010; BCCh XMAS). See docs/redesign_FSW_Atalay_memo.md and the SGU (2018) wrapper.
#
# Model process:  log(PVstar_t / PVstar_ss) = rho * log(PVstar_{t-1}/PVstar_ss)
#                                              + sigma * eps_t.
# We fit that AR(1) to the (log, HP-detrended) observed series.
#
# INPUT : a 2-column CSV in Data/ with a date column and a price/index column.
#         For Chile, natural anchors (best first):
#           - BCCh terms-of-trade index (Términos de intercambio), quarterly
#           - copper price (BCCh precio del cobre / LME / FRED PCOPPUSDM -> quarterly avg)
#           - a Chilean import price / import unit-value index
# OUTPUT: Data/external_shock_calibration.csv  (rho_pvstar, sigma_pvstar, uncond std)
# USAGE : EXT_FILE=tot_chile_quarterly.csv [EXT_COL=tot] julia --project=. compute_external_shock.jl
using CSV, DataFrames, Printf, Statistics, SparseArrays, LinearAlgebra

const SCRIPT_DIR = @__DIR__
const DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")
const LAMBDA     = 1600.0   # HP smoothing (quarterly), matches compute_data_moments.jl

EXT_FILE = get(ENV, "EXT_FILE", "tot_chile_quarterly.csv")  # put your series here
EXT_COL  = get(ENV, "EXT_COL",  "")                         # value column; "" = 2nd column

# HP cycle (identical to compute_data_moments.jl)
function hp_cycle(y::AbstractVector{<:Real}, λ::Real)
    T = length(y); T < 4 && return zeros(T)
    m = T-2; ri = vcat(1:m,1:m,1:m); ci = vcat(1:m,2:m+1,3:m+2)
    val = vcat(fill(1.0,m), fill(-2.0,m), fill(1.0,m)); D = sparse(ri,ci,val,m,T)
    H = sparse(1.0I,T,T) + λ*(D'D)
    collect(Float64,y) .- (H \ collect(Float64,y))
end

path = joinpath(DATA_DIR, EXT_FILE)
isfile(path) || error("""
  External series not found: $path
  Drop a 2-column CSV (date, value) into Data/. Suggested sources:
    - BCCh Base de Datos Estadísticos (si3.bcentral.cl) -> "Términos de intercambio" (quarterly)
    - copper price: BCCh (precio del cobre BML), FRED PCOPPUSDM (monthly -> quarterly avg),
      or World Bank Pink Sheet.
  Then re-run:  EXT_FILE=<file>.csv [EXT_COL=<col>] julia --project=. compute_external_shock.jl
  """)

df = CSV.read(path, DataFrame)
datecol = names(df)[1]
col     = isempty(EXT_COL) ? names(df)[2] : EXT_COL
# Restrict to the estimation sample (matches compute_data_moments.jl: 2006Q1–2023Q4).
# Override with EXT_START / EXT_END (ISO yyyy-mm-dd) to use a different window.
EXT_START = get(ENV, "EXT_START", "2006-01-01")
EXT_END   = get(ENV, "EXT_END",   "2023-12-31")
ds   = string.(df[!, datecol])
keep = (ds .>= EXT_START) .& (ds .<= EXT_END)
x    = Float64.(collect(skipmissing(df[keep, col])))
length(x) < 20 && error("Only $(length(x)) obs in [$EXT_START, $EXT_END] — widen EXT_START/EXT_END or check the file.")
all(x .> 0) || error("Series has non-positive values; expected a price/index level for log().")
@printf "  Sample window: %s .. %s  (%d quarters)\n" EXT_START EXT_END length(x)

# log -> HP-detrend -> AR(1)
cyc = hp_cycle(log.(x), LAMBDA)               # stationary cyclical component of the log price
y0  = cyc[2:end]; y1 = cyc[1:end-1]
rho    = (y1'y0) / (y1'y1)                     # OLS AR(1): cyc_t = rho*cyc_{t-1} + e_t
sigma  = std(y0 .- rho .* y1)                  # innovation std = sigma_pvstar
uncond = sigma / sqrt(max(1 - rho^2, 1e-8))

@printf "\n  External-shock calibration from %s (col '%s', %d quarters)\n" EXT_FILE col length(x)
@printf "  %s\n" repeat("-", 60)
@printf "  rho_pvstar   (AR1 persistence)   = %.4f\n" rho
@printf "  sigma_pvstar (innovation std)    = %.4f\n" sigma
@printf "  implied unconditional std        = %.4f\n" uncond
@printf "  (compare current SMM estimate:  rho ~ 0.93, sigma ~ 0.195, uncond ~ 0.54)\n\n"

CSV.write(joinpath(DATA_DIR, "external_shock_calibration.csv"),
          DataFrame(param = ["rho_pvstar","sigma_pvstar","uncond_std"],
                    value = [rho, sigma, uncond]))
@printf "  -> wrote Data/external_shock_calibration.csv\n"
@printf "  Next: set rho_pvstar/sigma_pvstar to these values in the calibration and\n"
@printf "        drop theta[31], theta[32] from the estimated set.\n\n"

# NOTE. Fitting AR(1) to the HP-cycle gives the persistence/vol of the cyclical
# (business-cycle-frequency) component — consistent with the HP-filtered moments
# the model targets. For an exact match one could instead target the HP-filtered
# model PVstar moment directly; the AR(1)-on-cycle fit is a transparent first pass.
