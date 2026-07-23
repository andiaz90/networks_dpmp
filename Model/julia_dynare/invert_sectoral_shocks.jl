# invert_sectoral_shocks.jl
# =========================
# Redesign Step 3 — LOADING INVERSION. Rescale the sectoral shock sizes so the
# model reproduces each sector's DATA output volatility. Run AFTER main_SOE_gap.jl
# (which writes tables/model_output_ex2_mfg_sectoral.csv with std_Y_m and y_d),
# then re-run main. Repeat 1-2x until std_Y_m ~ y_d.
# USAGE: julia --project=. invert_sectoral_shocks.jl
using CSV, DataFrames, Printf
const DATA_DIR = joinpath(@__DIR__, "..", "..", "Data")
const MOD_TBL  = joinpath(@__DIR__, "tables", "model_output_ex2_mfg_sectoral.csv")
isfile(MOD_TBL) || error("run main_SOE_gap.jl first (missing $MOD_TBL)")
m   = CSV.read(MOD_TBL, DataFrame)
cal_path = joinpath(DATA_DIR, "sectoral_shock_calibration.csv")
cal = CSV.read(cal_path, DataFrame)
n = nrow(cal)
@printf "\n  Sectoral loading inversion (rescale to hit data std_Y)\n"
@printf "  %-16s %9s %9s %8s\n" "sector" "std_Y_m" "y_d" "factor"
@printf "  %s\n" repeat("-", 46)
for i in 1:n
    mm = Float64(m.std_Y_m[i]); dd = Float64(m.y_d[i])
    if isnan(mm) || isnan(dd) || mm <= 0
        @printf "  %-16s   skip (bad std)\n" cal.name[i]; continue
    end
    f = clamp(dd / mm, 0.1, 10.0)
    cal.isigma_tfp_init[i] *= f
    cal.sigma_om_init[i]   *= f
    @printf "  %-16s %9.4f %9.4f %8.3f\n" cal.name[i] mm dd f
end
CSV.write(cal_path, cal)
@printf "\n  -> updated %s\n  Re-run main_SOE_gap.jl; repeat if std_Y_m still far from y_d.\n\n" cal_path
