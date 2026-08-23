# compute_sectoral_shocks.jl
# ==========================
# STAGE 1 of the data-driven calibration of the 12 sectoral SUPPLY (TFP) and
# DEMAND (taste) shock sizes, replacing their SMM estimation (redesign per
# FGI 2023 / FSW 2011 / Atalay 2017 — see docs/redesign_FSW_Atalay_memo.md).
#
#   Stage 1 (this file)  : analytic supply/demand variance split. THE BASELINE.
#   Stage 2 (calibrate_sectoral_shocks.jl) : model-inversion rescaling. RETIRED
#                          2026-08-21 — see the HISTORY note at the bottom.
#
# This file is now the last step before estimation. Run it, then
# run_smm_estimation.jl (which defaults to SHOCK_STAGE=1).
#
# IDEA. A sector-specific SUPPLY (TFP) shock raises output and LOWERS its price;
# a sector-specific DEMAND (taste) shock raises BOTH. So the sign of the sectoral
# output–price correlation identifies the supply/demand mix:
#     corr(Y_i, P_i) = -1  ->  pure supply        corr(Y_i, P_i) = +1  ->  pure demand
# We split each sector's output variance accordingly:
#     frac_demand_i = (1 + r_i)/2,   frac_supply_i = (1 - r_i)/2
#     sigma_supply_i = sqrt(frac_supply_i) * std(Y_i)
#     sigma_demand_i = sqrt(frac_demand_i) * std(Y_i)
#
# WHICH MOMENTS GO IN (changed 2026-08-21 — this is the substantive fix).
# The inputs are the IDIOSYNCRATIC moments std_Y_idio and corr_YPH_idio, i.e.
# computed on the sectoral cycles AFTER the common factor has been projected
# out (compute_data_moments.jl section 9c). They are NOT the total std_Y /
# corr_YPH, which is what this script used until today.
#
#   Using TOTAL volatility double-counts. Total sectoral volatility is
#   idiosyncratic + common. The model generates the common part from its OWN
#   aggregate and external shocks (pvstar, zeta, monetary), which are calibrated
#   separately. Sizing the SECTORAL shocks to reproduce TOTAL volatility
#   therefore asks the sectoral shocks to do work the aggregate shocks are
#   already doing, and the model over-predicts sectoral volatility.
#
#   Using the TOTAL correlation also mis-signs the mix. A sector whose output
#   and price comove positively because BOTH are driven by the aggregate cycle
#   is not a sector with a large idiosyncratic demand shock. Empirically this
#   matters a lot: projecting out one common factor flips corr(Y,PH) from
#   positive to negative for Construction, Trade/Hotels, Real Estate and
#   Business Services.
#
# INPUT : Data/sectoral_moments.csv   (std_Y_idio, corr_YPH_idio per sector;
#                                      falls back to std_Y/corr_YPH with a loud
#                                      warning if the _idio columns are absent)
# OUTPUT: Data/sectoral_shock_calibration.csv  AND an identical copy at
#         Data/sectoral_shock_calibration_stage1.csv, which is what
#         run_smm_estimation.jl reads under its default SHOCK_STAGE=1.
# USAGE : julia --project=. compute_sectoral_shocks.jl
using CSV, DataFrames, Printf, Statistics

const SCRIPT_DIR = @__DIR__
const DATA_DIR   = joinpath(abspath(SCRIPT_DIR, "..", ".."), "Data")

infile = joinpath(DATA_DIR, "sectoral_moments.csv")
isfile(infile) || error("sectoral_moments.csv not found — run compute_data_moments.jl first.")
sec = CSV.read(infile, DataFrame)

# --- Choose the inversion inputs: idiosyncratic if available --------------- #
const USE_IDIO = hasproperty(sec, :std_Y_idio) && hasproperty(sec, :corr_YPH_idio) &&
                 !all(isnan, sec.std_Y_idio) && !all(isnan, sec.corr_YPH_idio)

sY_in, r_in, src = if USE_IDIO
    Float64.(sec.std_Y_idio), Float64.(sec.corr_YPH_idio), "IDIOSYNCRATIC (common factor removed)"
else
    @printf "\n  %s\n" repeat("!", 70)
    @printf "  WARNING: sectoral_moments.csv has no usable std_Y_idio / corr_YPH_idio.\n"
    @printf "  Falling back to TOTAL volatility. This DOUBLE-COUNTS the common factor\n"
    @printf "  and mis-signs the supply/demand mix — see the header of this file.\n"
    @printf "  Fix: delete Data/sectoral_moments.csv and re-run compute_data_moments.jl.\n"
    @printf "  %s\n\n" repeat("!", 70)
    hasproperty(sec, :corr_YPH) || error("""
        sectoral_moments.csv has neither corr_YPH_idio nor corr_YPH. It was
        almost certainly written by bootstrap_csv.jl, which emits only
        sector/std_Y/std_PH/std_L. There is nothing to invert the supply/demand
        mix from. Run: julia --project=. compute_data_moments.jl
        """)
    Float64.(sec.std_Y), Float64.(sec.corr_YPH), "TOTAL (FALLBACK — see warning above)"
end

n = nrow(sec)

# bootstrap_csv.jl writes a minimal sectoral_moments.csv (sector/std_Y/std_PH/
# std_L only — no name, no corr_YPH, no _idio). Don't hard-fail on it.
sec_names = hasproperty(sec, :name) ? String.(sec.name) :
            ["Agriculture","Mining","Manufacturing","Utilities","Construction",
             "Trade/Hotels","Transport/Comm","Finance","Real Estate",
             "Business Serv.","Personal Serv.","Public Admin."][1:n]

frac_supply = fill(NaN, n); frac_demand = fill(NaN, n)
sig_supply  = fill(NaN, n); sig_demand  = fill(NaN, n)

for i in 1:n
    sY = sY_in[i]; r = r_in[i]
    (isnan(sY) || isnan(r)) && continue
    r  = clamp(r, -0.99, 0.99)
    fs = (1 - r) / 2; fd = (1 + r) / 2
    frac_supply[i] = fs;            frac_demand[i] = fd
    sig_supply[i]  = sqrt(fs) * sY; sig_demand[i]  = sqrt(fd) * sY
end

out = DataFrame(
    sector          = sec.sector,
    name            = sec_names,
    std_Y_data      = sY_in,          # the series actually inverted (idio or total)
    corr_YPH_data   = r_in,
    std_Y_total     = Float64.(sec.std_Y),   # kept for reference / the Stage-2 target
    frac_supply     = frac_supply,
    frac_demand     = frac_demand,
    sigma_supply    = sig_supply,   # target std of the sectoral SUPPLY (TFP) component of output
    sigma_demand    = sig_demand,   # target std of the sectoral DEMAND (taste) component of output
    # Shock stds for the .mod. These assume a unit shock->output loading, which
    # is what Stage 2 was built to correct — and, empirically, correcting it
    # made the fit worse rather than better (see HISTORY). lambda_A and
    # lambda_om absorb the residual scale instead, and lambda_A estimates at
    # 0.836, so the unit-loading assumption is not far off in aggregate.
    isigma_tfp_init = sig_supply,
    sigma_om_init   = sig_demand,
    # Provenance so a downstream reader can tell a Stage-1 file from a Stage-2 one.
    stage           = fill("1-analytic", n),
    moment_source   = fill(src, n),
)

outfile = joinpath(DATA_DIR, "sectoral_shock_calibration.csv")
CSV.write(outfile, out)
# The estimator defaults to SHOCK_STAGE=1 and reads the _stage1 copy, so write
# it here. Without this, a fresh checkout has no Stage-1 file and the default
# path errors — the _stage1 name used to be created only as a Stage-2 backup.
outfile_s1 = joinpath(DATA_DIR, "sectoral_shock_calibration_stage1.csv")
CSV.write(outfile_s1, out)

@printf "\n  Sectoral supply/demand shock calibration — STAGE 1 (analytic split)\n"
@printf "  Inversion inputs: %s\n\n" src
@printf "  %-16s %8s %9s %10s %10s\n" "sector" "std_Y" "corr_YP" "sig_supply" "sig_demand"
@printf "  %s\n" repeat("-", 58)
for i in 1:n
    @printf "  %-16s %8.4f %9.3f %10.4f %10.4f\n" sec_names[i] sY_in[i] r_in[i] sig_supply[i] sig_demand[i]
end
@printf "  %s\n" repeat("-", 58)
@printf "  demand share of sectoral variance: mean %.2f, range [%.2f, %.2f]\n" (
    mean(filter(!isnan, frac_demand))) (minimum(filter(!isnan, frac_demand))) (
    maximum(filter(!isnan, frac_demand)))
@printf "\n  -> wrote %s\n" outfile
@printf "  -> wrote %s  (read by SHOCK_STAGE=1, the default)\n" outfile_s1
@printf "  -> This is the BASELINE calibration. Stage 2 (calibrate_sectoral_shocks.jl)\n"
@printf "     is retired: it scored 14.75 against Stage 1's 11.46 and drove 10 of 24\n"
@printf "     sigmas onto their bounds. Kept only for the robustness appendix.\n"
@printf "  -> NEXT: julia --threads=1 --project=. run_smm_estimation.jl\n\n"

# ---------------------------------------------------------------------------
# HISTORY — why Stage 2 was built, and why it was retired.
#
# Until 2026-08-21 this script's output was used directly as the pinned
# sectoral shock sizes in the SMM (theta[7:18], theta[19:30]), and the two
# refinements listed here were never applied. The 2026-08-21 estimation run
# shows what that cost:
#
#   * lambda_A = 0.59 and lambda_om = 3.14 — the two GLOBAL shock-mix scales,
#     theta[37:38] — were the only parameters with meaningful sensitivity
#     (|d obj| ~ 10.3, against 0.45 for cl). They were doing Stage 2's job with
#     two scalars where 24 sector-specific ones are needed. lambda_A < 1 is the
#     estimator undoing the common-factor double-counting (refinement 1);
#     lambda_om = 3.14 is it fighting the unit-loading assumption (refinement 2).
#   * corr(Y_i,PH_i) was still 30% of the objective even though corr_YPH is the
#     exact statistic used to do the split — because the split is
#     partial-equilibrium and network spillovers break it in GE.
#
# Refinement 1 (common-factor removal) is now done upstream, in
# compute_data_moments.jl section 9c. Refinement 2 (model inversion) is Stage 2,
# calibrate_sectoral_shocks.jl. Once both have run, lambda_A and lambda_om
# should sit near 1; if they do not, the inversion has not converged and the
# Stage-2 diagnostics are the place to look.
# ---------------------------------------------------------------------------
#
# VERDICT 2026-08-21. Stage 2 was tried three ways and retired.
#   target = model TOTAL volatility   -> sigma -> 0 (spillover floor), obj 14.75
#   target = own-shock-only volatility -> sigma -> UB (std_Y_idio is not an
#                                         own-shock object; PC1 removes the
#                                         COMMON factor, not sectoral spillovers)
#   Stage 1 alone, no inversion        -> obj 11.46, every sigma interior
# Stage 1 wins on fit AND on credibility. lambda_A = 0.836 says the unit-loading
# assumption is close enough in aggregate that two global scales handle the rest.
# calibrate_sectoral_shocks.jl is kept, working, for the robustness appendix.
