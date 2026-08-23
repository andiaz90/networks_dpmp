"""
figs_SOE_gap.jl
===============
Julia equivalent of figs_SOE_gap.m + plot_manufacturing_shock.m +
plot_shock_effects.m — generates IRF figures and a steady-state table
for the NK-IOSOE Chile model.

Called from main_SOE_gap.jl (included at top level, so no world-age issue).
Reads dynare_irfs.csv (written by run_dynare_subprocess.jl).

PDF figures saved to DATA_DIR (requires Plots.jl: ] add Plots).
"""

using CSV, DataFrames, Printf, Statistics

# Load Plots.jl once at module level (graceful failure if not installed)
const _PLOTS_OK = Ref(false)
try
    @eval Main using Plots
    @eval Main gr()
    _PLOTS_OK[] = true
catch
end

# =========================================================================== #
#  HELPERS                                                                     #
# =========================================================================== #

# Accept AbstractString so InlineString types (String15, String31, etc.)
# from CSV.jl work without conversion — compare as strings
function get_irf(df::DataFrame, varname::AbstractString, shock::AbstractString; n_periods::Int=40)
    sub = filter(r -> String(r.variable) == String(varname) && String(r.shock) == String(shock), df)
    isempty(sub) && return zeros(n_periods)
    s = sort(sub, :period)
    n = min(n_periods, nrow(s))
    out = zeros(n_periods)
    out[1:n] .= s.value[1:n]
    return out
end

# Silent on success (a 32-figure loop printed 32 full paths and drowned the
# output, 2026-07-10); failures still print. _N_FIGS_SAVED lets the caller
# print one summary line ("32 figures → dir") instead.
const _N_FIGS_SAVED = Ref(0)
function _savefig_safe(p, path::String)
    try
        @eval Main savefig($p, $path)
        _N_FIGS_SAVED[] += 1
    catch e
        @printf "  WARNING: could not save %s — %s\n" path string(e)
    end
end

# =========================================================================== #
#  MAIN ENTRY POINT                                                            #
# =========================================================================== #

function generate_figures(;
    MOD_DIR::String,
    FIGURES_DIR::String,   # julia_dynare/figures/
    TABLES_DIR::String,    # julia_dynare/tables/
    SCRIPT_DIR::String,
    EXERCISE::Int,
    nsec::Int,
    names_vec,
    ss_results::NamedTuple,
    sec_results::DataFrame,
    exercise_label::String,
)
    @printf "\n--- Generating figures (Exercise %d) ---\n" EXERCISE

    # ---- Load IRF CSV ---------------------------------------------------- #
    irf_file = joinpath(MOD_DIR, "dynare_irfs.csv")
    if !isfile(irf_file) || filesize(irf_file) < 10
        @printf "  No IRF file found at %s\n" irf_file
        @printf "  (Run the model first to generate IRFs)\n"
        return
    end

    df_irf = CSV.read(irf_file, DataFrame)
    if nrow(df_irf) == 0
        @printf "  IRF file is empty — Dynare stoch_simul produced no IRFs.\n"
        @printf "  This usually means the model SS was not computed correctly.\n"
        return
    end

    all_shocks = unique(String.(df_irf.shock))
    n_periods  = 40   # show first 40 quarters in plots

    @printf "  IRFs loaded: %d shock(s), %d variable-shock pairs\n" length(all_shocks) length(unique(df_irf.variable))

    # For EXERCISE 2 (manufacturing), focus on epsA_3
    main_shock = if EXERCISE == 2
        s = filter(x -> contains(x, "epsA_3"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    elseif EXERCISE == 1
        s = filter(x -> contains(x, "eps_om") || contains(x, "eps_zeta"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    elseif EXERCISE == 3
        s = filter(x -> contains(x, "eps_i"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    else
        # Baseline: use epsA_3 (manufacturing) for sectoral plots
        s = filter(x -> contains(x, "epsA_3"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    end

    @printf "  Main shock for sectoral plots: %s\n" main_shock

    # ---- 1. AGGREGATE IRF TABLE ----------------------------------------- #
    # EInv is investment expenditure, the third leg of GDP = C + EInv + TB.
    agg_vars   = ["GDP", "EInv", "pi", "Q", "TB"]
    agg_labels = ["GDP (%%)", "Investment (%%)", "Inflation π (%%)",
                  "Real XR Q (%%)", "Trade balance TB (%%)"]

    @printf "\n  Aggregate IRFs — %s (period 1 impact, %% dev. from SS):\n" main_shock
    @printf "  %-8s" "Shock"
    for lbl in agg_labels; @printf "  %16s" lbl; end
    println()
    @printf "  %s\n" repeat("-", 72)

    for shock in all_shocks
        @printf "  %-8s" shock[1:min(8,length(shock))]
        for vn in agg_vars
            irf = get_irf(df_irf, vn, shock; n_periods=n_periods)
            @printf "  %16.4f" irf[1]
        end
        println()
    end

    # ---- 2. SECTORAL OUTPUT IMPACT TABLE --------------------------------- #
    # U_i is capital utilisation K_i/Kbar_i; RK_i the rental. Both are part of
    # the transmission: with nu = 0.4288 capital services respond to the rental,
    # so a sector's marginal cost rises less than under a pure fixed factor.
    @printf "\n  Sectoral output IRFs (period 1 impact, %% dev. from SS) — shock: %s\n" main_shock
    @printf "  %-6s  %-32s  %10s  %10s  %10s  %10s  %10s\n" "Sector" "Name" "ΔY" "ΔPH" "ΔL" "ΔU" "ΔRK"
    @printf "  %s\n" repeat("-", 98)
    for i in 1:nsec
        nm  = string(names_vec[i])[1:min(30, length(string(names_vec[i])))]
        dy  = get_irf(df_irf, "Y_$(i)",  main_shock; n_periods=n_periods)[1]
        dph = get_irf(df_irf, "PH_$(i)", main_shock; n_periods=n_periods)[1]
        dl  = get_irf(df_irf, "L_$(i)",  main_shock; n_periods=n_periods)[1]
        du  = get_irf(df_irf, "U_$(i)",  main_shock; n_periods=n_periods)[1]
        drk = get_irf(df_irf, "RK_$(i)", main_shock; n_periods=n_periods)[1]
        @printf "  %-6d  %-32s  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n" i nm dy dph dl du drk
    end

    # ---- 3. STEADY-STATE SUMMARY TABLE ---------------------------------- #
    @printf "\n  STEADY-STATE SUMMARY\n"
    @printf "  %-30s  %12s\n" "Variable" "Value"
    @printf "  %s\n" repeat("-", 45)
    ss_rows = [
        ("GDP",                     ss_results.GDP_ss),
        ("TB/GDP",                  ss_results.TB_ss / max(ss_results.GDP_ss, 1e-10)),
        ("Real exchange rate (Q)",  ss_results.Q_ss),
        ("Consumption (C)",         ss_results.C_ss),
        ("Labor (N)",               ss_results.N_ss),
        # NOMINAL, i.e. Σ pH_i·Y_i / GDP (2026-08-19). This printed the physical
        # sum Σ Y_i / GDP, so the same run reported two different gross-output
        # ratios — 1.830 in the steady-state accounting block against 2.227
        # here. Only the nominal one is comparable to the national accounts
        # (Chile 2021: 1.924). Matches ygdp_ss in main_SOE_gap.jl:766.
        ("Gross output/GDP",        sum(sec_results.pH_ss .* sec_results.Yi_ss) /
                                    max(ss_results.GDP_ss, 1e-10)),
        ("Foreign debt/GDP",        ss_results.Q_ss * abs(ss_results.Bstar_ss) / max(ss_results.GDP_ss, 1e-10)),
    ]
    # Factor shares and investment. The labour share of value added is the
    # headline number the paper quotes: 1.000 in a model without capital, 0.43
    # here, against 0.416 in the Chilean accounts.
    let
        _VA = sum(sec_results.pH_ss .* sec_results.Yi_ss) -
              sum(sec_results.PMi_ss .* sec_results.M_ss) -
              ss_results.PV_ss * sum(sec_results.Vi_ss)
        append!(ss_rows, [
            ("Investment (EInv)",       ss_results.EInv_ss),
            ("Investment/GDP",          ss_results.EInv_ss / max(ss_results.GDP_ss, 1e-10)),
            ("  data FBCF/VA",          0.1704),
            ("Labour share of VA",      sum(sec_results.PL_ss .* sec_results.L_ss) / max(_VA, 1e-10)),
            ("  data REM/VA",           0.416),
            ("Capital share of VA",     sum(sec_results.RKss_ss .* sec_results.Kbar_ss) / max(_VA, 1e-10)),
            ("  data EBE/VA",           0.568),
        ])
    end
    for (lbl, val) in ss_rows
        @printf "  %-30s  %12.4f\n" lbl val
    end
    @printf "\n  %-6s  %-30s  %8s  %8s  %8s\n" "Sector" "Name" "Y_ss" "L_ss" "PH_ss"
    @printf "  %s\n" repeat("-", 62)
    for i in 1:nsec
        r  = sec_results[i, :]
        nm = string(names_vec[i])[1:min(28, length(string(names_vec[i])))]
        @printf "  %-6d  %-30s  %8.4f  %8.4f  %8.4f\n" i nm r.Yi_ss r.L_ss r.pH_ss
    end

    # ---- 4. PDF FIGURES ------------------------------------------------- #
    if !_PLOTS_OK[]
        @printf "\n  Plots.jl not available — skipping PDF figures.\n"
        @printf "  Install with:  using Pkg; Pkg.add(\"Plots\")\n"
        @printf "\n--- Figures complete (text output only) ---\n"
        return
    end

    periods = 1:n_periods
    tag     = "ex$(EXERCISE)"

    # Human-readable description of each structural shock. Previously the
    # per-shock aggregate IRF title reused `exercise_label` (e.g. "TFP shock to
    # Manufacturing") for EVERY shock, which was misleading — each figure is the
    # model's response to its OWN shock. Map the shock name to a real label.
    shock_desc(shk) = begin
        if     shk == "eps_i";       "Monetary policy shock"
        elseif shk == "epschi";      "Labor-supply shock"
        elseif shk == "eps_pvstar";  "Import-price shock"
        elseif shk == "eps_postar";  "Oil-price shock"
        elseif shk == "eps_zeta";      "Aggregate demand (preference) shock"
        elseif startswith(shk, "epsA_")
            i = tryparse(Int, replace(shk, "epsA_" => ""))
            (i !== nothing && 1 <= i <= length(names_vec)) ? "TFP shock: $(names_vec[i])" : "TFP shock ($shk)"
        elseif startswith(shk, "eps_om_")
            i = tryparse(Int, replace(shk, "eps_om_" => ""))
            (i !== nothing && 1 <= i <= length(names_vec)) ? "Demand (taste) shock: $(names_vec[i])" : "Demand shock ($shk)"
        else
            shk
        end
    end

    # 4a. Aggregate IRFs — one PNG per shock.
    # The layout is DERIVED from length(agg_vars), not hardcoded. It used to be
    # (2,2), which silently broke on 2026-08-20 when EInv (investment) became a
    # fifth aggregate: five series into four slots made Plots draw the last panel
    # over the figure title. agg_vars feeds both the table above and this loop,
    # so anything added there must widen the grid here too — deriving it removes
    # the coupling.
    _nagg  = length(agg_vars)
    _ncol  = _nagg <= 4 ? 2 : 3
    _nrow  = ceil(Int, _nagg / _ncol)
    for shock in all_shocks
        irfs = [get_irf(df_irf, vn, shock; n_periods=n_periods) for vn in agg_vars]
        p = Plots.plot(layout=(_nrow,_ncol), size=(300*_ncol + 100, 300*_nrow),
                       titlefontsize=9,
                       plot_title="$shock — $(shock_desc(shock))")
        for (k, (irf, lbl)) in enumerate(zip(irfs, agg_labels))
            Plots.plot!(p, periods, irf, subplot=k, label="",
                       color=:steelblue, lw=2, title=lbl,
                       xlabel="Quarters", ylabel="% dev. from SS")
            Plots.hline!(p, [0.0], subplot=k,
                        color=:black, lw=0.8, ls=:dash, label="")
        end
        fname = joinpath(FIGURES_DIR, "irf_aggregate_$(tag)_$(shock).png")
        _savefig_safe(p, fname)
    end

    # 4b. Sectoral output IRF — 4×3 subplot layout for main shock
    p_Y = Plots.plot(layout=(4,3), size=(1200,900), titlefontsize=8,
                     plot_title="Sectoral Output — $main_shock")
    for i in 1:nsec
        irf = get_irf(df_irf, "Y_$(i)", main_shock; n_periods=n_periods)
        nm  = string(names_vec[i])[1:min(14, length(string(names_vec[i])))]
        Plots.plot!(p_Y, periods, irf, subplot=i, title=nm, label="",
                   color=:steelblue, lw=2, xlabel="Q", ylabel="% dev.")
        Plots.hline!(p_Y, [0.0], subplot=i,
                    color=:black, lw=0.5, ls=:dash, label="")
    end
    _savefig_safe(p_Y, joinpath(FIGURES_DIR, "irf_sectoral_Y_$(tag).png"))

    # 4c. Sectoral prices
    p_PH = Plots.plot(layout=(4,3), size=(1200,900), titlefontsize=8,
                      plot_title="Sectoral Prices — $main_shock")
    for i in 1:nsec
        irf = get_irf(df_irf, "PH_$(i)", main_shock; n_periods=n_periods)
        nm  = string(names_vec[i])[1:min(14, length(string(names_vec[i])))]
        Plots.plot!(p_PH, periods, irf, subplot=i, title=nm, label="",
                   color=:firebrick, lw=2, xlabel="Q", ylabel="% dev.")
        Plots.hline!(p_PH, [0.0], subplot=i,
                    color=:black, lw=0.5, ls=:dash, label="")
    end
    _savefig_safe(p_PH, joinpath(FIGURES_DIR, "irf_sectoral_PH_$(tag).png"))

    # 4d. Sectoral employment
    p_L = Plots.plot(layout=(4,3), size=(1200,900), titlefontsize=8,
                     plot_title="Sectoral Employment — $main_shock")
    for i in 1:nsec
        irf = get_irf(df_irf, "L_$(i)", main_shock; n_periods=n_periods)
        nm  = string(names_vec[i])[1:min(14, length(string(names_vec[i])))]
        Plots.plot!(p_L, periods, irf, subplot=i, title=nm, label="",
                   color=:darkgreen, lw=2, xlabel="Q", ylabel="% dev.")
        Plots.hline!(p_L, [0.0], subplot=i,
                    color=:black, lw=0.5, ls=:dash, label="")
    end
    _savefig_safe(p_L, joinpath(FIGURES_DIR, "irf_sectoral_L_$(tag).png"))

    # 4e. Output gap (Ygap variables — deviation of NK from flex-price)
    gap_vars_exist = any(r -> String(r.variable) == "Ygap_1", eachrow(df_irf))
    if gap_vars_exist
        p_gap = Plots.plot(layout=(4,3), size=(1200,900), titlefontsize=8,
                           plot_title="Output Gap — $main_shock")
        for i in 1:nsec
            irf = get_irf(df_irf, "Ygap_$(i)", main_shock; n_periods=n_periods)
            nm  = string(names_vec[i])[1:min(14, length(string(names_vec[i])))]
            Plots.plot!(p_gap, periods, irf, subplot=i, title=nm, label="",
                       color=:darkorange, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p_gap, [0.0], subplot=i,
                        color=:black, lw=0.5, ls=:dash, label="")
        end
        _savefig_safe(p_gap, joinpath(FIGURES_DIR, "irf_output_gap_$(tag).png"))
    end

    @printf "\n--- Figures complete: %d saved ---\n" _N_FIGS_SAVED[]
    @printf "  Figures → %s\n" FIGURES_DIR
    @printf "  Tables  → %s\n" TABLES_DIR
    _N_FIGS_SAVED[] = 0   # reset for a possible second call in the same session
end
