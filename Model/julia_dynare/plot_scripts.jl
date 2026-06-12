"""
plot_scripts.jl
===============
Julia translations of all MATLAB plotting scripts called from main_SOE_gap.m:

  figs_SOE_gap.m              → plot_output_gaps()
  plot_manufacturing_shock.m  → plot_manufacturing_shock()
  plot_figure7_manufacturing_shock.m → plot_figure7()
  plot_shock_effects.m        → plot_shock_effects()
  steady_state_table.m        → print_steady_state_table()

All PDF/PNG figures saved to FIGURES_DIR.
All CSV tables saved to TABLES_DIR.

Requires Plots.jl (] add Plots). Text output works without Plots.jl.
"""

using CSV, DataFrames, Printf, Statistics

# =========================================================================== #
#  HELPERS                                                                     #
# =========================================================================== #

const _PLOTS_LOADED = Ref(false)
function ensure_plots()
    _PLOTS_LOADED[] && return true
    try
        @eval Main using Plots
        @eval Main gr(dpi=200)
        _PLOTS_LOADED[] = true
        return true
    catch
        return false
    end
end

function get_irf2(df::DataFrame, varname::AbstractString, shock::AbstractString;
                  nT::Int=40)
    sub = filter(r -> String(r.variable) == String(varname) &&
                      String(r.shock) == String(shock), df)
    isempty(sub) && return zeros(nT)
    s = sort(sub, :period)
    out = zeros(nT)
    n = min(nT, nrow(s))
    out[1:n] .= Float64.(s.value[1:n])
    return out
end

function savefig2(p, path::String)
    try
        @eval Main savefig($p, $path)
        @printf "  → %s\n" basename(path)
    catch e
        @printf "  (could not save %s: %s)\n" basename(path) string(e)
    end
end

function hline0!(p; subplot=nothing, kw...)
    kwargs = isnothing(subplot) ? kw : (subplot=subplot, kw...)
    Plots.hline!(p, [0.0]; color=:black, lw=0.6, ls=:dash, label="", kwargs...)
end


# =========================================================================== #
#  1.  plot_output_gaps()  ←  figs_SOE_gap.m                                 #
#                                                                              #
#  Plots output gap (Ygap_i) time series for all 12 sectors.                  #
#  Equivalent to: aggregate_output_gap.png + sectoral_output_gaps.png         #
# =========================================================================== #

function plot_output_gaps(df_irf, main_shock, names_vec, FIGURES_DIR, tag; nT=40)
    periods = 1:nT
    pok = ensure_plots()
    sector_names = [string(n)[1:min(14,length(string(n)))] for n in names_vec]

    # ---- Aggregate gap (Ygap = NK output gap) ----------------------------
    @printf "\n  [figs_SOE_gap] Aggregate output gap\n"
    ygap = get_irf2(df_irf, "Ygap", main_shock; nT=nT)
    gdpgap = get_irf2(df_irf, "GDPgap", main_shock; nT=nT)
    @printf "  Ygap(peak)=%.4f  GDPgap(peak)=%.4f\n" maximum(abs.(ygap)) maximum(abs.(gdpgap))

    if pok
        p = Plots.plot(periods, [ygap gdpgap],
            label=["Output gap (Ygap)" "GDP gap (GDPgap)"],
            color=[:steelblue :firebrick], lw=2,
            xlabel="Quarters", ylabel="% dev. from SS",
            title="Aggregate Output Gap — $main_shock", legend=:topright)
        Plots.hline!(p, [0.0], color=:black, lw=0.8, ls=:dash, label="")
        savefig2(p, joinpath(FIGURES_DIR, "aggregate_output_gap_$(tag).png"))
    end

    # ---- Sectoral output gaps (4×3 panel) --------------------------------
    @printf "  [figs_SOE_gap] Sectoral output gaps\n"
    if pok
        p2 = Plots.plot(layout=(4,3), size=(1200,900),
                        plot_title="Sectoral Output Gaps — $main_shock",
                        titlefontsize=8)
        for i in 1:12
            irf = get_irf2(df_irf, "Ygap_$(i)", main_shock; nT=nT)
            Plots.plot!(p2, periods, irf, subplot=i,
                       title=sector_names[i], label="",
                       color=:steelblue, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p2, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p2, joinpath(FIGURES_DIR, "sectoral_output_gaps_$(tag).png"))
    end
end


# =========================================================================== #
#  2.  plot_manufacturing_shock()  ←  plot_manufacturing_shock.m              #
#                                                                              #
#  Four panels specific to the manufacturing TFP shock (epsA_3):             #
#    a) Key macro variables (3×3)                                             #
#    b) Aggregate gaps (1×3)                                                  #
#    c) Goods vs services (2×3)                                               #
#    d) Sectoral output gaps (4×3) — manufacturing in red                    #
# =========================================================================== #

function plot_manufacturing_shock(df_irf, names_vec, FIGURES_DIR, tag; nT=40)
    shock = "epsA_3"   # manufacturing TFP
    periods = 1:nT
    pok = ensure_plots()
    sector_names = [string(n)[1:min(14,length(string(n)))] for n in names_vec]

    @printf "\n  [plot_manufacturing_shock] Shock: %s\n" shock

    # ---- (a) Key macro variables 3×3 ------------------------------------
    macro_vars = [("GDP","GDP"),("pi","Inflation π"),("Q","Real XR (Q)"),
                  ("C","Consumption"),("N","Labor (N)"),("w","Wage (w)"),
                  ("TB","Trade balance"),("X","Exports"),("IMP","Imports")]
    @printf "  Panel (a): key macro variables\n"
    if pok
        p_a = Plots.plot(layout=(3,3), size=(1200,900),
                         plot_title="Manufacturing TFP Shock — Key Variables",
                         titlefontsize=8)
        for (k,(vn,lbl)) in enumerate(macro_vars)
            irf = get_irf2(df_irf, vn, shock; nT=nT)
            Plots.plot!(p_a, periods, irf, subplot=k, title=lbl, label="",
                       color=:steelblue, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p_a, [0.0], subplot=k, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p_a, joinpath(FIGURES_DIR, "manufacturing_macro_$(tag).png"))
    end

    # ---- (b) Aggregate output gaps 1×3 ----------------------------------
    @printf "  Panel (b): aggregate gaps\n"
    gap_vars = [("Ygap","Output gap"),("GDPgap","GDP gap"),("Ngap","Labor gap")]
    if pok
        p_b = Plots.plot(layout=(1,3), size=(1200,350),
                         plot_title="Aggregate Gaps — Manufacturing TFP",
                         titlefontsize=9)
        for (k,(vn,lbl)) in enumerate(gap_vars)
            irf = get_irf2(df_irf, vn, shock; nT=nT)
            Plots.plot!(p_b, periods, irf, subplot=k, title=lbl, label="",
                       color=:firebrick, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p_b, [0.0], subplot=k, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p_b, joinpath(FIGURES_DIR, "gap_aggregate_mfg_$(tag).png"))
    end

    # ---- (c) Goods vs services (2×3) ------------------------------------
    @printf "  Panel (c): goods vs services\n"
    gs_vars = [("Y_g","Output: Goods"),("Y_s","Output: Services"),
               ("Ygap_g","Gap: Goods"),("Ygap_s","Gap: Services"),
               ("C_g","Consumption: Goods"),("C_s","Consumption: Services")]
    if pok
        p_c = Plots.plot(layout=(2,3), size=(1200,600),
                         plot_title="Goods vs Services — Manufacturing TFP",
                         titlefontsize=8)
        for (k,(vn,lbl)) in enumerate(gs_vars)
            irf = get_irf2(df_irf, vn, shock; nT=nT)
            clr = contains(lbl,"Goods") ? :steelblue : :darkorange
            Plots.plot!(p_c, periods, irf, subplot=k, title=lbl, label="",
                       color=clr, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p_c, [0.0], subplot=k, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p_c, joinpath(FIGURES_DIR, "gap_goods_services_mfg_$(tag).png"))
    end

    # ---- (d) Sectoral output gaps — manufacturing highlighted in red -----
    @printf "  Panel (d): sectoral output gaps\n"
    if pok
        p_d = Plots.plot(layout=(4,3), size=(1200,900),
                         plot_title="Sectoral Output Gaps — Manufacturing TFP",
                         titlefontsize=8)
        for i in 1:12
            irf = get_irf2(df_irf, "Ygap_$(i)", shock; nT=nT)
            clr = (i == 3) ? :firebrick : :steelblue   # sector 3 = Manufacturing
            lw  = (i == 3) ? 2.5 : 1.5
            Plots.plot!(p_d, periods, irf, subplot=i,
                       title=sector_names[i], label="",
                       color=clr, lw=lw, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p_d, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p_d, joinpath(FIGURES_DIR, "gap_sectoral_mfg_$(tag).png"))
    end
end


# =========================================================================== #
#  3.  plot_figure7()  ←  plot_figure7_manufacturing_shock.m                  #
#                                                                              #
#  Figure 7 from the paper: three sub-figures (7a, 7b, 7c)                   #
#    7a: Manufacturing + open economy variables (2×4)                         #
#    7b: Aggregate macro effects (3×4)                                        #
#    7c: Sectoral imports — manufacturing highlighted (3×4)                   #
# =========================================================================== #

function plot_figure7(df_irf, names_vec, FIGURES_DIR, tag; nT=40)
    shock = "epsA_3"
    periods = 1:nT
    pok = ensure_plots()
    sector_names = [string(n)[1:min(14,length(string(n)))] for n in names_vec]

    @printf "\n  [plot_figure7] Figure 7: Manufacturing shock deep-dive\n"

    # ---- Figure 7a: Manufacturing sector + OE variables (2×4) -----------
    f7a_vars = [
        ("Y_3","Manufacturing Y"),("PH_3","Manufacturing PH"),
        ("L_3","Manufacturing L"),("M_3","Manufacturing M"),
        ("Q","Real XR (Q)"),("X","Exports"),
        ("IMP","Imports"),("TB","Trade balance"),
    ]
    if pok
        p7a = Plots.plot(layout=(2,4), size=(1400,700),
                         plot_title="Fig 7a: Manufacturing & Open Economy",
                         titlefontsize=8)
        for (k,(vn,lbl)) in enumerate(f7a_vars)
            irf = get_irf2(df_irf, vn, shock; nT=nT)
            clr = k <= 4 ? :firebrick : :darkgreen
            Plots.plot!(p7a, periods, irf, subplot=k, title=lbl, label="",
                       color=clr, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p7a, [0.0], subplot=k, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p7a, joinpath(FIGURES_DIR, "figure7a_$(tag).png"))
    end

    # ---- Figure 7b: Aggregate macro effects (3×4) -----------------------
    f7b_vars = [
        ("GDP","GDP"),("C","Consumption"),("N","Labor (N)"),("w","Wage"),
        ("pi","Inflation"),("r","Interest rate"),("Q","Real XR"),("TB","Trade bal."),
        ("Y_g","Goods output"),("Y_s","Services out."),("Ygap","Output gap"),("GDPgap","GDP gap"),
    ]
    if pok
        p7b = Plots.plot(layout=(3,4), size=(1400,900),
                         plot_title="Fig 7b: Aggregate Effects — Manufacturing TFP",
                         titlefontsize=8)
        for (k,(vn,lbl)) in enumerate(f7b_vars)
            irf = get_irf2(df_irf, vn, shock; nT=nT)
            Plots.plot!(p7b, periods, irf, subplot=k, title=lbl, label="",
                       color=:steelblue, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p7b, [0.0], subplot=k, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p7b, joinpath(FIGURES_DIR, "figure7b_$(tag).png"))
    end

    # ---- Figure 7c: Sectoral imports (4×3) — manufacturing in red -------
    if pok
        p7c = Plots.plot(layout=(4,3), size=(1200,900),
                         plot_title="Fig 7c: Sectoral Imports — Manufacturing TFP",
                         titlefontsize=8)
        for i in 1:12
            irf = get_irf2(df_irf, "M_$(i)", shock; nT=nT)
            clr = (i == 3) ? :firebrick : :steelblue
            lw  = (i == 3) ? 2.5 : 1.5
            Plots.plot!(p7c, periods, irf, subplot=i,
                       title=sector_names[i], label="",
                       color=clr, lw=lw, xlabel="Q", ylabel="% dev.")
            Plots.hline!(p7c, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
        end
        savefig2(p7c, joinpath(FIGURES_DIR, "figure7c_sectoral_imports_$(tag).png"))
    end
end


# =========================================================================== #
#  4.  plot_shock_effects()  ←  plot_shock_effects.m                          #
#                                                                              #
#  Comparison plot for preference (eps_om/eps_xi) and monetary (eps_i) shocks#
#  4×4 grid (15 panels) with two shocks overlaid per panel.                  #
# =========================================================================== #

function plot_shock_effects(df_irf, EXERCISE, FIGURES_DIR, tag; nT=40)
    periods = 1:nT
    pok = ensure_plots()

    # Select shock pair based on exercise
    if EXERCISE == 1
        shock_a = "eps_om"; shock_b = "eps_xi"
        label_a = "Goods-services pref."; label_b = "Preference (xi)"
    elseif EXERCISE == 3
        shock_a = "eps_i"; shock_b = "eps_xi"
        label_a = "Monetary policy"; label_b = "Preference (xi)"
    else   # Baseline: compare TFP vs preference
        shock_a = "epsA_3"; shock_b = "eps_pvstar"
        label_a = "Manufacturing TFP"; label_b = "Import price"
    end

    @printf "\n  [plot_shock_effects] %s vs %s\n" shock_a shock_b

    vars_4x4 = [
        ("GDP","GDP"),("pi","Inflation"),("Q","Real XR"),("TB","Trade bal."),
        ("C","Consumption"),("N","Labor"),("w","Wage"),("r","Int. rate"),
        ("Y_g","Goods output"),("Y_s","Services out."),("X","Exports"),("IMP","Imports"),
        ("Ygap","Output gap"),("GDPgap","GDP gap"),("Ngap","Labor gap"),
    ]

    if pok
        p = Plots.plot(layout=(4,4), size=(1400,1100),
                       plot_title="Shock Comparison: $label_a vs $label_b",
                       titlefontsize=7, legend=false)
        for (k,(vn,lbl)) in enumerate(vars_4x4)
            k > 15 && break
            irf_a = get_irf2(df_irf, vn, shock_a; nT=nT)
            irf_b = get_irf2(df_irf, vn, shock_b; nT=nT)
            Plots.plot!(p, periods, irf_a, subplot=k, title=lbl, label=label_a,
                       color=:steelblue, lw=2, xlabel="Q", ylabel="% dev.")
            Plots.plot!(p, periods, irf_b, subplot=k, label=label_b,
                       color=:firebrick, lw=2, ls=:dash)
            Plots.hline!(p, [0.0], subplot=k, color=:black, lw=0.5, ls=:dot, label="")
        end
        savefig2(p, joinpath(FIGURES_DIR, "shock_comparison_$(tag).png"))
    end
end


# =========================================================================== #
#  5.  print_steady_state_table()  ←  steady_state_table.m                   #
#                                                                              #
#  Prints and saves a formatted steady-state summary table.                   #
# =========================================================================== #

function print_steady_state_table(ss_results, sec_results, names_vec,
                                  TABLES_DIR, tag; ombar=0.57)
    @printf "\n============================================================\n"
    @printf "  STEADY-STATE SUMMARY\n"
    @printf "============================================================\n\n"

    # ---- Aggregate variables --------------------------------------------
    GDP  = ss_results.GDP_ss
    TB   = ss_results.TB_ss
    Q    = ss_results.Q_ss
    C    = ss_results.C_ss
    N    = ss_results.N_ss
    Bstar = ss_results.Bstar_ss
    w    = ss_results.w_ss

    @printf "  %-30s  %12s\n" "Aggregate Variable" "Value"
    @printf "  %s\n" repeat("-", 46)
    @printf "  %-30s  %12.4f\n" "GDP"                    GDP
    @printf "  %-30s  %12.4f\n" "Trade balance (TB)"     TB
    @printf "  %-30s  %12.4f\n" "TB/GDP (%%)"             100*TB/max(GDP,1e-10)
    @printf "  %-30s  %12.4f\n" "Real exchange rate (Q)" Q
    @printf "  %-30s  %12.4f\n" "Consumption (C)"        C
    @printf "  %-30s  %12.4f\n" "C/GDP (%%)"              100*C/max(GDP,1e-10)
    @printf "  %-30s  %12.4f\n" "Labor (N)"              N
    @printf "  %-30s  %12.4f\n" "Wage (w)"               w
    @printf "  %-30s  %12.4f\n" "Foreign debt (Bstar)"  Bstar
    @printf "  %-30s  %12.4f\n" "Q*Bstar/GDP (%%)"        100*Q*abs(Bstar)/max(GDP,1e-10)
    @printf "  %-30s  %12.4f\n" "Goods share (omG)"     ombar
    @printf "\n"

    # ---- Sectoral table -------------------------------------------------
    Yi_ss = Float64.(sec_results.Yi_ss)
    Y_tot = sum(Yi_ss)
    @printf "  %-6s  %-30s  %8s  %8s  %8s  %8s\n" "Sector" "Name" "Y_ss" "Y/Y_tot" "L_ss" "PH_ss"
    @printf "  %s\n" repeat("-", 74)
    for i in 1:length(Yi_ss)
        r  = sec_results[i, :]
        nm = string(names_vec[i])[1:min(28, length(string(names_vec[i])))]
        @printf "  %-6d  %-30s  %8.4f  %7.1f%%  %8.4f  %8.4f\n" i nm r.Yi_ss 100*r.Yi_ss/max(Y_tot,1e-10) r.L_ss r.pH_ss
    end
    @printf "\n  Gross output Y = %.4f   GDP = %.4f   Y/GDP = %.2f\n\n" Y_tot GDP Y_tot/max(GDP,1e-10)

    # ---- Save CSV -------------------------------------------------------
    df_ss_table = DataFrame(
        sector    = 1:length(Yi_ss),
        name      = [string(n) for n in names_vec],
        Y_ss      = Yi_ss,
        Y_share   = Yi_ss ./ max(Y_tot, 1e-10),
        L_ss      = Float64.(sec_results.L_ss),
        PH_ss     = Float64.(sec_results.pH_ss),
        GDP_ss    = fill(GDP, length(Yi_ss)),
        TB_ss     = fill(TB, length(Yi_ss)),
        Q_ss      = fill(Q,  length(Yi_ss)),
        C_ss      = fill(C,  length(Yi_ss)),
    )
    CSV.write(joinpath(TABLES_DIR, "steady_state_$(tag).csv"), df_ss_table)
    @printf "  Steady-state table saved to: tables/steady_state_%s.csv\n" tag
end


# =========================================================================== #
#  DISPATCH — called from main_SOE_gap.jl                                     #
# =========================================================================== #

"""
    run_all_plots(; df_irf, EXERCISE, names_vec, ss_results, sec_results,
                    FIGURES_DIR, TABLES_DIR, tag, ombar)

Run all plotting scripts appropriate for the given exercise.
Mirrors the if/elseif block in main_SOE_gap.m lines 968–994.
"""
function run_all_plots(;
    df_irf::DataFrame,
    EXERCISE::Int,
    names_vec,
    ss_results::NamedTuple,
    sec_results::DataFrame,
    FIGURES_DIR::String,
    TABLES_DIR::String,
    tag::String,
    ombar::Float64 = 0.57,
    nT::Int = 40,
)
    # Determine main shock
    all_shocks = unique(String.(df_irf.shock))
    main_shock = if EXERCISE == 2
        s = filter(x -> contains(x,"epsA_3"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    elseif EXERCISE == 1
        s = filter(x -> contains(x,"eps_om") || contains(x,"eps_xi"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    elseif EXERCISE == 3
        s = filter(x -> contains(x,"eps_i"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    else
        s = filter(x -> contains(x,"epsA_3"), all_shocks)
        isempty(s) ? first(all_shocks) : first(s)
    end

    @printf "\n%s\n  Generating figures  (exercise %d, shock: %s)\n%s\n" repeat("=",60) EXERCISE main_shock repeat("=",60)

    # 1. Output gaps — always generated
    plot_output_gaps(df_irf, main_shock, names_vec, FIGURES_DIR, tag; nT=nT)

    # 2. Exercise-specific plots
    if EXERCISE == 2
        plot_manufacturing_shock(df_irf, names_vec, FIGURES_DIR, tag; nT=nT)
        plot_figure7(df_irf, names_vec, FIGURES_DIR, tag; nT=nT)
    else
        plot_shock_effects(df_irf, EXERCISE, FIGURES_DIR, tag; nT=nT)
    end

    # 3. Steady-state table — always
    print_steady_state_table(ss_results, sec_results, names_vec,
                              TABLES_DIR, tag; ombar=ombar)

    if ensure_plots()
        @printf "\n  All figures saved to: %s\n" FIGURES_DIR
    else
        @printf "\n  (Plots.jl not installed — text output only)\n"
        @printf "  Install: using Pkg; Pkg.add(\"Plots\")\n"
    end
end
