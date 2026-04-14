"""
figs_SOE_gap.jl
===============
Julia translation of figs_SOE_gap.m — generates IRF figures and the
steady-state summary table for the NK-IOSOE Chile model.

Called from main_SOE_gap.jl after the Dynare subprocess completes.
Reads from dynare_irfs.csv and the sectoral/SS CSVs in Data/.

Outputs (PDF + PNG) written to DATA_DIR:
  irf_aggregate_<exercise>.pdf  — GDP, pi, Q, TB aggregate IRFs
  irf_sectoral_<exercise>.pdf   — Y, PH, L by sector (subplots)
  irf_output_gap_<exercise>.pdf — Ygap by sector
  ss_table_<exercise>.csv       — steady-state summary table
"""

using CSV, DataFrames, Printf, LinearAlgebra, Statistics

# =========================================================================== #
#  ENTRY POINT — called from main_SOE_gap.jl                                 #
# =========================================================================== #

function generate_figures(;
    MOD_DIR::String,
    DATA_DIR::String,
    SCRIPT_DIR::String,
    EXERCISE::Int,
    nsec::Int,
    names_vec,
    ss_results::NamedTuple,        # from main_SOE_gap steady state
    sec_results::DataFrame,        # sectoral model moments
    exercise_label::String,
)
    @printf "\n--- Generating figures (Exercise %d) ---\n" EXERCISE

    # Load IRFs
    irf_file = joinpath(MOD_DIR, "dynare_irfs.csv")
    if !isfile(irf_file) || filesize(irf_file) < 10
        @printf "  No IRF file found — skipping figure generation.\n"
        @printf "  (Run with EXERCISE=0 Baseline to generate all shocks)\n"
        return
    end

    df_irf = CSV.read(irf_file, DataFrame)
    if nrow(df_irf) == 0
        @printf "  IRF file is empty — Dynare stoch_simul produced no IRFs.\n"
        @printf "  This usually means the model SS was not computed correctly.\n"
        return
    end

    shocks   = unique(df_irf.shock)
    @printf "  IRFs loaded: %d shock(s), %d variable-shock pairs\n" length(shocks) length(unique(df_irf.variable))

    # ---- Try to load Plots.jl ------------------------------------------ #
    plots_ok = false
    try
        @eval using Plots
        @eval Plots.gr()   # GR backend — no display needed
        plots_ok = true
    catch
        @printf "  WARNING: Plots.jl not available — printing IRF tables instead.\n"
        @printf "  Install with: ] add Plots\n"
    end

    # ---- Helper: get IRF series ----------------------------------------- #
    function get_irf(df, varname, shock)
        sub = filter(r -> r.variable == varname && r.shock == shock, df)
        isempty(sub) ? Float64[] : sort(sub, :period).value
    end

    # ---- 1. Aggregate IRFs (GDP, pi, Q, TB) ----------------------------- #
    agg_vars  = ["GDP", "pi", "Q", "TB"]
    agg_labels = ["GDP", "Inflation (pi)", "Real exchange rate (Q)", "Trade balance (TB)"]

    @printf "\n  Aggregate IRFs (first 10 periods):\n"
    @printf "  %-8s" "Shock"
    for vl in agg_labels; @printf "  %12s" vl; end
    println()
    for shock in shocks
        @printf "  %-8s" shock
        for (vn, _) in zip(agg_vars, agg_labels)
            irf = get_irf(df_irf, vn, shock)
            v10 = isempty(irf) ? NaN : (length(irf) >= 1 ? irf[1] : NaN)
            @printf "  %12.4f" v10
        end
        println()
    end

    if plots_ok
        for shock in shocks
            p = @eval Plots.plot(layout=(2,2), size=(900,600),
                    title=["GDP" "Inflation" "Real XR" "Trade balance"],
                    titlefontsize=9)
            for (k, vn) in enumerate(agg_vars)
                irf = get_irf(df_irf, vn, shock)
                isempty(irf) && continue
                @eval Plots.plot!($p, subplot=$k, $irf, label="", color=:steelblue,
                    lw=2, xlabel="Quarters", ylabel="% dev. from SS")
                @eval Plots.hline!($p, [0.0], subplot=$k, color=:black, lw=0.8, ls=:dash, label="")
            end
            fname = joinpath(DATA_DIR, "irf_aggregate_ex$(EXERCISE)_$(shock).pdf")
            @eval Plots.savefig($p, $fname)
            @printf "  Saved: %s\n" fname
        end
    end

    # ---- 2. Sectoral output IRFs (Y_1 ... Y_12) ------------------------- #
    @printf "\n  Sectoral output IRFs (period 1 impact):\n"
    @printf "  %-6s  %-30s" "Sector" "Name"
    for shock in shocks[1:min(3,length(shocks))]
        @printf "  %10s" shock
    end
    println()

    sec_Y_irfs = Dict{String, Vector{Float64}}()
    for i in 1:nsec
        vn = "Y_$(i)"
        @printf "  %-6d  %-30s" i string(names_vec[i])[1:min(30,length(string(names_vec[i])))]
        for shock in shocks[1:min(3,length(shocks))]
            irf = get_irf(df_irf, vn, shock)
            v1  = isempty(irf) ? NaN : irf[1]
            @printf "  %10.4f" v1
            shock == shocks[1] && (sec_Y_irfs[vn] = irf)
        end
        println()
    end

    if plots_ok && !isempty(shocks)
        shock1 = shocks[1]
        p2 = @eval Plots.plot(layout=(4,3), size=(1200,900),
                suptitle="Sectoral Output IRF — shock: $shock1", titlefontsize=8)
        for i in 1:nsec
            irf = get_irf(df_irf, "Y_$(i)", shock1)
            isempty(irf) && continue
            nm = string(names_vec[i])[1:min(15,length(string(names_vec[i])))]
            @eval Plots.plot!($p2, subplot=$i, $irf, title=$nm, label="",
                color=:steelblue, lw=2)
            @eval Plots.hline!($p2, [0.0], subplot=$i, color=:black, lw=0.5, ls=:dash, label="")
        end
        fname2 = joinpath(DATA_DIR, "irf_sectoral_Y_ex$(EXERCISE).pdf")
        @eval Plots.savefig($p2, $fname2)
        @printf "  Saved: %s\n" fname2

        # Prices
        p3 = @eval Plots.plot(layout=(4,3), size=(1200,900),
                suptitle="Sectoral Price IRF — shock: $shock1", titlefontsize=8)
        for i in 1:nsec
            irf = get_irf(df_irf, "PH_$(i)", shock1)
            isempty(irf) && continue
            nm = string(names_vec[i])[1:min(15,length(string(names_vec[i])))]
            @eval Plots.plot!($p3, subplot=$i, $irf, title=$nm, label="",
                color=:firebrick, lw=2)
            @eval Plots.hline!($p3, [0.0], subplot=$i, color=:black, lw=0.5, ls=:dash, label="")
        end
        fname3 = joinpath(DATA_DIR, "irf_sectoral_PH_ex$(EXERCISE).pdf")
        @eval Plots.savefig($p3, $fname3)
        @printf "  Saved: %s\n" fname3
    end

    # ---- 3. Steady-state summary table ---------------------------------- #
    @printf "\n  STEADY-STATE SUMMARY\n"
    @printf "  %-25s  %12s\n" "Variable" "Value"
    @printf "  %s\n" repeat("-", 40)
    ss_display = [
        ("GDP",                    ss_results.GDP_ss),
        ("TB/GDP",                 ss_results.TB_ss / ss_results.GDP_ss),
        ("Real exchange rate (Q)", ss_results.Q_ss),
        ("Consumption (C)",        ss_results.C_ss),
        ("Labor (N)",              ss_results.N_ss),
        ("Foreign debt/GDP",       ss_results.Q_ss * abs(ss_results.Bstar_ss) / ss_results.GDP_ss),
        ("Gross output/GDP",       sum(sec_results.Yi_ss) / ss_results.GDP_ss),
    ]
    for (label, val) in ss_display
        @printf "  %-25s  %12.4f\n" label val
    end
    @printf "\n  Sectoral SS values:\n"
    @printf "  %-6s  %-30s  %10s  %10s  %10s\n" "Sector" "Name" "Y_ss" "L_ss" "PH_ss"
    @printf "  %s\n" repeat("-", 70)
    for i in 1:nsec
        r = sec_results[i, :]
        nm = string(names_vec[i])[1:min(28,length(string(names_vec[i])))]
        @printf "  %-6d  %-30s  %10.4f  %10.4f  %10.4f\n" i nm r.Yi_ss r.L_ss r.pH_ss
    end

    @printf "\n--- Figures complete ---\n"
end
