"""
shock_plots_common.jl
=====================
Shared figure / table / decomposition routines for the NK-IOSOE shock-analysis
scripts (oil, agriculture, manufacturing, mining).

This file is `include`d by each `*_shock_analysis.jl` script AFTER it has solved
the model and assembled the IRFs.  It centralises every figure and table so that
all shocks produce the SAME publication set, differing only in titles, filenames
and the shock-specific decomposition inputs.

Design
------
Each caller builds a NamedTuple `ctx` holding the IRF arrays, steady-state
scalars, label strings and a `fignames::Dict{Symbol,String}` filename map, then
calls `generate_shock_outputs(ctx)`.

Helper builders exposed for the callers:
  * `compute_inflation_aggregates(...)`  — sectoral & group home-price inflation
  * `leontief_decomp(...)`               — direct / network / total MC decomposition
  * `ge_mc_components(...)`              — 5-way marginal-cost decomposition (per horizon)
  * `standard_fignames(tag)`             — `<base>_<tag>_shock.pdf` naming scheme
  * `oil_fignames()`                     — legacy oil filenames (paper-compatible)

The plotting code resolves `Plots` / `groupedbar` from `Main` at call time, so
this file can be included even when Plots.jl is unavailable (figures are then
skipped by the caller via `ctx.HAS_PLOTS`).
"""

# Stdlib imports needed at load time (@printf/@sprintf expand at definition).
# `Plots`, `groupedbar`, `CSV`, `DataFrame(s)` are resolved from `Main` lazily at
# call time, so they are intentionally NOT imported here.
using Printf
using LinearAlgebra
using Logging

# =========================================================================== #
#  FILENAME SCHEMES                                                            #
# =========================================================================== #

"Standardised `<base>_<tag>_shock.pdf` filenames used by the TFP shock scripts."
function standard_fignames(tag::AbstractString)
    Dict{Symbol,String}(
        :agg              => "irf_aggregate_$(tag)_shock.pdf",
        :sec_Y            => "irf_sectoral_Y_$(tag)_shock.pdf",
        :sec_PH           => "irf_sectoral_PH_$(tag)_shock.pdf",
        :sec_inflation    => "irf_sectoral_inflation_$(tag)_shock.pdf",
        :agg_inflation    => "irf_aggregate_inflation_$(tag)_shock.pdf",
        :gs_inflation     => "irf_goods_vs_services_inflation_$(tag)_shock.pdf",
        :sec_MC           => "irf_sectoral_MC_$(tag)_shock.pdf",
        :labor_agg        => "irf_labor_aggregate_$(tag)_shock.pdf",
        :sec_L            => "irf_sectoral_L_$(tag)_shock.pdf",
        :gdpgap_agg       => "irf_gdpgap_aggregate_$(tag)_shock.pdf",
        :sec_Ygap         => "irf_sectoral_Ygap_$(tag)_shock.pdf",
        :exposure         => "$(tag)_intensity_exposure.pdf",
        :decomp_mc_impact => "decomposition_mc_inflation_$(tag)_baseline.pdf",
        :decomp_mc_6m     => "decomposition_mc_inflation_6m_$(tag)_baseline.pdf",
        :decomp_mc_12m    => "decomposition_mc_inflation_12m_$(tag)_baseline.pdf",
        # 3-way grouped version (direct / indirect / others)
        :decomp_mc3_impact => "decomposition_mc3_inflation_$(tag)_baseline.pdf",
        :decomp_mc3_6m     => "decomposition_mc3_inflation_6m_$(tag)_baseline.pdf",
        :decomp_mc3_12m    => "decomposition_mc3_inflation_12m_$(tag)_baseline.pdf",
        :decomp_output    => "decomposition_output_$(tag)_baseline.pdf",
    )
end

"Legacy oil filenames (kept verbatim so existing \\includegraphics in the paper resolve)."
function oil_fignames()
    Dict{Symbol,String}(
        :agg              => "irf_aggregate_oil_shock.pdf",
        :sec_Y            => "irf_sectoral_Y_oil_shock.pdf",
        :sec_PH           => "irf_sectoral_PH_oil_shock.pdf",
        :sec_inflation    => "irf_sectoral_inflation_oil.pdf",
        :agg_inflation    => "irf_aggregate_inflation_oil.pdf",
        :gs_inflation     => "irf_goods_vs_services_inflation_oil.pdf",
        :sec_MC           => "irf_sectoral_MC_oil_shock.pdf",
        :labor_agg        => "irf_labor_aggregate_oil_shock.pdf",
        :sec_L            => "irf_sectoral_L_oil_shock.pdf",
        :gdpgap_agg       => "irf_gdpgap_aggregate_oil_shock.pdf",
        :sec_Ygap         => "irf_sectoral_Ygap_oil_shock.pdf",
        :exposure         => "oil_intensity_exposure.pdf",
        :decomp_mc_impact => "decomposition_mc_inflation_oil_baseline.pdf",
        :decomp_mc_6m     => "decomposition_mc_inflation_6m_oil_baseline.pdf",
        :decomp_mc_12m    => "decomposition_mc_inflation_12m_oil_baseline.pdf",
        # 3-way grouped version (direct / indirect / others)
        :decomp_mc3_impact => "decomposition_mc3_inflation_oil_baseline.pdf",
        :decomp_mc3_6m     => "decomposition_mc3_inflation_6m_oil_baseline.pdf",
        :decomp_mc3_12m    => "decomposition_mc3_inflation_12m_oil_baseline.pdf",
        :decomp_output    => "decomposition_output_oil_baseline.pdf",
    )
end

# Standard sector display names in Spanish (consistent across every figure/table).
const SHOCK_BAR_NAMES = [
    "Agric. y Pesca", "Minería", "Manufactura", "Electricidad",
    "Construcción", "Comercio y Hot.", "Transp. y Com.", "Finanzas",
    "Inmobiliario", "Serv. Empres.", "Serv. Pers.", "Adm. Pública",
]
const SHOCK_PANEL_NAMES = [
    "Agricultura y Pesca", "Minería", "Manufactura", "Electricidad y Agua",
    "Construcción", "Comercio y Hoteles", "Transporte y Com.", "Finanzas",
    "Inmobiliario", "Serv. Empresariales", "Serv. Personales", "Adm. Pública",
]
const SHOCK_TAB_NAMES = [
    "Agricultura y Pesca", "Minería", "Manufactura", "Electricidad",
    "Construcción", "Comercio y Hoteles", "Transporte y Com.", "Finanzas",
    "Inmobiliario", "Serv.\\ Empresariales", "Serv.\\ Personales", "Adm.\\ Pública",
]


# =========================================================================== #
#  COMPUTE HELPERS (shock-agnostic)                                            #
# =========================================================================== #

"""
    compute_inflation_aggregates(ph_irf_mat, pi_irf, C_gi_ss, C_si_ss, nsec, n_irf)

Sectoral nominal home-price inflation (annualised, pp) and consumption-weighted
goods / services / aggregate inflation, plus impact / 6-month / 12-month
cumulative home-price inflation.  Mirrors the oil-shock construction.
Returns a NamedTuple.
"""
function compute_inflation_aggregates(ph_irf_mat, pi_irf, C_gi_ss, C_si_ss, nsec, n_irf)
    pi_sec_mat = zeros(nsec, n_irf)
    for i in 1:nsec
        pi_sec_mat[i, 1] = (ph_irf_mat[i, 1] + pi_irf[1]) * 4
        for h in 2:n_irf
            pi_sec_mat[i, h] = (ph_irf_mat[i, h] - ph_irf_mat[i, h-1] + pi_irf[h]) * 4
        end
    end
    cons_ss_all = C_gi_ss .+ C_si_ss
    w_agg_c = cons_ss_all ./ sum(cons_ss_all)
    w_g_c   = C_gi_ss ./ sum(C_gi_ss)
    w_s_c   = C_si_ss ./ sum(C_si_ss)
    pi_agg_irf   = [sum(w_agg_c .* pi_sec_mat[:, h]) for h in 1:n_irf]
    pi_goods_irf = [sum(w_g_c   .* pi_sec_mat[:, h]) for h in 1:n_irf]
    pi_serv_irf  = [sum(w_s_c   .* pi_sec_mat[:, h]) for h in 1:n_irf]
    # Cumulative home-price inflation (level deviations) at h = 1, 2, 4
    ph_impact = ph_irf_mat[:, 1]
    infl_irf_impact = ph_impact .+ pi_irf[1]
    infl_6m  = ph_irf_mat[:, 2] .+ sum(pi_irf[1:2])
    infl_12m = ph_irf_mat[:, 4] .+ sum(pi_irf[1:4])
    return (pi_sec_mat=pi_sec_mat, pi_agg_irf=pi_agg_irf,
            pi_goods_irf=pi_goods_irf, pi_serv_irf=pi_serv_irf,
            infl_irf_impact=infl_irf_impact, infl_6m=infl_6m, infl_12m=infl_12m)
end

"""
    leontief_decomp(direct_mc, modalpha, modbeta, Yi_ss)

Split a vector of *direct* marginal-cost impulses into the direct component and
the network amplification propagated through the IO matrix via the (partial-
equilibrium) Leontief inverse `(I - diag(α_m) Γ)^{-1}`.  Works for any shock:
the caller supplies the `direct_mc` vector (oil: per-sector oil cost push; TFP:
the shocked sector's own cost push, zero elsewhere).
Returns a NamedTuple with per-sector and output-weighted aggregates.
"""
function leontief_decomp(direct_mc, modalpha, modbeta, Yi_ss)
    nsec = length(direct_mc)
    Leontief = inv(I(nsec) - Diagonal(modalpha) * modbeta)
    total_mc = Leontief * direct_mc
    network_mc = total_mc .- direct_mc
    amp_ratio = total_mc ./ max.(direct_mc, 1e-10)
    yi_w = Yi_ss ./ sum(Yi_ss)
    agg_direct  = sum(yi_w .* direct_mc)
    agg_total   = sum(yi_w .* total_mc)
    agg_network = agg_total - agg_direct
    agg_amp     = agg_total / max(agg_direct, 1e-10)
    return (direct_mc=direct_mc, network_mc=network_mc, total_mc=total_mc,
            amp_ratio=amp_ratio, agg_direct=agg_direct, agg_network=agg_network,
            agg_total=agg_total, agg_amp=agg_amp)
end

"""
    sectoral_va_irf(get_irf, nsec, n_irf, Yi_ss, M_ss, Vi_ss, pH_ss, PMi_ss, PV_ss)

Real (double-deflated, base-SS-price) sectoral VALUE-ADDED IRF, in % deviation
from steady state.  Value added of sector `i` is gross output minus intermediate
materials and imported inputs:
    VA_i = P^H_i Y_i − P^M_i M_i − P^V V_i,
all valued at steady-state prices so the response is a pure quantity (real) move.
This is the GDP-consistent object: a negative TFP shock raises intermediate use
per unit of output, so gross output `Y_i` can rise mechanically while value added
(and hence GDP) falls.  Reporting value added removes that double-counting.
Returns an `nsec × n_irf` matrix of percent deviations.
"""
function sectoral_va_irf(get_irf, nsec, n_irf, Yi_ss, M_ss, Vi_ss, pH_ss, PMi_ss, PV_ss)
    va = zeros(nsec, n_irf)
    for i in 1:nsec
        yi = get_irf("Y_$(i)"); mi = get_irf("M_$(i)"); vi = get_irf("V_$(i)")
        VA_ss_i = pH_ss[i]*Yi_ss[i] - PMi_ss[i]*M_ss[i] - PV_ss*Vi_ss[i]
        for t in 1:n_irf
            dVA = pH_ss[i]*(yi[t]/100*Yi_ss[i]) -
                  PMi_ss[i]*(mi[t]/100*M_ss[i]) -
                  PV_ss*(vi[t]/100*Vi_ss[i])
            va[i, t] = abs(VA_ss_i) > 1e-12 ? 100.0*dVA/VA_ss_i : 0.0
        end
    end
    return va
end

"""
    ge_mc_components(kind, h, ph_irf_mat, mc_irf_mat, w_irf, modalphaV, modalpha,
                     alpha_L_vec, modbeta; po_irf, pv_irf, modalphaOil, a_irf_mat)

5-way decomposition of the sectoral marginal-cost IRF at horizon `h` into
intuitive cost-push channels.  Two `kind`s:

  :price (oil)  → [direct oil, network via prices, non-oil imports, labor, residual]
  :tfp          → [own TFP, network via prices, imports, labor, residual]

Returns `(comps = (c1,c2,c3,c4,c5), mc = mc_h)` plus colour/label vectors via the
companion `ge_component_style(kind)`.
"""
function ge_mc_components(kind::Symbol, h::Int, ph_irf_mat, mc_irf_mat, w_irf,
                          modalphaV, modalpha, alpha_L_vec, modbeta;
                          po_irf=nothing, pv_irf=nothing, modalphaOil=nothing,
                          a_irf_mat=nothing)
    net_h = modalpha .* (modbeta * ph_irf_mat[:, h])
    lab_h = alpha_L_vec .* w_irf[h]
    mc_h  = mc_irf_mat[:, h]
    if kind === :price
        dir_h = modalphaV .* modalphaOil .* po_irf[h]
        imp_h = modalphaV .* (1.0 .- modalphaOil) .* pv_irf[h]
    elseif kind === :tfp
        dir_h = -a_irf_mat[:, h]                      # own productivity push
        imp_h = modalphaV .* pv_irf[h]                # composite imports
    else
        error("ge_mc_components: unknown kind $kind")
    end
    res_h = mc_h .- dir_h .- net_h .- imp_h .- lab_h
    return (comps=(copy(dir_h), copy(net_h), copy(imp_h), copy(lab_h), copy(res_h)),
            mc=copy(mc_h))
end

"Colours and legend labels for the 5-way MC decomposition bars."
function ge_component_style(kind::Symbol)
    colors = [:darkorange, :steelblue, :forestgreen, :crimson, :gray60]
    labels = kind === :price ?
        ["Petróleo directo (pp)", "Red vía precios (pp)", "Import. no-petróleo (pp)",
         "Trabajo (pp)", "Residual (pp)"] :
        ["PTF propia (pp)", "Red vía precios (pp)", "Importaciones (pp)",
         "Trabajo (pp)", "Residual (pp)"]
    return colors, labels
end

"""
    group3_components(comps5)

Collapse the 5-way GE marginal-cost decomposition into three groups for the
coarser "direct / indirect / others" figure:
  * Directo   = own cost push (component 1: own TFP, or direct oil)
  * Indirecto = network propagation via input prices (component 2)
  * Otros     = imports + labour + residual (components 3 + 4 + 5)
`comps5` is the 5-tuple of `nsec`-vectors returned in `ge_mc_components(...).comps`.
Returns a 3-tuple `(direct, indirect, others)`.
"""
function group3_components(comps5)
    direct   = copy(comps5[1])
    indirect = copy(comps5[2])
    others   = comps5[3] .+ comps5[4] .+ comps5[5]
    return (direct, indirect, others)
end

"Colours and legend labels for the 3-way (direct / indirect / others) MC bars."
function group3_style()
    colors = [:darkorange, :steelblue, :gray60]
    labels = ["Directo (pp)", "Indirecto / red (pp)", "Otros (pp)"]
    return colors, labels
end


# =========================================================================== #
#  SAVE HELPERS                                                                #
# =========================================================================== #

# Suppress PlotUtils "No strict ticks" warnings and GR/Qt C-level stderr noise.
function shock_savefig_quiet(p, path)
    try
        redirect_stderr(devnull) do
            with_logger(NullLogger()) do
                Main.Plots.savefig(p, path)
            end
        end
    catch e
        @warn "savefig failed for $path" exception=(e, catch_backtrace())
        Main.Plots.savefig(p, path)   # retry without suppression so errors surface
    end
end

# Save a figure locally and mirror it to the Overleaf folder if available.
function shock_save_fig(p, fname, ctx)
    local_path = joinpath(ctx.FIGURES_DIR, fname)
    shock_savefig_quiet(p, local_path)
    @printf "  Saved: %s\n" fname
    if ctx.overleaf_ok
        ol_path = joinpath(ctx.OVERLEAF_FIG_DIR, fname)
        cp(local_path, ol_path; force=true)
        @printf "  → Overleaf: %s\n" ol_path
    end
end

function shock_sync_table(fname, ctx)
    ctx.overleaf_ok || return
    src = joinpath(ctx.TABLES_DIR, fname)
    dst = joinpath(ctx.OVERLEAF_TAB_DIR, fname)
    isfile(src) || return
    cp(src, dst; force=true)
    @printf "  → Overleaf: %s\n" dst
end

# Solid rectangle for manual stacked-bar plotting.
shock_bar_rect(x, y0, y1, w=0.65) = Main.Plots.Shape(
    [x - w/2, x + w/2, x + w/2, x - w/2],
    [y0,      y0,      y1,      y1])


# =========================================================================== #
#  STACKED GE-DECOMPOSITION FIGURE                                             #
# =========================================================================== #

function make_ge_decomp_fig(comps, comp_colors, comp_labels, infl_h, infl_label,
                            title_str, ylabel_str, bar_names, nsec)
    P = Main.Plots
    comp_mat = hcat(comps...)                 # nsec × 5
    pos_mat  = max.(comp_mat, 0.0)
    neg_mat  = min.(comp_mat, 0.0)
    pos_tops = vec(sum(pos_mat, dims=2))
    neg_tops = vec(sum(neg_mat, dims=2))
    all_vals = vcat(pos_tops, neg_tops, infl_h, [0.0])
    # Additive padding around the true data range, always including zero so the
    # negative (deflationary) part of the response is never clipped.
    dmin = minimum(all_vals); dmax = maximum(all_vals)
    span = max(dmax - dmin, 1e-6)
    ylo = min(dmin, 0.0) - 0.10 * span
    yhi = max(dmax, 0.0) + 0.14 * span
    p = P.plot(
        xticks=(1:nsec, bar_names), xrotation=55,
        ylabel=ylabel_str, title=title_str, titlefontsize=12,
        size=(1400, 660), legend=:topright, ylims=(ylo, yhi),
        bottom_margin=26P.mm, left_margin=14P.mm,
        right_margin=8P.mm, top_margin=3P.mm,
        xlims=(0.3, nsec + 0.7))
    P.hline!(p, [0.0], color=:black, lw=0.6, ls=:dash, label="")
    ncomp = length(comps); bw = 0.65
    for k in 1:ncomp
        pos_bot = k > 1 ? vec(sum(pos_mat[:, 1:k-1], dims=2)) : zeros(nsec)
        neg_bot = k > 1 ? vec(sum(neg_mat[:, 1:k-1], dims=2)) : zeros(nsec)
        lbl_used = false
        for i in 1:nsec
            lbl = lbl_used ? "" : comp_labels[k]
            if pos_mat[i, k] > 1e-10
                s = shock_bar_rect(i, pos_bot[i], pos_bot[i] + pos_mat[i, k], bw)
                P.plot!(p, s, color=comp_colors[k], label=lbl, alpha=0.85,
                        linecolor=:white, linewidth=0.3)
                lbl_used = true; lbl = ""
            end
            if abs(neg_mat[i, k]) > 1e-10
                s = shock_bar_rect(i, neg_bot[i] + neg_mat[i, k], neg_bot[i], bw)
                P.plot!(p, s, color=comp_colors[k], label=lbl, alpha=0.85,
                        linecolor=:white, linewidth=0.3)
                lbl_used = true
            end
        end
        if !lbl_used
            P.plot!(p, [NaN], [NaN], color=comp_colors[k], label=comp_labels[k], lw=4)
        end
    end
    P.scatter!(p, 1:nsec, infl_h, color=:black, markershape=:diamond,
               markersize=7, markerstrokewidth=1, label=infl_label)
    return p
end


# =========================================================================== #
#  MAIN FIGURE SET                                                             #
# =========================================================================== #

function generate_shock_figures(ctx)
    P = Main.Plots
    fn = ctx.fignames
    periods = 1:ctx.nT
    nsec = ctx.nsec
    short_names = SHOCK_PANEL_NAMES   # already short; byte-slicing breaks on accents
    bar_names = SHOCK_BAR_NAMES
    base_lbl = ctx.baseline_label
    base_col = ctx.baseline_color
    base_lw  = ctx.baseline_lw
    epsY_str = "$(round(ctx.epsY_baseline, digits=2))"

    # Padded y-limits that always include zero so the negative part is never clipped.
    function padlims(v)
        dmn = minimum(v); dmx = maximum(v); sp = max(dmx - dmn, 1e-6)
        (min(dmn, 0.0) - 0.10 * sp, max(dmx, 0.0) + 0.14 * sp)
    end

    @printf "\n--- Generando figuras (%s, εY base = %.2f) ---\n" ctx.title_short ctx.epsY_baseline

    # ---- Figura 1: IRFs agregadas (2×5) ---- #
    agg_panels = [
        (ctx.gdp_irf,      "PIB",                    "% desv.",   1.0),
        (ctx.pi_irf,       "Inflación IPC",          "pp anual",  4.0),
        (ctx.q_irf,        "Tipo de Cambio Real",    "% desv.",   1.0),
        (ctx.tb_irf,       "Balanza Comercial",      "pp del PIB",ctx.TB_ss/ctx.GDP_ss),
        (ctx.cg_irf,       "Consumo de Bienes",      "% desv.",   1.0),
        (ctx.r_irf,        "Tasa de Interés Nominal","pp anual",  4.0),
        (ctx.w_irf,        "Salario Real",           "% desv.",   1.0),
        (ctx.pi_goods_irf, "Inflación de Bienes",    "pp anual",  1.0),
        (ctx.pi_serv_irf,  "Inflación de Servicios", "pp anual",  1.0),
        (ctx.cs_irf,       "Consumo de Servicios",   "% desv.",   1.0),
    ]
    p_agg = P.plot(layout=(2,5), size=(2000,700),
        plot_title="$(ctx.title_long) — Respuestas Agregadas",
        titlefontsize=12, margin=5P.mm)
    for (k, (irf_v, ttl, yl, scl)) in enumerate(agg_panels)
        P.plot!(p_agg, periods, irf_v .* scl, subplot=k,
            label=(k==1 ? base_lbl : ""), color=base_col, lw=base_lw,
            title=ttl, ylabel=(k∈[1,6] ? yl : ""), xlabel=(k>5 ? "Trimestres" : ""))
        P.hline!(p_agg, [0.0], subplot=k, color=:black, lw=0.6, ls=:dash, label="")
    end
    shock_save_fig(p_agg, fn[:agg], ctx)

    # ---- Helper para paneles sectoriales 4×3 ---- #
    function sectoral_panel(mat, title_suffix, fname)
        p = P.plot(layout=(4,3), size=(1200,900),
            plot_title="$(ctx.title_short) — $(title_suffix)", titlefontsize=10)
        for i in 1:nsec
            P.plot!(p, periods, mat[i, :], subplot=i, label="",
                color=base_col, lw=base_lw, title=short_names[i], titlefontsize=9)
            P.hline!(p, [0.0], subplot=i, color=:black, lw=0.5, ls=:dash, label="")
        end
        shock_save_fig(p, fname, ctx)
    end

    sectoral_panel(ctx.y_irf_mat,    "Producto Sectorial",                    fn[:sec_Y])
    sectoral_panel(ctx.ph_irf_mat,   "Precios Internos Sectoriales",          fn[:sec_PH])
    sectoral_panel(ctx.pi_sec_mat,   "Inflación de Precios Internos Sectorial", fn[:sec_inflation])
    sectoral_panel(ctx.mc_irf_mat,   "Costo Marginal Sectorial",              fn[:sec_MC])
    sectoral_panel(ctx.l_irf_mat,    "Empleo Sectorial",                      fn[:sec_L])
    sectoral_panel(ctx.ygap_irf_mat, "Brecha de Producto Sectorial",          fn[:sec_Ygap])

    # ---- Inflación de precios internos agregada ---- #
    p_pi_agg = P.plot(size=(1000, 500),
        title="$(ctx.title_short) — Inflación de Precios Internos Agregada",
        titlefontsize=11, xlabel="Trimestres", ylabel="Desv. pp anual del EE",
        legend=:topright, margin=5P.mm)
    P.plot!(p_pi_agg, periods, ctx.pi_agg_irf, label="Agregada", color=base_col, lw=base_lw)
    P.hline!(p_pi_agg, [0.0], color=:black, lw=0.5, ls=:dash, label="")
    shock_save_fig(p_pi_agg, fn[:agg_inflation], ctx)

    # ---- Inflación: bienes vs servicios ---- #
    p_gs = P.plot(size=(800, 500),
        title="$(ctx.title_short) — Inflación de Precios Internos: Bienes vs Servicios",
        titlefontsize=11, xlabel="Trimestres", ylabel="Desv. pp anual del EE",
        legend=:topright, margin=5P.mm)
    P.plot!(p_gs, periods, ctx.pi_agg_irf,   label="Agregada", color=:black, lw=2.5)
    P.plot!(p_gs, periods, ctx.pi_goods_irf, label="Bienes (sectores 1–5)", color=:steelblue, lw=2, ls=:dash)
    P.plot!(p_gs, periods, ctx.pi_serv_irf,  label="Servicios (sectores 6–12)", color=:firebrick, lw=2, ls=:dot)
    P.hline!(p_gs, [0.0], color=:black, lw=0.5, ls=:dash, label="")
    shock_save_fig(p_gs, fn[:gs_inflation], ctx)

    # ---- Mercado laboral agregado (1×2) ---- #
    p_lab = P.plot(layout=(1,2), size=(1000, 400),
        plot_title="$(ctx.title_short) — Mercado Laboral Agregado",
        titlefontsize=11, margin=5P.mm)
    P.plot!(p_lab, periods, ctx.n_irf_v, subplot=1, label=base_lbl, color=base_col, lw=base_lw,
        xlabel="Trimestres", ylabel="% desv. del EE", title="Empleo Agregado")
    P.plot!(p_lab, periods, ctx.w_irf, subplot=2, label=base_lbl, color=base_col, lw=base_lw,
        xlabel="Trimestres", ylabel="% desv. del EE", title="Salario Real")
    P.hline!(p_lab, [0.0], subplot=1, color=:black, lw=0.6, ls=:dash, label="")
    P.hline!(p_lab, [0.0], subplot=2, color=:black, lw=0.6, ls=:dash, label="")
    shock_save_fig(p_lab, fn[:labor_agg], ctx)

    # ---- Brecha del PIB ---- #
    p_gap = P.plot(size=(700, 400), title="$(ctx.title_short) — Brecha del PIB",
        titlefontsize=11, margin=5P.mm, xlabel="Trimestres", ylabel="% desv. del EE")
    P.plot!(p_gap, periods, ctx.gdpgap_irf, label=base_lbl, color=base_col, lw=base_lw)
    P.hline!(p_gap, [0.0], color=:black, lw=0.6, ls=:dash, label="")
    shock_save_fig(p_gap, fn[:gdpgap_agg], ctx)

    # ---- Mapa de exposición (barra específica del shock) ---- #
    p_exp = P.bar(1:nsec, ctx.exposure_vec,
        xticks=(1:nsec, bar_names),
        xrotation=45, label=ctx.exposure_label, color=:darkorange,
        ylabel=ctx.exposure_ylabel, title=ctx.exposure_title,
        titlefontsize=14, size=(900, 450), bottom_margin=10P.mm)
    shock_save_fig(p_exp, fn[:exposure], ctx)

    # ---- Descomposición EG del costo marginal (impacto / 6m / 12m) ---- #
    cc, cl = ctx.ge_colors, ctx.ge_labels
    p1 = make_ge_decomp_fig(ctx.ge_h1.comps, cc, cl, ctx.infl_irf_impact,
        "Inflación precios internos (pp t-a-t, t=1)",
        "$(ctx.title_short) — Descomposición EG del Costo Marginal (Impacto, εY=$(epsY_str))",
        "Puntos porcentuales (t = 1)", bar_names, nsec)
    shock_save_fig(p1, fn[:decomp_mc_impact], ctx)
    p2 = make_ge_decomp_fig(ctx.ge_h2.comps, cc, cl, ctx.infl_6m,
        "Inflación precios internos (acum. 6 meses, pp)",
        "$(ctx.title_short) — Descomposición EG del Costo Marginal (6 meses, εY=$(epsY_str))",
        "Puntos porcentuales", bar_names, nsec)
    shock_save_fig(p2, fn[:decomp_mc_6m], ctx)
    p4 = make_ge_decomp_fig(ctx.ge_h4.comps, cc, cl, ctx.infl_12m,
        "Inflación precios internos (acum. 12 meses, pp)",
        "$(ctx.title_short) — Descomposición EG del Costo Marginal (12 meses, εY=$(epsY_str))",
        "Puntos porcentuales", bar_names, nsec)
    shock_save_fig(p4, fn[:decomp_mc_12m], ctx)

    # ---- Descomposición agrupada en 3: Directo / Indirecto / Otros ---- #
    # Coarser version of the GE marginal-cost decomposition above, grouping the
    # five channels into the direct cost push, the network (indirect) propagation
    # via input prices, and everything else (imports + labour + residual).
    g3c, g3l = group3_style()
    p1g = make_ge_decomp_fig(group3_components(ctx.ge_h1.comps), g3c, g3l, ctx.infl_irf_impact,
        "Inflación precios internos (pp t-a-t, t=1)",
        "$(ctx.title_short) — Descomposición del Costo Marginal: Directo/Indirecto/Otros (Impacto, εY=$(epsY_str))",
        "Puntos porcentuales (t = 1)", bar_names, nsec)
    shock_save_fig(p1g, fn[:decomp_mc3_impact], ctx)
    p2g = make_ge_decomp_fig(group3_components(ctx.ge_h2.comps), g3c, g3l, ctx.infl_6m,
        "Inflación precios internos (acum. 6 meses, pp)",
        "$(ctx.title_short) — Descomposición del Costo Marginal: Directo/Indirecto/Otros (6 meses, εY=$(epsY_str))",
        "Puntos porcentuales", bar_names, nsec)
    shock_save_fig(p2g, fn[:decomp_mc3_6m], ctx)
    p4g = make_ge_decomp_fig(group3_components(ctx.ge_h4.comps), g3c, g3l, ctx.infl_12m,
        "Inflación precios internos (acum. 12 meses, pp)",
        "$(ctx.title_short) — Descomposición del Costo Marginal: Directo/Indirecto/Otros (12 meses, εY=$(epsY_str))",
        "Puntos porcentuales", bar_names, nsec)
    shock_save_fig(p4g, fn[:decomp_mc3_12m], ctx)

    # ---- Producto (valor agregado) e inflación sectorial en el impacto ---- #
    # Report sectoral VALUE ADDED under the "Y_i" label: gross output double-counts
    # intermediates and can rise under a negative TFP shock even as GDP falls.
    # Falls back to gross output if the caller did not supply a VA vector.
    out_impact = get(ctx, :va_irf_impact, ctx.y_irf_impact)
    pi_sec_impact = ctx.pi_sec_mat[:, 1]
    p_out = Main.groupedbar(
        [out_impact pi_sec_impact],
        xticks=(1:nsec, bar_names), xrotation=55,
        label=["Producto Y_i (% desv.)" "Inflación π_i (pp anual)"], color=[:steelblue :firebrick],
        ylabel="Respuesta en el impacto (t = 1)",
        title="$(ctx.title_short) — Impacto: Producto e Inflación Sectorial (εY=$(epsY_str))",
        titlefontsize=12, size=(1400, 600), legend=:topright,
        ylims=padlims(vcat(out_impact, pi_sec_impact)),
        bottom_margin=24P.mm, left_margin=14P.mm, right_margin=5P.mm)
    shock_save_fig(p_out, fn[:decomp_output], ctx)

    @printf "  Figuras de descomposición guardadas (base)\n"
    return nothing
end


# =========================================================================== #
#  TABLES + CSV                                                                #
# =========================================================================== #

function generate_shock_tables(ctx)
    nsec = ctx.nsec
    tab_names = SHOCK_TAB_NAMES
    out_impact = get(ctx, :va_irf_impact, ctx.y_irf_impact)   # VA (GDP-consistent) if supplied
    @printf "\n--- Generating LaTeX tables (%s) ---\n" ctx.title_short

    # ---- Table 1: Sectoral decomposition ---- #
    decomp_path = joinpath(ctx.TABLES_DIR, "$(ctx.tag)_shock_decomposition.tex")
    open(decomp_path, "w") do f
        println(f, "\\begin{table}[htbp]")
        println(f, "\\centering")
        println(f, "\\caption{$(ctx.decomp_caption)}")
        println(f, "\\label{tab:$(ctx.tag)_decomp}")
        println(f, "\\begin{threeparttable}")
        println(f, "\\small")
        println(f, "\\begin{tabular}{@{}l c c c c c c@{}}")
        println(f, "\\toprule")
        println(f, " & $(ctx.decomp_share_head) & \\multicolumn{3}{c}{Aumento CM (\\%)} & Amplif. & Producto \\\\")
        println(f, "\\cmidrule(lr){3-5}")
        println(f, "Sector & $(ctx.decomp_share_sub) & Directo & Red & Total & Razón & IRF (\\%) \\\\")
        println(f, "\\midrule")
        for i in 1:nsec
            gs = ctx.goods[i] ? "G" : "S"
            dval = abs(ctx.direct_mc[i]) < 1e-9 ? 0.0 : ctx.direct_mc[i]   # avoid "-0.000"
            amp_str = dval > 1e-8 ? @sprintf("%5.2f", ctx.amp_ratio[i]) : "---"
            @printf(f, "%s (%s) & %4.1f & %.3f & %.3f & %.3f & %s & \$%+.3f\$ \\\\\n",
                tab_names[i], gs, ctx.share_vec[i]*100, dval,
                ctx.network_mc[i], ctx.total_mc[i], amp_str, out_impact[i])
        end
        println(f, "\\midrule")
        @printf(f, "Agregado (ponderado por producto) & --- & %.3f & %.3f & %.3f & %5.2f & \$%+.3f\$ \\\\\n",
            ctx.agg_direct, ctx.agg_network, ctx.agg_total, ctx.agg_amp, ctx.gdp_irf[1])
        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\begin{tablenotes}[flushleft]")
        println(f, "\\footnotesize")
        println(f, "\\item \\textit{Notas:} $(ctx.decomp_notes)")
        println(f, "\\end{tablenotes}")
        println(f, "\\end{threeparttable}")
        println(f, "\\end{table}")
    end
    @printf "  Saved: %s\n" basename(decomp_path)
    shock_sync_table(basename(decomp_path), ctx)

    # ---- Table 2: Aggregate responses ---- #
    agg_path = joinpath(ctx.TABLES_DIR, "$(ctx.tag)_shock_aggregates.tex")
    open(agg_path, "w") do f
        println(f, "\\begin{table}[htbp]")
        println(f, "\\centering")
        println(f, "\\caption{$(ctx.agg_caption)}")
        println(f, "\\label{tab:$(ctx.tag)_agg}")
        println(f, "\\begin{threeparttable}")
        println(f, "\\begin{tabular}{@{}l c c c@{}}")
        println(f, "\\toprule")
        println(f, "Variable & Impacto & Mín/Máx & Trimestre \\\\")
        println(f, "\\midrule")
        vars_agg = [
            ("PIB",                ctx.gdp_irf,                          "\\% desv."),
            ("Inflación IPC",      ctx.pi_irf .* 4,                      "pp anual"),
            ("Consumo",            ctx.c_irf,                            "\\% desv."),
            ("Empleo",             ctx.n_irf_v,                          "\\% desv."),
            ("Tipo de Cambio Real",ctx.q_irf,                            "\\% desv."),
            ("Balanza Comercial",  ctx.tb_irf .* (ctx.TB_ss/ctx.GDP_ss), "pp del PIB"),
            ("Tasa de Política",   ctx.r_irf .* 4,                       "pp anual"),
        ]
        for (label, irf_v, _) in vars_agg
            impact = irf_v[1]
            ext_val = abs(minimum(irf_v)) > abs(maximum(irf_v)) ? minimum(irf_v) : maximum(irf_v)
            ext_q   = abs(minimum(irf_v)) > abs(maximum(irf_v)) ? argmin(irf_v) : argmax(irf_v)
            @printf(f, "%s & %+.3f & %+.3f & %d \\\\\n", label, impact, ext_val, ext_q)
        end
        println(f, "\\bottomrule")
        println(f, "\\end{tabular}")
        println(f, "\\begin{tablenotes}[flushleft]")
        println(f, "\\footnotesize")
        println(f, "\\item \\textit{Notas:} $(ctx.agg_notes)")
        println(f, "\\end{tablenotes}")
        println(f, "\\end{threeparttable}")
        println(f, "\\end{table}")
    end
    @printf "  Saved: %s\n" basename(agg_path)
    shock_sync_table(basename(agg_path), ctx)

    # ---- Optional hand-written section text ---- #
    let sec_src = joinpath(ctx.TABLES_DIR, "$(ctx.tag)_shock_section.tex")
        if isfile(sec_src)
            @printf "  Saved: %s (local)\n" basename(sec_src)
            shock_sync_table(basename(sec_src), ctx)
        end
    end

    # ---- IRF CSV (long format) ---- #
    let
        vars_to_save = ["GDP", "pi", "Q", "TB", "C", "N", "w", "r"]
        for i in 1:nsec
            push!(vars_to_save, "Y_$(i)", "PH_$(i)", "MC_$(i)", "L_$(i)")
        end
        rows = []
        for vn in vars_to_save
            irf_v = ctx.get_irf(vn)
            for h in 1:ctx.nT
                push!(rows, (period=h, variable=vn, value=irf_v[h]))
            end
        end
        df_out = Main.DataFrame(rows)
        Main.CSV.write(joinpath(ctx.TABLES_DIR, "$(ctx.tag)_shock_irfs.csv"), df_out)
        @printf "  Saved: %s_shock_irfs.csv\n" ctx.tag
    end

    # ---- Decomposition CSV ---- #
    df_decomp = Main.DataFrame(
        sector=1:nsec, name=ctx.names_vec, share=ctx.share_vec,
        alphaV=ctx.modalphaV, direct_mc=ctx.direct_mc, network_mc=ctx.network_mc,
        total_mc=ctx.total_mc, amp_ratio=ctx.amp_ratio, mc_irf=ctx.mc_irf_impact,
        y_irf=ctx.y_irf_impact, ph_irf=ctx.ph_irf_impact, yi_ss=ctx.Yi_ss)
    Main.CSV.write(joinpath(ctx.TABLES_DIR, "$(ctx.tag)_shock_decomposition.csv"), df_decomp)
    @printf "  Saved: %s_shock_decomposition.csv\n" ctx.tag
    return nothing
end


# =========================================================================== #
#  ORCHESTRATOR                                                                #
# =========================================================================== #

function generate_shock_outputs(ctx)
    if ctx.HAS_PLOTS
        generate_shock_figures(ctx)
    else
        @printf "\n  [Plots.jl not available — skipping figures]\n"
    end
    generate_shock_tables(ctx)
    return nothing
end
