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
  * `standard_fignames(tag)`             — `<base>_<tag>_shock.png` naming scheme
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
#  EXOGENOUS-SHOCK COLUMN LOOKUP (name -> ghu column)                          #
# =========================================================================== #
# ghu / g1_3 columns follow the exogenous DECLARATION order, which is exactly the
# modfile.json "exogenous" array. Resolving the shock column BY NAME keeps the
# IRF scripts correct regardless of the .mod shock ordering (e.g. after adding the
# 12 sectoral demand shocks, eps_postar moved from column 5 to column 4).

const _SHOCK_EXO_NAMES = Ref{Vector{String}}(String[])

function exo_shock_names(mod_dir::AbstractString)
    isempty(_SHOCK_EXO_NAMES[]) || return _SHOCK_EXO_NAMES[]
    mf = joinpath(mod_dir, "NK_SOE_lev_gap2", "model", "json", "modfile.json")
    names = String[]
    if isfile(mf)
        try
            raw = read(mf, String)
            key = findfirst("\"exogenous\"", raw)   # exact key (not exogenous_deterministic)
            if key !== nothing
                lb = findnext('[', raw, key[end] + 1)
                depth = 0; rb = lb
                for p in lb:lastindex(raw)
                    c = raw[p]
                    c == '[' && (depth += 1)
                    if c == ']'
                        depth -= 1
                        depth == 0 && (rb = p; break)
                    end
                end
                names = [String(m.captures[1]) for m in eachmatch(r"\"name\"\s*:\s*\"([^\"]+)\"", raw[lb:rb])]
            end
        catch; end
    end
    _SHOCK_EXO_NAMES[] = names
    return names
end

"""
    exo_col(shock_name, mod_dir) -> Int

Column of `shock_name` in ghu (= position in the model's exogenous declaration
order). Errors clearly if the shock is absent (usually means the context is stale
— rebuild with main_SOE_gap.jl).
"""
function exo_col(shock_name::AbstractString, mod_dir::AbstractString)
    exo = exo_shock_names(mod_dir)
    k = findfirst(==(shock_name), exo)
    k === nothing && error(
        "Shock '$shock_name' not found in the model's exogenous list " *
        "($(length(exo)) shocks: $(isempty(exo) ? "<none — modfile.json missing?>" : join(exo, ", "))). " *
        "Recompile the unified model: julia --project=. main_SOE_gap.jl")
    return k
end

# =========================================================================== #
#  FILENAME SCHEMES                                                            #
# =========================================================================== #

"Standardised `<base>_<tag>_shock.png` filenames used by the TFP shock scripts."
function standard_fignames(tag::AbstractString)
    Dict{Symbol,String}(
        :agg              => "irf_aggregate_$(tag)_shock.png",
        :sec_Y            => "irf_sectoral_Y_$(tag)_shock.png",
        :sec_PH           => "irf_sectoral_PH_$(tag)_shock.png",
        :sec_inflation    => "irf_sectoral_inflation_$(tag)_shock.png",
        :agg_inflation    => "irf_aggregate_inflation_$(tag)_shock.png",
        :gs_inflation     => "irf_goods_vs_services_inflation_$(tag)_shock.png",
        :affected_inflation => "irf_affected_vs_other_inflation_$(tag)_shock.png",
        :sec_MC           => "irf_sectoral_MC_$(tag)_shock.png",
        :labor_agg        => "irf_labor_aggregate_$(tag)_shock.png",
        :sec_L            => "irf_sectoral_L_$(tag)_shock.png",
        :gdpgap_agg       => "irf_gdpgap_aggregate_$(tag)_shock.png",
        :sec_Ygap         => "irf_sectoral_Ygap_$(tag)_shock.png",
        :exposure         => "$(tag)_intensity_exposure.png",
        :decomp_mc_impact => "decomposition_mc_inflation_$(tag)_baseline.png",
        :decomp_mc_6m     => "decomposition_mc_inflation_6m_$(tag)_baseline.png",
        :decomp_mc_12m    => "decomposition_mc_inflation_12m_$(tag)_baseline.png",
        # 3-way grouped version (direct / indirect / others)
        :decomp_mc3_impact => "decomposition_mc3_inflation_$(tag)_baseline.png",
        :decomp_mc3_6m     => "decomposition_mc3_inflation_6m_$(tag)_baseline.png",
        :decomp_mc3_12m    => "decomposition_mc3_inflation_12m_$(tag)_baseline.png",
        # sectoral-INFLATION (first-difference) versions of the price decompositions
        :decomp_pi_impact  => "decomposition_pi_$(tag)_baseline.png",
        :decomp_pi_6m      => "decomposition_pi_6m_$(tag)_baseline.png",
        :decomp_pi_12m     => "decomposition_pi_12m_$(tag)_baseline.png",
        :decomp_pi3_impact => "decomposition_pi3_$(tag)_baseline.png",
        :decomp_pi3_6m     => "decomposition_pi3_6m_$(tag)_baseline.png",
        :decomp_pi3_12m    => "decomposition_pi3_12m_$(tag)_baseline.png",
        # simplified 2-way MC version (direct / encadenamientos only)
        :decomp_mc2_impact => "decomposition_mc2_$(tag)_baseline.png",
        :decomp_mc2_6m     => "decomposition_mc2_6m_$(tag)_baseline.png",
        :decomp_mc2_12m    => "decomposition_mc2_12m_$(tag)_baseline.png",
        :decomp_output    => "decomposition_output_$(tag)_baseline.png",
    )
end

"Legacy oil filenames (kept verbatim so existing \\includegraphics in the paper resolve)."
function oil_fignames()
    Dict{Symbol,String}(
        :agg              => "irf_aggregate_oil_shock.png",
        :sec_Y            => "irf_sectoral_Y_oil_shock.png",
        :sec_PH           => "irf_sectoral_PH_oil_shock.png",
        :sec_inflation    => "irf_sectoral_inflation_oil.png",
        :agg_inflation    => "irf_aggregate_inflation_oil.png",
        :gs_inflation     => "irf_goods_vs_services_inflation_oil.png",
        :affected_inflation => "irf_affected_vs_other_inflation_oil.png",
        :sec_MC           => "irf_sectoral_MC_oil_shock.png",
        :labor_agg        => "irf_labor_aggregate_oil_shock.png",
        :sec_L            => "irf_sectoral_L_oil_shock.png",
        :gdpgap_agg       => "irf_gdpgap_aggregate_oil_shock.png",
        :sec_Ygap         => "irf_sectoral_Ygap_oil_shock.png",
        :exposure         => "oil_intensity_exposure.png",
        :decomp_mc_impact => "decomposition_mc_inflation_oil_baseline.png",
        :decomp_mc_6m     => "decomposition_mc_inflation_6m_oil_baseline.png",
        :decomp_mc_12m    => "decomposition_mc_inflation_12m_oil_baseline.png",
        # 3-way grouped version (direct / indirect / others)
        :decomp_mc3_impact => "decomposition_mc3_inflation_oil_baseline.png",
        :decomp_mc3_6m     => "decomposition_mc3_inflation_6m_oil_baseline.png",
        :decomp_mc3_12m    => "decomposition_mc3_inflation_12m_oil_baseline.png",
        # sectoral-INFLATION (first-difference) versions of the price decompositions
        :decomp_pi_impact  => "decomposition_pi_oil_baseline.png",
        :decomp_pi_6m      => "decomposition_pi_6m_oil_baseline.png",
        :decomp_pi_12m     => "decomposition_pi_12m_oil_baseline.png",
        :decomp_pi3_impact => "decomposition_pi3_oil_baseline.png",
        :decomp_pi3_6m     => "decomposition_pi3_6m_oil_baseline.png",
        :decomp_pi3_12m    => "decomposition_pi3_12m_oil_baseline.png",
        # simplified 2-way MC version (direct / encadenamientos only)
        :decomp_mc2_impact => "decomposition_mc2_oil_baseline.png",
        :decomp_mc2_6m     => "decomposition_mc2_6m_oil_baseline.png",
        :decomp_mc2_12m    => "decomposition_mc2_12m_oil_baseline.png",
        :decomp_output    => "decomposition_output_oil_baseline.png",
    )
end

# Standard sector display names in Spanish (consistent across every figure/table).
const SHOCK_BAR_NAMES = [
    "Agric. y Pesca", "Minería", "Manufactura", "EGA",
    "Construcción", "Comercio y Hot.", "Transp. y Com.", "Finanzas",
    "Inmobiliario", "Serv. Empres.", "Serv. Pers.", "Adm. Pública",
]
const SHOCK_PANEL_NAMES = [
    "Agricultura y Pesca", "Minería", "Manufactura", "EGA",
    "Construcción", "Comercio y Hoteles", "Transporte y Com.", "Finanzas",
    "Inmobiliario", "Serv. Empresariales", "Serv. Personales", "Adm. Pública",
]
const SHOCK_TAB_NAMES = [
    "Agricultura y Pesca", "Minería", "Manufactura", "Electricidad",
    "Construcción", "Comercio y Hoteles", "Transporte y Com.", "Finanzas",
    "Inmobiliario", "Serv.\\ Empresariales", "Serv.\\ Personales", "Adm.\\ Pública",
]

# --------------------------------------------------------------------------- #
#  IPoM (BCCh) chart palette — sampled from IPoM Diciembre 2025, Gráfico I.8   #
#  ("Contribuciones a la variación anual del IPC total", incidence bars).      #
#  Black diamond = total, as in the IPoM incidence charts.                     #
# --------------------------------------------------------------------------- #
const IPOM_NAVY      = "#1F2E58"   # dark navy   (Energía volátiles)
const IPOM_LIGHTBLUE = "#4DB3E9"   # light blue  (Resto volátiles)
const IPOM_RED       = "#C21E24"   # red         (Servicios sin volátiles)
const IPOM_ORANGE    = "#ED6331"   # orange      (Alimentos volátiles)
const IPOM_GREEN     = "#3BB957"   # green       (Bienes sin volátiles)
const IPOM_YELLOW    = "#E7EB13"   # yellow      (Alimentos sin volátiles)


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
            infl_irf_impact=infl_irf_impact, infl_6m=infl_6m, infl_12m=infl_12m,
            cons_ss_all=cons_ss_all)
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

6-way decomposition of the sectoral relative home-price IRF at horizon `h` into
cost-push channels plus the markup.  Two `kind`s:

  :price (oil)  → [direct oil, network via prices, non-oil imports, labor, residual, markup]
  :tfp          → [own TFP, network via prices, imports, labor, residual, markup]

The first five components sum exactly to the marginal-cost IRF (the residual is
defined as the gap to m̂c).  Adding the markup component μ̂ = p̂H − m̂c, the full
stack sums exactly to the relative home price p̂H, returned as `price`.

Returns `(comps = (c1,...,c6), mc = mc_h, price = ph_h)` plus colour/label
vectors via the companion `ge_component_style(kind)`.
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
    ph_h  = ph_irf_mat[:, h]
    mu_h  = ph_h .- mc_h                              # markup: p̂H − m̂c
    return (comps=(copy(dir_h), copy(net_h), copy(imp_h), copy(lab_h),
                   copy(res_h), copy(mu_h)),
            mc=copy(mc_h), price=copy(ph_h))
end

"Colours (IPoM palette) and legend labels for the 6-way price decomposition bars."
function ge_component_style(kind::Symbol)
    colors = [IPOM_NAVY, IPOM_LIGHTBLUE, IPOM_GREEN, IPOM_ORANGE, IPOM_YELLOW, IPOM_RED]
    labels = kind === :price ?
        ["Petróleo directo", "Red vía precios", "Import. no-petróleo",
         "Trabajo", "Residual", "Margen"] :
        ["PTF propia", "Red vía precios", "Importaciones",
         "Trabajo", "Residual", "Margen"]
    return colors, labels
end

"""
    group3_components(comps6)

Collapse the 6-way GE price decomposition into four groups for the coarser
"direct / indirect / markup / others" figure:
  * Directo   = own cost push (component 1: own TFP, or direct oil)
  * Indirecto = network propagation via input prices (component 2)
  * Margen    = markup p̂H − m̂c (component 6)
  * Otros     = imports + labour + residual (components 3 + 4 + 5)
`comps6` is the 6-tuple of `nsec`-vectors returned in `ge_mc_components(...).comps`.
Returns a 4-tuple `(direct, indirect, markup, others)`.
"""
function group3_components(comps6)
    direct   = copy(comps6[1])
    indirect = copy(comps6[2])
    markup   = copy(comps6[6])
    others   = comps6[3] .+ comps6[4] .+ comps6[5]
    return (direct, indirect, markup, others)
end

"Colours (IPoM palette) and legend labels for the grouped (direct / indirect / markup / others) bars."
function group3_style()
    colors = [IPOM_NAVY, IPOM_LIGHTBLUE, IPOM_RED, IPOM_ORANGE]
    labels = ["Directo", "Encadenamientos", "Márgenes + Expectativas", "Otros Costos Marginales"]
    return colors, labels
end

"Colours and labels for the simplified 2-way (direct / encadenamientos) bars."
function group2_style()
    colors = [IPOM_NAVY, IPOM_LIGHTBLUE]
    labels = ["Efectos directos del shock", "Efectos de encadenamiento"]
    return colors, labels
end

"""
    ge_inflation_components(geh, gehm1, pi_h, nsec)

First-difference (INFLATION) version of the GE price decomposition at one
horizon.  Sectoral nominal home-price inflation (annualised, pp) is
    π_i[h] = 4·(p̂H_i[h] − p̂H_i[h−1] + π[h]),
where p̂H is the relative home price and π is aggregate CPI inflation.  Since
the price components sum exactly to p̂H at each horizon, differencing each of
the six components (×4) preserves additivity.  The common nominal-drift term
4·π[h] (identical across sectors) is folded into the RESIDUAL component rather
than shown as an explicit bar, so the stack still sums exactly to sectoral
inflation π_i[h] and the figure keeps the same 6 bars as the price version.
Pass `gehm1 = nothing` for h = 1 (pre-shock components are zero).
Returns a 6-tuple of `nsec`-vectors (order as in `ge_mc_components`).
"""
function ge_inflation_components(geh, gehm1, pi_h, nsec)
    base = gehm1 === nothing ?
        [4.0 .* c for c in geh.comps] :
        [4.0 .* (c .- cm) for (c, cm) in zip(geh.comps, gehm1.comps)]
    base[5] = base[5] .+ 4.0 * pi_h     # aggregate CPI inflation → residual
    return Tuple(base)
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

"""
    save_gs_inflation_excel(periods, ctx, fig_fname)

Write an Excel workbook with the exact series plotted in the goods-vs-services
home-price inflation figure (`fig_fname`). Columns: horizon plus aggregate,
goods and services inflation (deviations in annual pp from steady state). The
file is saved alongside the figure in `ctx.FIGURES_DIR` with the same base name
(`.png` → `.xlsx`). `XLSX`/`DataFrame` are resolved lazily from `Main`; if XLSX
is unavailable the data are written as CSV so nothing is lost.
"""
function save_gs_inflation_excel(periods, ctx, fig_fname)
    hor   = collect(periods)
    agg   = collect(Float64, ctx.pi_agg_irf[periods])
    goods = collect(Float64, ctx.pi_goods_irf[periods])
    serv  = collect(Float64, ctx.pi_serv_irf[periods])
    base  = replace(fig_fname, r"\.png$" => "")
    cols  = ["Trimestre", "Inflacion_Agregada", "Inflacion_Bienes", "Inflacion_Servicios"]
    data  = Any[hor, agg, goods, serv]
    if isdefined(Main, :XLSX)
        xlsx_path = joinpath(ctx.FIGURES_DIR, "$(base).xlsx")
        Main.XLSX.writetable(xlsx_path, data, cols; overwrite=true,
                             sheetname="inflacion_bienes_servicios")
        @printf "  Saved: %s\n" basename(xlsx_path)
    else
        csv_path = joinpath(ctx.FIGURES_DIR, "$(base).csv")
        Main.CSV.write(csv_path, Main.DataFrame(data, cols))
        @printf "  Saved: %s (XLSX no disponible, se guardó CSV)\n" basename(csv_path)
    end
    return nothing
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
                            title_str, ylabel_str, bar_names, nsec;
                            axis2_sectors=Int[])
    P = Main.Plots
    # Drop Public Administration from the decomposition figures.
    keep = [i for i in 1:nsec if !occursin("blica", bar_names[i])]
    comps = Tuple(c[keep] for c in comps)
    infl_h = infl_h[keep]
    bar_names = bar_names[keep]
    nsec = length(keep)
    comp_mat = hcat(comps...)                 # nsec × ncomp
    # ---- Optional secondary axis -------------------------------------------- #
    # The sectors in `axis2_sectors` (original 1:12 indices, e.g. the sectors
    # directly hit by the shock) are MOVED TO THE LEFTMOST bars and read on the
    # LEFT axis (true values).  The remaining sectors are drawn rescaled UP by a
    # round factor `s_ax2` and read on the right-hand axis (= left limits ÷
    # s_ax2), so their much smaller responses stay visible.  A dashed separator
    # and group labels (added below) mark which sectors read on which axis.
    ax2pos = findall(i -> i in axis2_sectors, keep)
    s_ax2 = 1.0
    if !isempty(ax2pos) && length(ax2pos) < nsec
        ord = vcat(ax2pos, setdiff(1:nsec, ax2pos))   # shocked sectors → leftmost
        comp_mat  = comp_mat[ord, :]
        infl_h    = infl_h[ord]
        bar_names = bar_names[ord]
        ax2pos = collect(1:length(ax2pos))
        oth    = setdiff(1:nsec, ax2pos)
        extreme(idx) = maximum(vcat(
            vec(sum(max.(comp_mat[idx, :], 0.0), dims=2)),
            .-vec(sum(min.(comp_mat[idx, :], 0.0), dims=2)),
            abs.(infl_h[idx]), [1e-12]))
        ratio = extreme(ax2pos) / extreme(oth)
        if ratio > 1.5
            # round the scale factor up to 1–2–5×10^k
            m = 10.0^floor(log10(ratio)); v = ratio / m
            s_ax2 = (v <= 1.0 ? 1.0 : v <= 2.0 ? 2.0 : v <= 5.0 ? 5.0 : 10.0) * m
            comp_mat[oth, :] .*= s_ax2     # draw the small group on the left scale
            infl_h[oth] .*= s_ax2          # right axis shows their true values
        end
    end
    pos_mat  = max.(comp_mat, 0.0)
    neg_mat  = min.(comp_mat, 0.0)
    pos_tops = vec(sum(pos_mat, dims=2))
    neg_tops = vec(sum(neg_mat, dims=2))
    all_vals = vcat(pos_tops, neg_tops, infl_h, [0.0])
    # Additive padding around the true data range, always including zero so the
    # negative (deflationary) part of the response is never clipped.  Extra top
    # padding when the dual-axis group labels are drawn inside the plot.
    dmin = minimum(all_vals); dmax = maximum(all_vals)
    span = max(dmax - dmin, 1e-6)
    ylo = min(dmin, 0.0) - 0.10 * span
    yhi = max(dmax, 0.0) + (s_ax2 > 1.0 ? 0.28 : 0.14) * span
    # IPoM-style: title + units line in navy (no ylabel), frameless legend on
    # top, no grid, solid thin zero line — mirrors the BCCh incidence charts.
    p = P.plot(
        xticks=(1:nsec, bar_names), xrotation=55,
        title="$(title_str)\n($(ylabel_str))",
        titlefontsize=12, titlefontcolor=IPOM_NAVY, titlelocation=:left,
        size=(1400, 760), legend=:outertop, legend_columns=3,
        legendfontsize=9, foreground_color_legend=nothing,
        background_color_legend=nothing, grid=false,
        ylims=(ylo, yhi),
        bottom_margin=20P.mm, left_margin=10P.mm,
        right_margin=(s_ax2 > 1.0 ? 14 : 8)P.mm, top_margin=4P.mm,
        xlims=(0.3, nsec + 0.7))
    P.hline!(p, [0.0], color=:black, lw=0.8, label="")
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
    # Right-hand axis for the rescaled (non-shocked) sectors: same axis area,
    # ylims divided by the scale factor, so those bars read true values on the
    # right.  Dashed separator + group labels show which axis applies.
    if s_ax2 > 1.0
        xsep = last(ax2pos) + 0.5
        P.vline!(p, [xsep], color=:gray40, ls=:dash, lw=1.2, label="")
        P.annotate!(p, (1 + last(ax2pos)) / 2, yhi - 0.04 * (yhi - ylo),
            P.text("Sectores afectados\n(eje izquierdo)", 9, :gray30, :center))
        P.annotate!(p, (xsep + 0.5 + nsec) / 2, yhi - 0.04 * (yhi - ylo),
            P.text("Resto de sectores (eje derecho)", 9, :gray30, :center))
        sp2 = P.twinx(p)
        P.plot!(sp2, [NaN], [NaN], label="", grid=false,
                ylims=(ylo / s_ax2, yhi / s_ax2), xlims=(0.3, nsec + 0.7))
    end
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

    # Sectoral "output" panel reports VALUE ADDED (GDP-consistent), matching the
    # decomposition_output bars: gross output Y_i double-counts intermediates and
    # can rise under a negative TFP shock even as VA and GDP fall.  Falls back to
    # gross output if the caller did not supply va_irf_mat.
    sectoral_panel(get(ctx, :va_irf_mat, ctx.y_irf_mat), "Producto Sectorial (Valor Agregado)", fn[:sec_Y])
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
    P.plot!(p_gs, periods, ctx.pi_goods_irf, label="Bienes (sectores 1–5)", color=IPOM_NAVY, lw=2, ls=:dash)
    P.plot!(p_gs, periods, ctx.pi_serv_irf,  label="Servicios (sectores 6–12)", color=IPOM_RED, lw=2, ls=:dot)
    P.hline!(p_gs, [0.0], color=:black, lw=0.5, ls=:dash, label="")
    shock_save_fig(p_gs, fn[:gs_inflation], ctx)
    # Excel con los datos exactos graficados en la figura bienes vs servicios.
    save_gs_inflation_excel(periods, ctx, fn[:gs_inflation])

    # ---- Inflación: sector(es) afectado(s) vs resto ---- #
    # Three lines: aggregate home-price inflation, the consumption-weighted
    # inflation of the sector(s) directly hit by the shock, and the
    # consumption-weighted inflation of all remaining ("other") sectors. Weights
    # are steady-state consumption shares, matching the goods-vs-services figure.
    aff = collect(Int, get(ctx, :affected_sectors, Int[]))
    if !isempty(aff)
        wgt = get(ctx, :cons_ss_all, ones(nsec))
        oth = setdiff(1:nsec, aff)
        waff = sum(wgt[aff]) > 1e-12 ? wgt[aff] ./ sum(wgt[aff]) : fill(1/length(aff), length(aff))
        pi_aff = [sum(waff .* ctx.pi_sec_mat[aff, h]) for h in periods]
        pi_oth = if isempty(oth)
            zeros(length(periods))
        else
            woth = sum(wgt[oth]) > 1e-12 ? wgt[oth] ./ sum(wgt[oth]) : fill(1/length(oth), length(oth))
            [sum(woth .* ctx.pi_sec_mat[oth, h]) for h in periods]
        end
        aff_lbl = get(ctx, :affected_label, "Sector(es) afectado(s)")
        # Explicit y-limits with extra bottom padding so the affected-sector
        # name footnote fits inside the axis without touching the lines.
        allv = vcat(ctx.pi_agg_irf[periods], pi_aff, pi_oth)
        dmn = minimum(allv); dmx = maximum(allv); sp = max(dmx - dmn, 1e-6)
        ylo = min(dmn, 0.0) - 0.22 * sp
        yhi = max(dmx, 0.0) + 0.12 * sp
        p_aff = P.plot(size=(800, 500),
            title="$(ctx.title_short) — Inflación de Precios Internos: Afectado vs Resto",
            titlefontsize=11, xlabel="Trimestres", ylabel="Desv. pp anual del EE",
            legend=:topright, ylims=(ylo, yhi), margin=5P.mm)
        P.plot!(p_aff, periods, ctx.pi_agg_irf, label="Agregada", color=:black, lw=2.5)
        P.plot!(p_aff, periods, pi_aff, label="$(aff_lbl) (afectado)", color=IPOM_RED, lw=2, ls=:dash)
        P.plot!(p_aff, periods, pi_oth, label="Otros sectores", color=IPOM_LIGHTBLUE, lw=2, ls=:dot)
        P.hline!(p_aff, [0.0], color=:black, lw=0.5, ls=:dash, label="")
        # Footnote inside the axis naming the affected sectors explicitly
        # (essential for the oil shock, where the set is data-determined).
        aff_names = join(SHOCK_BAR_NAMES[aff], ", ")
        P.annotate!(p_aff, (first(periods) + last(periods)) / 2, ylo + 0.05 * (yhi - ylo),
            P.text("Sectores afectados: $(aff_names)", 9, :gray30, :center))
        shock_save_fig(p_aff, fn[:affected_inflation], ctx)
    end

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

    # ---- Mapa de exposición (barra específica del shock, estilo IPoM) ---- #
    # Same house style as the GE-decomposition bars (make_ge_decomp_fig): navy
    # left-aligned title with the units on a second line, frameless top legend,
    # no grid, solid thin black zero line, IPoM-navy bars with thin white edges.
    # Public Administration is dropped so the x-axis matches the decomposition
    # figures (decomposition_pi3_*), with which this map is read side by side.
    exp_keep  = [i for i in 1:nsec if !occursin("blica", bar_names[i])]
    exp_vals  = ctx.exposure_vec[exp_keep]
    exp_names = bar_names[exp_keep]
    nexp      = length(exp_keep)
    exp_dmax  = maximum(vcat(exp_vals, 0.0)); exp_span = max(exp_dmax, 1e-6)
    p_exp = P.plot(
        xticks=(1:nexp, exp_names), xrotation=55,
        title="$(ctx.exposure_title)\n($(ctx.exposure_ylabel))",
        titlefontsize=12, titlefontcolor=IPOM_NAVY, titlelocation=:left,
        size=(1400, 760), legend=:outertop, legend_columns=3,
        legendfontsize=9, foreground_color_legend=nothing,
        background_color_legend=nothing, grid=false,
        ylims=(0.0, exp_dmax + 0.14 * exp_span),
        bottom_margin=20P.mm, left_margin=10P.mm, right_margin=8P.mm,
        top_margin=4P.mm, xlims=(0.3, nexp + 0.7))
    P.hline!(p_exp, [0.0], color=:black, lw=0.8, label="")
    lbl_used = false
    for i in 1:nexp
        s = shock_bar_rect(i, 0.0, exp_vals[i], 0.65)
        P.plot!(p_exp, s, color=IPOM_NAVY, label=(lbl_used ? "" : ctx.exposure_label),
                alpha=0.85, linecolor=:white, linewidth=0.3)
        lbl_used = true
    end
    shock_save_fig(p_exp, fn[:exposure], ctx)

    # ---- Descomposición EG del precio interno (impacto / 6m / 12m) ---- #
    # Bars: 5 marginal-cost channels + markup (p̂H − m̂c); they sum exactly to
    # the relative home-price IRF p̂H, shown as the black diamond.
    cc, cl = ctx.ge_colors, ctx.ge_labels
    p1 = make_ge_decomp_fig(ctx.ge_h1.comps, cc, cl, ctx.ge_h1.price,
        "Precio interno relativo al impacto",
        "$(ctx.title_short) — Descomposición EG del Precio Interno (Impacto)",
        "desviación del EE, puntos porcentuales", bar_names, nsec)
    shock_save_fig(p1, fn[:decomp_mc_impact], ctx)
    p2 = make_ge_decomp_fig(ctx.ge_h2.comps, cc, cl, ctx.ge_h2.price,
        "Precio interno relativo a 6 meses",
        "$(ctx.title_short) — Descomposición EG del Precio Interno (6 meses)",
        "desviación del EE, puntos porcentuales", bar_names, nsec)
    shock_save_fig(p2, fn[:decomp_mc_6m], ctx)
    p4 = make_ge_decomp_fig(ctx.ge_h4.comps, cc, cl, ctx.ge_h4.price,
        "Precio interno relativo a 12 meses",
        "$(ctx.title_short) — Descomposición EG del Precio Interno (12 meses)",
        "desviación del EE, puntos porcentuales", bar_names, nsec)
    shock_save_fig(p4, fn[:decomp_mc_12m], ctx)

    # ---- Descomposición agrupada: Directo / Indirecto / Margen / Otros ---- #
    # Coarser version of the GE price decomposition above, grouping the channels
    # into the direct cost push, the network (indirect) propagation via input
    # prices, the markup, and everything else (imports + labour + residual).
    # Same dual-axis layout as the inflation decompositions: the shocked
    # sector(s) (ctx.pi_decomp_axis2) read on the left axis, the rest rescaled
    # onto a right-hand axis so their smaller responses stay legible.
    ax2 = collect(Int, get(ctx, :pi_decomp_axis2, Int[]))
    g3c, g3l = group3_style()
    p1g = make_ge_decomp_fig(group3_components(ctx.ge_h1.comps), g3c, g3l, ctx.ge_h1.price,
        "Precio interno relativo al impacto",
        "$(ctx.title_short) — Descomposición del Precio Interno (Impacto)",
        "desviación del EE, puntos porcentuales", bar_names, nsec; axis2_sectors=ax2)
    shock_save_fig(p1g, fn[:decomp_mc3_impact], ctx)
    p2g = make_ge_decomp_fig(group3_components(ctx.ge_h2.comps), g3c, g3l, ctx.ge_h2.price,
        "Precio interno relativo a 6 meses",
        "$(ctx.title_short) — Descomposición del Precio Interno (6 meses)",
        "desviación del EE, puntos porcentuales", bar_names, nsec; axis2_sectors=ax2)
    shock_save_fig(p2g, fn[:decomp_mc3_6m], ctx)
    p4g = make_ge_decomp_fig(group3_components(ctx.ge_h4.comps), g3c, g3l, ctx.ge_h4.price,
        "Precio interno relativo a 12 meses",
        "$(ctx.title_short) — Descomposición del Precio Interno (12 meses)",
        "desviación del EE, puntos porcentuales", bar_names, nsec; axis2_sectors=ax2)
    shock_save_fig(p4g, fn[:decomp_mc3_12m], ctx)

    # ---- Descomposición EG de la INFLACIÓN sectorial (impacto / 6m / 12m) ---- #
    # First-difference version of the price decompositions above: the bars stack
    # the change in each cost channel (annualised, ×4); the common aggregate CPI
    # term 4·π[h] is folded into the residual bar, so the stack sums exactly to
    # sectoral home-price inflation π_i (black diamond = ctx.pi_sec_mat) with the
    # same 6 bars and colours as the price figures.  For the sectoral TFP shocks
    # the caller supplies `pi_decomp_axis2` (the shocked sector(s)), which are
    # rescaled onto a right-hand axis so the smaller responses stay legible.
    ge_fun = get(ctx, :ge_fun, nothing)
    if ge_fun !== nothing
        ge_h3 = ge_fun(3)                       # needed for the 12-month difference
        g2c, g2l = group2_style()
        for (h, geh, gehm1, k6, k4, k2, lab, dlbl) in (
                (1, ctx.ge_h1, nothing,    :decomp_pi_impact, :decomp_pi3_impact, :decomp_mc2_impact, "Impacto",  "Inflación Sectorial al Impacto"),
                (2, ctx.ge_h2, ctx.ge_h1,  :decomp_pi_6m,     :decomp_pi3_6m,     :decomp_mc2_6m,     "6 meses",  "Inflación Sectorial a 6 Meses"),
                (4, ctx.ge_h4, ge_h3,      :decomp_pi_12m,    :decomp_pi3_12m,    :decomp_mc2_12m,    "12 meses", "Inflación Sectorial a 12 Meses"))
            comps6 = ge_inflation_components(geh, gehm1, ctx.pi_irf[h], nsec)
            pih    = ctx.pi_sec_mat[:, h]
            p6 = make_ge_decomp_fig(comps6, cc, cl, pih, dlbl,
                "$(ctx.title_short) — Descomposición EG de la Inflación Sectorial ($(lab))",
                "desviación del EE, puntos porcentuales anualizados", bar_names, nsec; axis2_sectors=ax2)
            shock_save_fig(p6, fn[k6], ctx)
            comps4 = group3_components(comps6)
            p4c = make_ge_decomp_fig(comps4, g3c, g3l, pih, dlbl,
                "$(ctx.title_short) — Descomposición de la Inflación Sectorial ($(lab))",
                "desviación del EE, puntos porcentuales anualizados", bar_names, nsec; axis2_sectors=ax2)
            shock_save_fig(p4c, fn[k4], ctx)
            # Simplified 2-way MC version: only the direct cost push and the
            # network (encadenamientos) bars from the MARGINAL-COST
            # decomposition (levels, not first differences); the black diamond
            # is the relative home-price response p̂H ("Efecto del shock en
            # precios"), which sits near zero because prices are sticky —
            # large cost shocks, small price pass-through.  Single shared axis
            # for all sectors (no right-hand rescaled axis).
            p2c = make_ge_decomp_fig((geh.comps[1], geh.comps[2]), g2c, g2l, geh.price,
                "Efecto del shock en precios",
                "$(ctx.title_short) — Efectos Directos y de Encadenamiento en el Costo Marginal Sectorial ($(lab))",
                "desviación del EE, puntos porcentuales", bar_names, nsec)
            shock_save_fig(p2c, fn[k2], ctx)
        end
    end

    # ---- Producto (valor agregado) e inflación sectorial en el impacto ---- #
    # Report sectoral VALUE ADDED under the "Y_i" label: gross output double-counts
    # intermediates and can rise under a negative TFP shock even as GDP falls.
    # Falls back to gross output if the caller did not supply a VA vector.
    out_impact = get(ctx, :va_irf_impact, ctx.y_irf_impact)
    pi_sec_impact = ctx.pi_sec_mat[:, 1]
    p_out = Main.groupedbar(
        [out_impact pi_sec_impact],
        xticks=(1:nsec, bar_names), xrotation=55,
        label=["Producto Y_i (% desv.)" "Inflación π_i (pp anual)"], color=[IPOM_NAVY IPOM_RED],
        ylabel="Respuesta en el impacto (t = 1)",
        title="$(ctx.title_short) — Impacto: Producto e Inflación Sectorial",
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
