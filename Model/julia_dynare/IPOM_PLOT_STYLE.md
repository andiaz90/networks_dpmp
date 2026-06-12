---
name: ipom-policy-plots
description: >
  Agustín's house style for POLICY figures at the Central Bank of Chile —
  the IPoM / BCCh chart format used across the Networks-DPMP (NK-IOSOE) project
  (e.g. decomposition_pi3_oil_baseline). Use this skill EVERY time Agustín asks
  for a plot, chart, or figure "for policy purposes", "for the bank", "IPoM
  style", "BCCh style", "central-bank style", "presentation-ready", or wants a
  NEW figure to be CONSISTENT with the existing project figures, or wants an
  existing figure restyled to match them. Triggers: policy plot, IPoM, IPoM
  Diciembre, BCCh chart, gráfico, incidencias / contribuciones bars, sectoral
  decomposition bars, exposure map, IRF figure, "make it consistent with the
  other plots", "use our colors", terms-of-trade / oil / manufacturing shock
  figures. Applies to Julia (Plots.jl) figures in the NK-IOSOE codebase and to
  quick matplotlib previews.
metadata:
  type: reference
---

# IPoM / BCCh policy-figure house style

This is the chart format used throughout the **Networks-DPMP (NK-IOSOE)** project
and intended to mirror the **Banco Central de Chile IPoM** incidence/contribution
charts (sampled from *IPoM Diciembre 2025, Gráfico I.8*). Use it for any figure
that may end up in a policy note, IPoM box, presentation, or the paper.

When restyling or creating a project figure, **edit the canonical Julia source**
(`Model/julia_dynare/shock_plots_common.jl` and the per-shock `*_shock_analysis.jl`
drivers) rather than producing a one-off. The helpers below already exist there —
prefer reusing them over re-implementing.

## Palette (do not substitute)

```
IPOM_NAVY      = "#1F2E58"   # dark navy   — direct/own channel, single-series bars, all titles
IPOM_LIGHTBLUE = "#4DB3E9"   # light blue  — network / encadenamientos
IPOM_RED       = "#C21E24"   # red         — markup / "Márgenes + Expectativas"
IPOM_ORANGE    = "#ED6331"   # orange      — "Otros costos marginales"
IPOM_GREEN     = "#3BB957"   # green
IPOM_YELLOW    = "#E7EB13"   # yellow
# Totals/aggregate marker: BLACK diamond (markershape=:diamond), as in IPoM incidence charts.
```

Established component orders/labels (already coded in `shock_plots_common.jl`):

- **6-way price decomposition:** `[NAVY, LIGHTBLUE, GREEN, ORANGE, YELLOW, RED]`
- **group3 (Directo / Encadenamientos / Márgenes / Otros):** `[NAVY, LIGHTBLUE, RED, ORANGE]`
- **group2 (Efectos directos / Efectos de encadenamiento):** `[NAVY, LIGHTBLUE]`
- **Single-series bar (e.g. exposure map):** `NAVY`

## Format checklist (every policy figure)

1. **Title:** left-aligned, navy, two lines — name on line 1, **units in parentheses**
   on line 2: `title="$(name)\n($(units))"`, `titlefontsize=12`,
   `titlefontcolor=IPOM_NAVY`, `titlelocation=:left`. No separate y-axis label;
   the units line replaces it.
2. **Legend:** frameless, on top — `legend=:outertop`, `legend_columns=3`,
   `legendfontsize=9`, `foreground_color_legend=nothing`,
   `background_color_legend=nothing`. Single-series charts may still carry one
   legend entry for the units, or set `legend=false`.
3. **Grid:** off — `grid=false`. White background.
4. **Zero line:** solid thin black — `hline!([0.0], color=:black, lw=0.8, label="")`.
   Always keep zero in `ylims` so deflationary/negative parts are never clipped.
5. **Bars:** manual rectangles via `shock_bar_rect(x, y0, y1, 0.65)`, `alpha=0.85`,
   `linecolor=:white`, `linewidth=0.3`. Stack positives upward, negatives downward.
6. **Totals:** black diamond markers, `markersize=7`, `markerstrokewidth=1`.
7. **x-axis:** sector names rotated — `xrotation=55`, `xlims=(0.3, nsec+0.7)`.
   **Drop "Administración Pública"** so the sector axis lines up with the
   `decomposition_pi3_*` figures the chart is read beside
   (`keep = [i for i in 1:nsec if !occursin("blica", bar_names[i])]`).
8. **Size & margins:** `size=(1400, 760)`, `bottom_margin=20mm`, `left_margin=10mm`,
   `top_margin=4mm`. Pad y-limits ~10% below / 14% above the data range.
9. **Language:** Spanish labels throughout. Units phrasings in use:
   `"desviación del EE, puntos porcentuales anualizados"`,
   `"% desv. del EE"`, `"pp del PIB"`, `"Trimestres"` (x-axis on IRF lines).
10. **Line IRFs:** baseline in `IPOM_NAVY`, `lw=2.5`; `xlabel="Trimestres"`,
    `ylabel="% desv. del EE"` (or `"Desv. pp anual del EE"`); dashed black zero line.

## Canonical Julia recipe (bars)

```julia
P = Main.Plots
keep  = [i for i in 1:nsec if !occursin("blica", bar_names[i])]
vals  = data_vec[keep]; names = bar_names[keep]; n = length(keep)
dmax  = maximum(vcat(vals, 0.0)); span = max(dmax, 1e-6)
p = P.plot(
    xticks=(1:n, names), xrotation=55,
    title="$(title_str)\n($(units_str))",
    titlefontsize=12, titlefontcolor=IPOM_NAVY, titlelocation=:left,
    size=(1400, 760), legend=:outertop, legend_columns=3, legendfontsize=9,
    foreground_color_legend=nothing, background_color_legend=nothing,
    grid=false, ylims=(0.0, dmax + 0.14*span),
    bottom_margin=20P.mm, left_margin=10P.mm, top_margin=4P.mm,
    xlims=(0.3, n + 0.7))
P.hline!(p, [0.0], color=:black, lw=0.8, label="")
lbl_used = false
for i in 1:n
    s = shock_bar_rect(i, 0.0, vals[i], 0.65)
    P.plot!(p, s, color=IPOM_NAVY, label=(lbl_used ? "" : series_label),
            alpha=0.85, linecolor=:white, linewidth=0.3)
    lbl_used = true
end
```

For stacked multi-component decompositions and the dual-axis (affected vs. rest)
variant, **reuse `make_ge_decomp_fig(...)`** in `shock_plots_common.jl` instead
of reimplementing; pass colors/labels from `ge_component_style`, `group3_style`,
or `group2_style`.

## Quick matplotlib preview (when Julia isn't available)

For a fast visual check outside the Julia pipeline, mirror the same rules:
navy two-line left title, `loc="lower center", bbox_to_anchor=(0.5,1.0), ncol=3,
frameon=False` legend, hide top/right/left spines, `ax.grid(False)`,
`axhline(0, color="black", lw=0.8)`, bars `color="#1F2E58", alpha=0.85,
edgecolor="white", linewidth=0.3`, `xticklabels(rotation=55, ha="right")`.
A preview is **not** a deliverable — the real figure must be regenerated from the
edited `.jl` source so renderer, fonts, and PDF output stay consistent.

## Guardrails

- Never invent a new palette or relabel established components; reuse the
  constants and `*_style` helpers above.
- Don't drop or reorder sectors except the documented "Administración Pública"
  drop for axis alignment.
- Don't bake numbers into a figure — figures read from the model output
  (`*_shock_analysis.jl` context), never hardcoded.
- Keep everything in Spanish for policy/IPoM outputs.
