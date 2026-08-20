# Atalay (2017) — How Important Are Sectoral Shocks?

**Enghin Atalay**, *American Economic Journal: Macroeconomics*, 9(4): 254–280, 2017.
Working-paper version dated 22 June 2017 (UW-Madison), retrieved from
<https://enghinatalay.github.io/elasticities.pdf>.

**PDF status:** not in `Literatura/literature pdf/` — this note is built from the
full text of the author's working paper. Download the PDF from the link above and
drop it in as `atalay_2017_aej_macro.pdf`; the online appendices B–F (which
contain the plant-level ε_Q estimates and the country-by-country regressions) are
NOT in the working-paper text and must be fetched separately from the AEJ page.

**Why this paper matters to NK-IOSOE:** it is the source of our two production
elasticities, `epsM = 0.1` and `epsY = 0.8`. It is also — contrary to what the
project has assumed — a model **with capital**, and its factor shares are built
on exactly the gross-output denominator we adopted on 2026-08-19. Both facts bear
directly on the capital nest.

---

## 1. Technology

Equation (2), CRS, per industry J:

$$Q_{tJ}=A_{tJ}\Big[(1-\gamma_J)^{1/\epsilon_Q}\big((K_{tJ}/\delta_J)^{\delta_J}(L_{tJ}/(1-\delta_J))^{1-\delta_J}\big)^{\frac{\epsilon_Q-1}{\epsilon_Q}}+\gamma_J^{1/\epsilon_Q}M_{tJ}^{\frac{\epsilon_Q-1}{\epsilon_Q}}\Big]^{\frac{\epsilon_Q}{\epsilon_Q-1}}$$

> "Each industry produces a quantity (Q_tJ) of good J at date t using capital
> (K_tJ), labor (L_tJ), and intermediate inputs (M_tJ) according to the following
> constant-returns-to-scale production function."

Nesting:

| Nest | Inputs | Elasticity |
|---|---|---|
| Outer CES | value-added composite vs. intermediate bundle M | **ε_Q** |
| Inner (value added) | **capital and labour, Cobb–Douglas** | **1**, by functional form |
| Intermediate bundle, eq. (5) | goods from all N sectors | **ε_M** |
| Investment bundle, eq. (4) | goods from all N sectors | ε_X |

$$M_{tJ}=\Big(\sum_{I=1}^{N}(\Gamma^M_{IJ})^{1/\epsilon_M}M_{t,I\to J}^{\frac{\epsilon_M-1}{\epsilon_M}}\Big)^{\frac{\epsilon_M}{\epsilon_M-1}}$$

**Capital is sector-specific and endogenously accumulated**, eq. (3):
`K_{t+1,J} = (1-δ_K) K_tJ + X_tJ`, with an industry-specific investment good
X_tJ that is itself a CES bundle of other industries' output (weights Γ^X). There
is no reallocation of K across sectors. δ_K = 0.10, β = 0.96, annual frequency.

> "To emphasize, the parameters Γ^M_IJ, Γ^X_IJ, γ_J, and δ_J are time invariant."

## 2. Elasticity estimates — Table 1

Estimating equation (13), pooled; instruments are three military-spending demand
shifters from Acemoglu, Akcigit & Kerr (2016).

| | (1) OLS | (2) OLS + yr FE | (3) IV | (4) IV + yr FE |
|---|---|---|---|---|
| **ε_M** | −0.07 (0.04) | −0.13 (0.04) | −0.13 (0.19) | −0.11 (0.20) |
| **ε_Q** | 1.18 (0.06) | 1.27 (0.06) | **0.84 (0.44)** | **0.88 (0.35)** |
| N | 4800 | 4800 | 4592 | 4592 |

> "For ε_Q, the OLS estimates result in an estimate of 1.2-1.3; the IV
> specifications produce estimates closer to 0.8-0.9. In these specifications,
> the standard errors for ε_Q are substantially larger: unit elasticities — as
> used previously in the literature — cannot be rejected."

**The paper's benchmark:**

> "In the following section, I will refer to **ε_M = 0.1, ε_Q = 1, and ε_D = 1**
> as my benchmark set of parameter values."

Footnote 19, the justification, and the origin of 4/5:

> "It is true that there is some weak evidence in favor of ε_Q < 1. However,
> given the large variability of the estimates of ε_Q, I will choose the
> conventional value of 1 for ε_Q for the benchmark parameter configuration, and
> consider a **secondary specification with ε_Q = 4/5** in many of the other
> robustness checks. ... I choose 0.1 as my benchmark value for ε_M as it lies in
> the middle of positive portion of the 90 percent confidence interval."

Other estimates: plant-level (Online App. B) gives ε_Q in **0.4–0.8**; Figure 2's
average slope implies ε_Q ≈ 0.6; the six-country WIOT exercise uses **ε_M = 1/3**.
Rotemberg–Woodford (1996) get 0.7 for manufacturing; Oberfield–Raval (2015)
0.6–0.9.

### ⚠ Consequence for our calibration

- `epsM = 0.1` **is** Atalay's benchmark. Cite it as such. ✔
- `epsY = 0.8` is **not** his benchmark — it is his *secondary specification*
  (footnote 19) and sits on top of his IV point estimates of 0.84/0.88. The paper
  text currently implies 0.8 is "the" Atalay value; it should say "Atalay's
  secondary specification, and close to his IV point estimates of 0.84–0.88 (s.e.
  0.44, 0.35)". Do not claim it is estimated precisely — unity cannot be rejected.
- Atalay's own headline: **ε_Q barely matters, ε_M is everything.**
  > "within the range of elasticities that I have estimated in Section 3,
  > complementarities among intermediate inputs are important for assessing the
  > role of aggregate fluctuations, but the elasticity of substitution between
  > value added and intermediate inputs is not."

  That is a useful defensive citation for us: it says the 0.8-vs-1.0 choice on
  ε_Y is second-order, so a referee attacking it is attacking a margin the
  originating paper says is flat.

## 3. Factor shares — the convention

Appendix A, Table 7:

> "The third through fifth columns of Table 7 give the cost shares of capital,
> labor, and intermediate inputs. These are computed from the BEA GDP by Industry
> dataset. **The intermediate input cost share is computed as the ratio of
> intermediate input expenditures relative to total gross output. The labor share
> is the ratio of labor compensation to total gross output. The remainder defines
> the capital cost share.**"

This is **exactly our post-2026-08-19 convention**, and exactly the construction
behind the new `alpha_K` column:

| | Atalay | NK-IOSOE (Chile 2021) |
|---|---|---|
| intermediate share | intermediates / GO | `alpha + alpha_V` = (dom + imp) / GO |
| labour share | labour compensation / GO | `1 - alpha - alpha_V - alpha_K` ≈ REM/GO |
| capital share | residual | `alpha_K` = EBE / GO |

Footnote 21: *"When ε_Q = 1, the intermediate input cost share and γ_J are equal
to one another."*

**US capital shares from Table 7, for comparison with ours:** range 0.08
(Primary Metals, Textiles) to 0.51 (Utilities, FIRE); Agriculture 0.32,
Mining 0.23, Oil & Gas 0.40. Our Chilean values run 0.117 (Financiera) to 0.704
(Inmobiliario), with Minería 0.576 — **materially more capital-intensive at the
top than the US**, which is exactly why the no-capital assumption bites harder in
Chile. Atalay's intermediate shares average ≈0.55 of GO, i.e. GO/VA ≈ 2.2 against
our 1.92.

Γ^M and Γ^X come from the 1997 IO Table and Capital Flows Table; 35-sector
Jorgenson KLEMS is used for output quantities only, not for shares. He adds
35% to the diagonal of Γ^X for maintenance and repair (McGrattan–Schmitz 1999),
partly for invertibility (fn. 31).

## 4. Headline result

Measure (17): R² = v′Σ_ind v / v′Σ_all v, v = value-added shares.

> "With our benchmark configuration — (ε_D, ε_M, ε_Q) equal to (1, 1/10, 1) —
> **83 percent** of the variation of aggregate output is due to sectoral shocks."

versus **21 percent** at (1,1,1). Table 2 shows the result is driven almost
entirely by ε_M and ε_D, not ε_Q. Robustness (Table 4): 0.82 with 9 industries,
1.00 with the 1972 IO table. Six countries with ε_M = 1/3: at least half of
aggregate volatility from sectoral shocks in five of six (Japan the exception,
0.30).

Important caveat about the *channel* — the elasticities do not change Hulten:

> "For this class of models, the aggregate impact of shocks to an individual
> sector is only a function of the sector's gross output share; to a first-order,
> the elasticities of substitution do not matter. Instead, the elasticities
> matter because they alter the way in which co-movement in fundamental shocks
> map to co-movement in observable data."

## 5. Notes for the NK-IOSOE paper

1. **Atalay has capital and we (baseline) do not.** Citing him for ε_Y while
   omitting capital is a gap a referee can see: his ε_Q is the elasticity between
   a *value-added composite that contains capital* and intermediates, whereas our
   ε_Y is the elasticity between labour alone and intermediates. The capital nest
   closes that gap.
2. **But our nest is flatter than his.** Our four-limb CES imposes
   σ(K,L) = ε_Y = 0.8; Atalay imposes σ(K,L) = 1 (Cobb–Douglas), Baqaee–Farhi
   (2022) use 0.6. Ours sits between them, which is defensible — but say so
   explicitly rather than leaving it implicit. Adding a separate value-added nest
   would match Atalay exactly at the cost of one more elasticity.
3. **His capital accumulates; ours is fixed.** We are at the LPR φ→∞ limit. In
   Atalay capital moves through investment at annual frequency; at our quarterly
   frequency a fixed factor is the more defensible short-run approximation, but
   it is an approximation and belongs in the text.
4. **ε_M = 0.1 is safe; ε_Y = 0.8 needs the footnote-19 citation.**

## Related notes

- [[fgi_jme_2023]] — same CES structure, no capital, hiring costs as the stand-in
- [[econometrica_2023_rubbo_networks_phillips_curves_and_monetary_policy]] — Remark 3, why labour must absorb all VA
- [[luttini_pasten_rubbo_2024]] — Chilean semi-fixed capital, our nest's parent
- [[baqaee-farhi-2022-supply-and-demand-in-disaggregated-keynesian-economies-with-an-application-to-the-covid-19-crisis]] — sector-specific capital, σ(K,L) = 0.6
