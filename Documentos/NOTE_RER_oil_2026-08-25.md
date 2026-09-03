# The real exchange rate and the oil shock — diagnostic note

**Date:** 2026-08-25 · **Model:** `NK_SOE_lev_gap2.mod` · **θ:** current baseline

---

## 1. How the exchange rate is defined

The CPI is the numeraire (`p_g = exp(om_g)*C/C_g` at l. 278 together with the Cobb–Douglas
aggregator normalises the consumption-basket price to 1), so every price in the model —
`PH_i`, `PV`, `PO`, `w` — is CPI-deflated and the exchange-rate variable is the
**CPI-based real exchange rate**

$$Q_t \equiv \mathcal{E}_t P^*_t / P_t, \qquad Q\uparrow \;=\; \text{depreciation}.$$

There is **no nominal exchange rate level** in the model. Only its growth rate is recovered,
residually, at l. 891: `pi_e = Q/Q(-1)*pi/Pipstar`.

Sign convention verified from three equations:

| Equation | Line | Implication |
|---|---|---|
| `PO = Q*POstar` | 978 | world oil price in domestic real terms rises with `Q` |
| `PH_2 = Q*Pcstar` | 732 | copper likewise |
| `X = omegaX*(PX/Q)^(-etastar)*Ystar` | 918 | exports rise with `Q` |

`Q` is pinned by three equations: the Euler/UIP condition (l. 874), the Schmitt-Grohé–Uribe
debt-elastic premium `r_star = Rworld*exp(-phi_b*(bbar - Q_ss*Bstar/GDP_ss))` (l. 885), and the
BOP accumulation `Q*Bstar = -TB + r_star(-1)*Q*Bstar(-1)/Pipstar` (l. 896).

### The Backus–Smith structure

Substituting `pi_e` into the Euler equation, domestic inflation cancels and what remains is a
*foreign-real* Euler equation:

$$\frac{Q_{t+1}}{Q_t}=\frac{\Pi^*}{\beta R^*_t}\left(\frac{Z^c_{t+1}}{Z^c_t}\right)^{\gamma}$$

With `phi_b = 0.0024`, $R^*_t \approx R^w = \Pi^*/\beta$, the prefactor is ≈ 1 and the model's
RER is a near-exact **Backus–Smith relation**: `Q` tracks the path of domestic marginal utility.
Any shock that lowers consumption on impact and lets it recover must produce a real
**appreciation** on impact. The trade-balance / net-foreign-asset channel is present but, at
this `phi_b`, second-order.

---

## 2. What the model delivers for a 10% oil shock

Current run (`tables/oil_shock_irfs.csv`, 2026-08-25, ρ = 0.9):

| Quarter | Q | GDP | π | C | TB |
|---|---|---|---|---|---|
| 1 | **−0.211** | −0.194 | +0.036 | −0.155 | −2.550 |
| 2 | −0.215 | −0.191 | +0.027 | −0.164 | −1.874 |
| 3 | −0.198 | −0.182 | +0.017 | −0.164 | −1.482 |

i.e. a **real appreciation of ~0.21%**, exactly as the Backus–Smith algebra predicts:
consumption falls and recovers, so `Q` must undershoot and climb back.

### Sensitivity to the SOE closure

`phi_b` is **not estimated**. It is hard-coded at `phi_b_val = 0.0024` (`main_SOE_gap.jl` l. 513,
"XMAS posterior mean"; `mod/params_jl.mod` l. 22) and appears nowhere in the 38-element θ vector
(`CSV_PARAM_NAMES`, `utils.jl` l. 534).

Impact response of `Q` across the existing `PHIB_OVERRIDE` runs:

| φ_b | Q on impact | vintage |
|---|---|---|
| 0.001 | −0.496 | 2026-06-02 |
| **0.0024 (baseline)** | **−0.211** | 2026-08-25 |
| 0.010 | −0.031 | 2026-07-08 |
| 0.050 | +0.093 | 2026-07-08 |

The sign flips at roughly **φ_b ≈ 0.015–0.02** — six to eight times the calibrated value and
about twenty times the Schmitt-Grohé–Uribe (2003) value. The closure cannot be nudged into a
depreciation within any defensible parameter range.

> *Caveat:* the 0.001/0.01/0.05 runs predate the copper sector, the capital nest, the government
> block, GHH preferences and sticky wages. They are indicative of the direction, not the level.
> The sweep could not be re-run here — Julia and Dynare are not installed in this environment.

---

## 3. Is the appreciation counterfactual? — apparently **not**

The textbook prior for Chile (net oil importer, terms-of-trade deterioration) is a
depreciation. The data do not say that.

HP-filtered (λ = 1600) quarterly logs, BIS RBCL real broad REER (sign-flipped to the model's
convention: up = depreciation) against Brent, with the real relative price of copper in Chile
(BCCh copper-mining deflator ÷ GDP deflator) as a control:

**Full sample 2006Q1–2023Q4, n = 72**

| Specification | β(oil) | β(copper) | R² | 10% oil ⇒ |
|---|---|---|---|---|
| Q ~ oil | −0.0824 (0.0190) | — | 0.212 | **−0.82% Q** |
| Q ~ copper | — | −0.0685 (0.0291) | 0.073 | |
| Q ~ oil + copper | −0.0747 (0.0204) | −0.0291 (0.0289) | 0.224 | **−0.75% Q** |

**Excluding COVID (2020Q1–2021Q2), n = 66**

| Specification | β(oil) | β(copper) | R² | 10% oil ⇒ |
|---|---|---|---|---|
| Q ~ oil | −0.0686 (0.0232) | — | 0.120 | **−0.69% Q** |
| Q ~ oil + copper | −0.0733 (0.0269) | +0.0113 (0.0326) | 0.122 | **−0.73% Q** |

corr(oil, copper) = +0.37 (full) / +0.50 (ex-COVID), so the two commodity prices do co-move —
but controlling for copper leaves the oil coefficient essentially unchanged and significant
(t ≈ 3.7 full, 2.7 ex-COVID).

**Conclusion: the model has the right sign and undershoots the magnitude by a factor of 3–4**
(−0.21% vs −0.7 to −0.8%). Raising `phi_b` would move the model *away* from the data, not
toward it.

### Two caveats before this goes in the paper

1. This is a correlation with the **unconditional** Brent price, not with an identified
   oil-supply shock. The likely global driver is the well-documented negative oil–dollar
   correlation: a broad USD depreciation mechanically appreciates the CLP in a *broad* REER
   index. Settling this properly needs Kilian or Baumeister–Hamilton oil-supply shocks.
2. The copper control is the Chilean mining deflator relative to the GDP deflator, which itself
   embeds the exchange rate. It is a same-source proxy, not a clean world copper price.

---

## 4. What the SMM does and does not see

Current fit of the RER block (`tables/moment_fit_baseline.txt`, moments 41–43):

| # | Moment | Data | Model | Weight | % of objective |
|---|---|---|---|---|---|
| 41 | std(Q) | 0.0429 | 0.0367 | 1/d² ≈ 543 | 0.3% |
| 42 | autocorr(Q) | 0.7168 | 0.4018 | 2.0 | 2.9% |
| 43 | **corr(GDP,Q)** | **+0.0289** | **−0.1146** | **1.0 (default)** | **0.3%** |

Three observations.

- **std(Q) is not the binding constraint.** It is nearly matched, and `phi_b` is not estimated
  anyway, so nothing in the objective is "holding `phi_b` down" — it was never free.
- **`w[43]` is never assigned.** `build_weighting_matrix` sets weights for 37–42, 44–46, 59–65,
  78, 83–84, but index 43 keeps the `ones(N_MOMENTS)` default. This is the **Backus–Smith
  moment** — the single statistic that measures the mechanism in §1 — and it carries the
  smallest weight in the aggregate block, by omission rather than by decision.
- **The model has the Backus–Smith sign and the data do not.** Data corr(GDP,Q) ≈ 0
  (acyclical); model −0.11 (appreciates in booms). The miss is real but contributes 0.3% of the
  objective, so the estimator has no reason to fix it.

- **autocorr(Q) is badly missed** (0.72 vs 0.40): the model's RER is far too transitory,
  which is the same near-complete-markets pathology seen from a different angle — `Q` inherits
  the persistence of consumption rather than that of a slow-moving net-foreign-asset state.

---

## 5. Suggested next steps

1. **Set `w[43]` deliberately.** Either give it the 3.0 used for the other bounded correlations
   (`corr(GDP,pi)`, `rbar`) or document explicitly that the Backus–Smith moment is untargeted
   and reported as validation. Leaving it at an unassigned default is the one option that is
   hard to defend to a referee.
2. **Add a country-premium shock,** as the code comment at l. 869 already anticipates:
   `exp(-phi_b*(bbar - ...) + nu_t)`. This is the standard device (Itskhoki–Mukhin 2021) for
   breaking Backus–Smith without distorting the debt dynamics, and it is the natural route to
   both a larger oil-shock RER response and a higher autocorr(Q). It must be named a country
   risk-premium shock, not a preference shock.
3. **Re-run the φ_b sweep on the current model** (0.0024 / 0.005 / 0.01 / 0.02 / 0.05) so the
   sensitivity table in §2 is a single vintage. As a robustness exhibit only — the data in §3
   say this is not the margin to move.
4. **Get an identified oil-supply shock series** before claiming the §3 elasticity as a target.

---

*Sources: `Model/julia_dynare/mod/NK_SOE_lev_gap2.mod`, `main_SOE_gap.jl`, `utils.jl`,
`smm_estimation.jl`, `tables/oil_shock_irfs.csv`, `tables/oil_chiib0p0{1,5}_shock_irfs.csv`,
`tables/moment_fit_baseline.txt`, `Data/raw/reer_chile_bis.xlsx`,
`Data/raw/brent_crude_quarterly.csv`, `Data/deflactor_pib.csv`.*
