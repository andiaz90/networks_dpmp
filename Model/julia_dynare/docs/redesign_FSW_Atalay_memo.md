# Estimation Redesign — Consistency with FSW (2011) / Atalay (2017) and Adequacy for the NK-IOSOE (Chile) Model

*Working memo, 2026-07-22. Companion to the literature wrappers `foerster_sarte_watson_2011_jpe.md` and `atalay_2017_aej_macro.md`.*

## 1. Why the current SMM fails (recap)

Estimating 36 parameters — including 12 sectoral TFP variances, 12 sectoral demand variances, and the two production elasticities — on 60 second moments is under-identified and mis-targeted. Symptoms from the 1.5h run: near-uniform sectoral output/labor volatility (~0.30) versus the data's heterogeneous 0.01–0.16; rank correlations ≈ 0.2 vs target 1; `corr(Y_i,PH_i)` positive (demand-driven) where the data are negative (supply-driven); `corr(GDP,Q)=0.99`; nine parameters pinned at bounds; a rank-deficient inference Jacobian (the SE = 0 / t = Inf artifact); J-test rejection. The diagnosis: the cycle is driven by a large common/external shock that swamps the sector-specific shocks, and the estimator was free to raise `epsM` toward substitutability — which homogenizes the network — because the objective weights volatility *magnitudes* (95%) over comovement/dispersion.

## 2. The redesign (three moves)

1. **Fix the production elasticities at complementarity values** — `epsM ≈ 0.10`, `epsY ≈ 0.80` — rather than estimating them. *(Implemented: `utils.jl` bounds pinned.)*
2. **Measure the sector-specific shock variances from the data** and fix them, instead of SMM-estimating `isigma_tfp_i` / `sigma_om_i`. *(Pending: a decomposition script — §5.)*
3. **Shrink the SMM to the transmission block** — labor adjustment cost, wage rigidity `kappaw`, and the external/premium parameters — estimated on dispersion- and comovement-focused moments, not 60 magnitude targets. *(Pending: θ-reduction + reweighting.)*

## 3. Consistency with the literature

**Ferrante–Graves–Iacoviello (2023).** They *feed in* measured sectoral productivity shocks and calibrate the reallocation shock to an observed expenditure shift, and SMM-estimate only three transmission parameters (hiring cost, σ_m, σ_Y) on 12 dispersion moments. Our moves 2–3 are the stochastic analog of exactly this division of labor: measure the shocks, estimate only transmission.

**Atalay (2017).** He estimates ε_m ≈ 0.1 (complementarity) and ε_Q ≈ 0.84, and shows the elasticity *is* the answer: with ε_m = 0.1, 83% of aggregate volatility is sector-specific; at ε = 1, only 21%. Our move 1 fixes `epsM`, `epsY` at exactly these values. This is directly consistent — and it is the precise antidote to our failure, which was the estimator pushing `epsM` to 0.5 (i.e. toward the 21% world that homogenizes sectors).

**Foerster–Sarte–Watson (2011).** Their structural factor decomposition of sectoral output + IO into common vs sector-specific components is the template for move 2 with the data we actually have.

**Where we deviate — stated honestly.** Atalay *estimates* the elasticities from KLEMS with a military-spending IV; we *borrow* his point estimates as calibration because Chile lacks the intermediate-input/price time series his method needs. FGI feed a specific historical episode; we instead use an FSW-style stochastic decomposition. So we are consistent with the *family* of approaches — measure or borrow the shocks and elasticities, estimate only transmission — while adapting the specific technique to our data. That adaptation should be flagged in the paper's estimation section, not hidden.

## 4. Adequacy of these references for our model

**What transfers cleanly.** (a) The variance-decomposition logic (sector-specific vs common) is model-agnostic. (b) Atalay's elasticities are production-network technology parameters and the model's production nest is the same CES-over-intermediates form, so the values are portable. (c) FSW's measurement needs only output + IO — inputs we have. (d) The complementarity → heterogeneity mechanism is general.

**What does not transfer — adequacy limits to respect.**

- **External block.** FSW/Atalay are *closed economies*; their "common factor" is a domestic aggregate. Our dominant driver is the import-price / real-exchange-rate (external) block, which they have no counterpart for. Implication: an FSW decomposition on Chilean output will fold external shocks into the "common" factor. That is fine for *isolating the sector-specific idiosyncratic component* (all we need to calibrate `isigma_tfp_i`), but we must not read the "common" share as domestic aggregate TFP.
- **Commodity SOE.** Atalay's 83%-sector-specific result is a *US* number. For Chile — mining-heavy, terms-of-trade-exposed — external/terms-of-trade shocks plausibly account for a large share, so we should **not target or expect 83% sector-specific**. Borrow the *method and the elasticities*, not the quantitative decomposition result.
- **Nominal side.** Both papers are real. They say nothing about our price/inflation moments, the sectoral Phillips curves, or `corr(Y_i,PH_i)`. Those belong to the NK-network references already in the folder (Rubbo 2023, Baqaee–Farhi 2022, Pastén et al. 2020). Net: FSW/Atalay discipline the **quantity / shock-composition** side of the fit; the **nominal / comovement** side needs the NK-network literature. Use them for what they cover.
- **Frequency.** Atalay/FSW are annual; our model and data are quarterly. The elasticity values are technology parameters and are frequency-robust, but the shock-variance measurement (move 2) must be done at quarterly frequency on our series.

**Verdict.** Adequate as the methodological backbone for moves 1–2 (elasticity calibration and data-measured sectoral shocks), with the explicit caveat that the SOE/external and nominal dimensions of our model lie outside their scope and are governed by the NK-network references. They are the right citations for *why we stop estimating sectoral shock sizes and fix the elasticities*; they are not a complete template for the whole model.

## 5. Data feasibility — and a refinement that uses what we have

We have quarterly sectoral **output** (`pib_sectorial_bc`), **employment** (`count_workers_by_sector`), and **deflators** (`deflactor_pib`), 2006–2023, plus the 2021 IO table. We do **not** have sectoral **capital** or **intermediate-input time series**, so a direct Solow-residual sectoral TFP (à la Atalay/FGI) is **not** feasible.

But we have something FSW's output-only method does not exploit: **sectoral prices alongside quantities.** The model has both a supply (TFP) and a demand (taste) shock per sector, and they move prices oppositely — a supply shock raises `Y_i` and lowers `PH_i`; a demand shock raises both. So the sector-level joint behavior of `(Y_i, PH_i)` **identifies the supply/demand mix per sector**, which is exactly the object that is currently mis-estimated (our `corr(Y_i,PH_i)` has the wrong sign). Concretely, move 2 becomes: from the HP-filtered sectoral `(Y_i, PH_i)`, back out sector-specific supply and demand components (via the model's own two-shock mapping, or a sign-restricted decomposition), and calibrate `isigma_tfp_i` and `sigma_om_i` to reproduce them. This is more feasible *and* more disciplined than an output-only FSW factor model, and it targets the comovement failure directly.

## 6. Status and next steps

- **Done:** `epsM`, `epsY` fixed at Atalay complementarity values (`utils.jl`). Literature wrappers written (`foerster_sarte_watson_2011_jpe.md`, `atalay_2017_aej_macro.md`) and indexed.
- **Next (needs a script run on the Chilean data — I can write it, you run it):** the sectoral `(Y_i, PH_i)` supply/demand decomposition (move 2) → outputs a `sectoral_shock_calibration.csv` with `isigma_tfp_i`, `sigma_om_i`.
- **Then:** reduce the estimated θ from 36 to the transmission block (drop the 24 sectoral shock variances and the two elasticities), and reweight the objective so comovement/dispersion moments bind (§3 of the earlier inference audit). This also dissolves the corner-solution / rank-deficient-Jacobian problem, so the standard errors become meaningful.
- **SOE caveat to carry:** keep the external/premium block estimated (or measured), since Chile's external shocks are first-order and outside the FSW/Atalay scope.
