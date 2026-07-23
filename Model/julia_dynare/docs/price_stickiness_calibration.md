# Sectoral Price-Stickiness Calibration (κ_i): Pastén → Chile, and the Frequency Units

*Working documentation, 2026-07. Covers the full pipeline from Pastén, Schoenle &
Weber (2020) US price-change frequencies to the model's sectoral Rotemberg
adjustment costs κ_i, the US→Chile rescaling, and a frequency-unit inconsistency
in the current code.*

## 1. The pipeline

The sector-specific Rotemberg price-adjustment cost κ_i is built in four conceptual
steps. Steps 1–2 are baked into the data files; steps 3–4 happen in
`main_SOE_gap.jl`.

| Step | Object | Where |
|------|--------|-------|
| 1 | Pastén et al. (2020) **US monthly** frequency of price adjustment, aggregated to the 12 sectors | `Data/raw/fpa_vector_few_industries_chile.csv` |
| 2 | Rescaled to **Chilean monthly** frequency (US→Chile level adjustment) | `Data/fpa_vector_few_industries_chile.csv` |
| 3 | Time-aggregate **monthly → quarterly**: `θ_q = 1 − (1 − θ_m)^3` | `main_SOE_gap.jl` (toggle `SMM_FPA_QUARTERLY`) |
| 4 | Calvo → Rotemberg: `κ_i = θ̄(ε−1)/[(1−θ̄)(1−θ̄β)]`, `θ̄ = 1 − θ_q` | `main_SOE_gap.jl:~247` |

`θ` = frequency of price adjustment (fraction of firms resetting per period);
`θ̄ = 1 − θ` = Calvo probability of NOT adjusting; higher κ_i = stickier.

## 2. Step 1 — Pastén et al. (2020) US frequencies

Source: Pastén, Schoenle & Weber (2020), *JME* 116, "The Propagation of Monetary
Policy Shocks in a Heterogeneous Production Economy." Frequencies come from
confidential **BLS producer-price (PPI) microdata, 2005–2011** (~25,000
establishments; services included from 2005). Reported US facts: **mean monthly
frequency 0.19** (mean duration 4.6 months), ranging from 0.0401 (semiconductor
mfg, 24.4 months) to 0.9375 (dairy, 0.36 months). These are **monthly**.

Note: the PPI is a *producer*-price survey and does **not** cover government /
public administration, so Pastén et al. have **no public-administration
frequency**. Sector 12's value therefore does not descend from a Pastén
public-admin number — it comes from whatever BEA-industry aggregation produced the
raw vector, and should be treated with caution.

## 3. Step 2 — US → Chile rescaling

`Data/fpa_vector_few_industries_chile.csv` = the raw US vector rescaled toward
Chilean price stickiness. Empirically the transform is close to a **constant
scaling of durations by ≈1.44×** (equivalently a ≈0.73 multiplier on
frequencies):

| statistic | value |
|-----------|-------|
| frequency ratio (CL/US) | mean 0.729, cv 1.9% |
| duration ratio (CL/US)  | mean 1.443, cv 1.1% |
| US mean frequency (this aggregation) | 0.275 |
| Chile mean frequency | 0.202 |

So Chilean prices are calibrated as ~44% stickier (longer-lived) than the US
Pastén cross-section, holding the cross-sectional *pattern* fixed. The Chilean
mean frequency 0.20 remains **monthly** in magnitude (a quarterly frequency would
average ≈0.47 = 1−(1−0.19)^3).

**Chilean anchor (confirmed):** the Chilean manufacturing price-change frequency
comes from **Albagli, Grigoli, Luttini, Quevedo & Rojas (2026), "Beyond Costs: The
Dominant Role of Strategic Complementarities in Pricing," Central Bank of Chile**
(draft June 10 2026), Appendix B. They compute the frequency of price adjustment
from Chilean firm-to-firm electronic-invoice microdata (2014–2023) — a
PRODUCER-price basis, consistent with Pastén's US PPI basis. Reported monthly
frequencies (Tables B.1–B.2), all MONTHLY:

| measure | median freq | duration |
|---------|-------------|----------|
| Observed, all firms, sales-weighted | 40.5% | 2.0 mo |
| Regular (sales filtered), all firms, sales-weighted | 31.6% | 2.7 mo |
| Observed, PPI products | 30.9% | 2.7 mo |
| Regular, PPI products | 25.2% | 3.5 mo |

The model's Manufacturing value θ_m = 0.366 (duration 2.19 mo) matches the OBSERVED,
all-firms headline (~40% / "roughly every two months"). The US Pastén Manufacturing
frequency (raw) is 0.488 (1.49 mo); the ratio 0.366/0.488 ≈ 0.75 is applied to the
whole cross-section (only manufacturing is observed for Chile, so the Chile/US
stickiness ratio is assumed constant across sectors — this is why we rescale).

> **Two choices to pin down for the paper:**
> 1. **Observed vs regular.** θ_m = 0.366 uses the OBSERVED frequency (includes
>    temporary sales). Calvo/Rotemberg nominal rigidity should arguably use the
>    REGULAR (sales-filtered) frequency, and Pastén's US PPI is itself closer to
>    regular. The regular Chilean anchor is 31.6% (all firms) or 25.2% (PPI),
>    i.e. STICKIER (2.7–3.5 mo) — using it would raise all κ_i.
> 2. **All-firms vs PPI.** For a manufacturing/producer anchor consistent with the
>    US PPI source, the PPI figure (obs 30.9% / reg 25.2%) is the cleaner match
>    than the all-firms 40%.

## 4. Steps 3–4 — the frequency-unit inconsistency (and fix)

The model is **quarterly** (β = 0.986 quarterly; data HP-filtered at λ = 1600,
2006Q1–2023Q4). The Chilean frequencies from step 2 are **monthly**. The current
code (`main_SOE_gap.jl`) feeds them into the Calvo→Rotemberg map **without the
monthly→quarterly aggregation** (step 3 omitted), i.e. it treats a monthly
frequency as if it were quarterly. This overstates stickiness by roughly the
monthly/quarterly gap: a sector that adjusts prices every 4.6 *months* is encoded
as adjusting every 4.6 *quarters*.

Effect on κ_i and implied durations (ε = 10, β = 0.986):

| sector | θ (monthly) | κ legacy (θ as quarterly) | dur | κ fixed (θ_q) | dur |
|--------|-------------|---------------------------|-----|---------------|-----|
| Agriculture   | 0.410 | 31   | 1.9 q | 2.9  | 0.6 q |
| Manufacturing | 0.366 | 42   | 2.2 q | 4.1  | 0.7 q |
| Personal Serv.| 0.081 | 1088 | 11.8 q| 133  | 3.9 q |
| Public Admin. | 0.066 | 1612 | 14.6 q| 202  | 4.9 q |

Under the fix, durations run 0.6–4.9 quarters (≈2–15 months) — inside Pastén's US
range and preserving the sector ordering — instead of the legacy 2–15 *quarters*.

Corroboration: the legacy over-stickiness matches the model's misfit — sectoral
price volatilities are too *low* (model std(PH) ≈ 0.02 vs data 0.03–0.17) and
output volatilities too *high*, both symptomatic of prices that are too sticky.

**Fix (implemented):** `main_SOE_gap.jl` applies `θ_q = 1 − (1 − θ_m)^3` when
`SMM_FPA_QUARTERLY=1` (default OFF, preserving legacy behavior). Setting it to 1
makes the calibration internally consistent with the quarterly model. This is a
paper-wide change (all κ_i move; all IRFs and shock-decomposition results shift).

## 5. What the paper currently says (gap)

`paper.tex` documents only that "the price adjustment cost κ_i is taken from
\citet{pasten2021sectoral}," and the "Frequencies of price adjustment" subsection
is red-flagged ("[necesitamos usar microdatos… tenemos??]", "[cómo hacen esta
agregación?]"). The paper does **not** document (a) that the values are US Pastén
frequencies **rescaled to Chile** (~1.44× duration), nor (b) the monthly→quarterly
handling. Both belong in the calibration section; §3–§4 above are drafted to fill
those gaps once the rescaling target/source (§3 OPEN) is confirmed.
