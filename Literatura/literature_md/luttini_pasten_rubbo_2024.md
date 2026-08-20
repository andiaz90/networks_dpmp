# Luttini, Pastén & Rubbo (2024) — Measuring the Redistributive Effects of Monetary Policy: An Application to the Chilean Economy

**Emiliano Luttini** (World Bank), **Ernesto Pastén** (BCCh), **Elisa Rubbo**
(Chicago Booth). **This version: January 2024.**
PDF in `Literatura/literature pdf/Luttini_etal_Libro (1).pdf`.

⚠ **Cite the title exactly as above** — "Measuring the Redistributive Effects",
not "Understanding the Heterogeneous Effects", which is how it circulates in
seminar listings and in the BCCh conference programme.

⚠ **This is a short applied companion to Rubbo (2023), not a self-contained
quantitative paper.** Section IV ("Calibration") runs two pages and contains
**not a single numerical parameter value**. No parameter table, no φ, no
elasticities, no N_f, no ownership matrices. Budget your expectations
accordingly: it gives you the algebra and the data construction, nothing else.

**Why this is the single most important paper for our capital extension:** it is
Rubbo's (2023) network framework applied to *Chile*, by *our own institution*, and
it is the only paper in the folder that carries **both** capital and a segmented
labour market. Our new capital nest is its φ→∞ limit.

---

## Capital block (Section A.2, "Investment")

There are N_f capital assets, indexed f. Each is produced by combining a **fixed
endowment** K̄_f with an investment good I_f, eq. (4):

$$K_f=\big[(1+\phi_f)I_f\big]^{\frac{1}{1+\phi_f}}\bar K_f$$

> "For convenience, we assume that the investment component I_f fully depreciates
> from one period to the next, while the endowment component K̄_f never
> depreciates."

Investment is produced CRS from labour, capital and intermediates, eq. (5), and
sold at marginal cost P^I_f to **capital retailers**, who buy the endowments from
agents, combine them with the investment good, and sell capital services to firms
at a competitive rental R_f. Profit maximisation gives the **capital supply
curve**, eq. (6):

$$U_f^{\phi_f}=\frac{R_f}{P_f}\bar K_f,\qquad U_f\equiv\big[(1+\phi_f)I_f\big]^{\frac{1}{1+\phi_f}}$$

> "[U_f] can be interpreted as a measure of capital utilization."

Profits are `φ_f/(1+φ_f) · R_f K̄_f`, investment expenditure `1/(1+φ_f) · R_f K̄_f`.

**φ_f is the one parameter that nests everything.** φ_f → ∞ gives a pure fixed
factor (K_f = K̄_f, no investment margin); φ_f = 0 gives perfectly elastic capital
supply at a constant rental.

### ⚠ They never assign φ_f a value

Section IV calibrates six objects — employment shares, labour shares, IO
linkages, capital shares, consumption shares, price rigidity. φ is not among
them, and neither are the substitution elasticities, N_f, the α^K mapping, K̄_f,
or the ownership matrices Ξ and Z. Nor do they say whether φ is common or
asset-specific. If you need a number, it is not in this paper; the only lead is
Rubbo (2023), "Monetary Non-Neutrality in the Cross-Section", cited without a
link.

**We recover it from the Chilean accounts instead**, using their own eq. (9):
investment expenditure is the fraction 1/(1+φ) of capital income, and both sides
are observed. Capital income Σα_K,i·GO_i = 121,483 (MIP sheet 1) and FBCF =
36,458 (Cuadro 20, col. 7), so **1/(1+φ) = 0.3001, φ = 2.332, ν = 1/φ = 0.4288**.
Implemented in `Data/build_sector_calibration.py`; see
`Model/julia_dynare/docs/capital_nest.md`.

They also never calibrate the investment bundle. We use Cuadro 20's FBCF by
producing activity: Construcción 67.0%, Servicios empresariales 15.4%, Transporte
7.5%, Manufactura 5.9%, and zero for minería, utilities, financiera, personales
and administración pública.

**N_f is never stated either.** Table 1 defines α^K ∈ R^{N×F}, but Figure 8
reports exactly one capital share per industry and the calibration text produces
only an industry-level residual — so the implicit specification is F = N = 111
with α^K diagonal. That is an inference, not a statement in the paper. We use one
asset per sector, F = N = 12.

## Production and shares

All firms z in industry i, eq. (10): `Y_iz = G_i({L_ihz}, {K_ifz}, {X_ijz})` —
CRS, general aggregator, so the K–L elasticity is unrestricted (they specialise
it in the quantitative section). Input subsidies eliminate markups, eq. (13):
`1 − τ*_i = (ε_i − 1)/ε_i` — same device as Rubbo (2023) fn. 6.

Table 1 definitions:

| Object | Definition |
|---|---|
| Labour shares | `α^L ∈ R^{N×H}, α_ih = W_h L_ih /(MC_i Y_i)` |
| **Capital shares** | `α^K ∈ R^{N×F}, α_if = R_f K_if /(MC_i Y_i)` |
| Input-output matrix | `Ω_in = P_n X_in /(MC_i Y_i)` |

Calibration of the capital share — **exactly our construction**:

> "**Industry-level capital shares.** These shares are computed as the difference
> between one and the sum of the industry-level labor and intermediate inputs
> shares described above."

IO linkages from the 2017 Chilean IO matrix, 111 industries, normalised so each
buying industry's shares sum to its cost share of intermediates in national
accounting, exports and imports excluded. Price rigidity from **confidential VAT
electronic invoices, 2015–2022, transaction level** (a data source we do not have
and cannot replicate — worth knowing before a referee asks why we use Pastén,
Schoenle & Weber's US frequencies instead).

## MINING — there is nothing here

`mining`, `minería`, `copper`, `cobre`, `CODELCO`, `mineral`, `extract*`: **zero
occurrences in the body text.** The only "Mining" in the PDF is a hand-placed
x-axis annotation inside the raster images of Figures 5 and 8.

And the reason is structural — they delete the external sector:

> "As the model is a closed economy, we adjust totals by excluding exports and
> imports." (Section IV.A)
> "These shares are normalized such that their sum for a given buying industry
> coincides with the share in costs of intermediate inputs reported in National
> accounting **once exports and imports are excluded**."

Chilean copper is therefore absent from their calibration and mining is treated
exactly like the other 110 industries. **There is no mining treatment to borrow.**

One readable fact from Figure 5: the single most price-flexible industry in the
economy, δ ≈ 0.949 against ≤0.86 for everything else, sits directly under the
"Mining" annotation — almost certainly copper, a world-price commodity, though
the paper never labels it.

## Figure 8 — the 111 capital shares

The text gives only the range:

> "While capital intensity varies significantly across industries (with capital
> shares ranging from **0.1 to 0.8**), the average employer's capital share is
> similar across demographic groups."

Pixel extraction of the 111 bars (600 dpi, ±0.005) gives **min 0.017, max ≥0.800
(bar 98 may be clipped), mean 0.302, median 0.276** — so the "0.1 to 0.8" in the
text understates the low end; several industries are below 0.04.

**Our GO-weighted α_K is 0.295 against their 0.302 mean** — an independent check
that our EBE/GO construction is in the right place for this economy.

Two corrections to the paper: Figure 8's caption says "sorted in ascending order"
but **the bars are not sorted** — they are in natural industry order, pixel-identical
to Figure 5's x-positions. Figure 9 *is* sorted, and runs a very narrow 0.274 to
0.303 across the 50 clusters, which is their actual point.

## Errata to carry forward

| Location | Issue |
|---|---|
| eq. (6) | `P_f` is undefined; the retailer FOC shows it means **P^I_f**, the investment-bundle price. Definition 3 has the same typo, writing `P_f^{*C}`. |
| eqs. (8)–(9) | written with K̄_f; should be **K_f** (they coincide only at U = 1). The budget constraint (2) correctly uses K_f. |
| eq. (22), matrix W | top-right block reads `Z'(I+Φ^K)^{-1}`; since the rebated retailer profit share is φ/(1+φ), it should plausibly be `Z'Φ^K(I+Φ^K)^{-1}`. **Resolve against Rubbo (2023) before implementing.** |
| §III worker clusters | "50 clusters" = 2 genders × 5 quintiles × 5 ages, but only **four** age bins are listed (35–44 missing). |

## Results relevant to us

> "capital shares ranging from 0.1 to 0.8"

which brackets our Chilean 0.117–0.704 well — an independent check that the
`alpha_K` column is in the right range for this economy.

The headline finding on capital:

> "the presence of semi-fixed capital assets amplifies the cross-sectional
> dispersion of employment [responses]"

and, from the conclusions to Figure 13 ("Employment responses with and without
capital assets"):

> "ignoring the presence of semi-fixed capital assets would lead to
> underestimating the cross-sectional range of [employment responses]"

They also report that monetary policy has a larger effect on workers in sectors
with **more flexible prices or larger capital shares** when labour and capital
are strong complements, and that accounting for heterogeneous capital intensity
would *reduce* the cross-sectional range by a stated amount (p. 1, abstract).

## How to use this in the NK-IOSOE paper

1. **Cite it as the justification for the fixed factor**, not Baqaee–Farhi. It is
   the same country, the same IO vintage family, and the same residual
   construction of the capital share. BF's fn. 8 (McKenzie replication) is the
   theoretical licence; LPR is the empirical precedent.
2. **Our nest is their φ→∞ special case — say so.** That framing turns an
   apparent simplification into a nested restriction, and gives us a natural
   robustness dial: estimating or calibrating a finite φ softens the extreme DRS
   our fixed-factor version implies in mining (RTS 0.42) and housing (0.30).
3. **They find capital matters for the CROSS-SECTION of employment responses.**
   That is directly on top of our sectoral-volatility fit problem — it is prior
   evidence that the nest should move the moments we cannot currently match, and
   in the right direction (more cross-sectional dispersion).
4. **Their labour block is richer than ours** (demographic types h, each with its
   own wage rigidity and labour supply elasticity). We have a single labour type
   with sectoral reallocation costs. If a referee asks why not heterogeneous
   labour, this is the paper to point at as complementary rather than competing.

## Related notes

- [[econometrica_2023_rubbo_networks_phillips_curves_and_monetary_policy]] — the parent framework; Remark 3 is why capital is needed at all
- [[atalay_2017_aej_macro]] — sector-specific accumulated capital, σ(K,L)=1
- [[baqaee-farhi-2022-supply-and-demand-in-disaggregated-keynesian-economies-with-an-application-to-the-covid-19-crisis]] — fn. 8, McKenzie replication argument
- [[fgi_jme_2023]] — the no-capital alternative, hiring costs as stand-in
