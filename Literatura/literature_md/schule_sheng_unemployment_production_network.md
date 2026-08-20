# Schüle & Sheng — Unemployment in a Production Network

**Finn Schüle, Haoyu Sheng.** Working paper, December 2023 (page last modified
November 2024). <https://www.finnschule.com/papers/labor_network/>
PDF: <https://www.finnschule.com/papers/labor_network/schule_sheng_upn.pdf>

**PDF status:** not yet in `Literatura/literature pdf/` — download from the link
above and save as `schule_sheng_upn.pdf`. This note is built from the abstract
and the paper page only; the body has not been read.

**Why it is in this folder:** it is the closest existing answer to "a production
network model built around the labour market". Everything else in our literature
either has frictionless labour (Rubbo, Comin et al., Ghassibe, Cesà-Bianchi &
Ferrero) or a reduced-form reallocation cost (FGI, and us).

---

## What they do

Abstract, verbatim:

> "We model a production network economy with sectoral and occupational
> unemployment by incorporating matching between job seekers across various
> occupations and employers in different production sectors. We derive
> expressions that unpack how the impact of microeconomic shocks on output and
> unemployment depends on the interaction between the network linkages, search
> costs, and changes in labor market tightness. When labor markets are slack, our
> model predicts larger output and employment responses because the
> network-adjusted labor productivity gain outweighs search costs. Calibrating
> our model to the U.S. economy, we demonstrate that our model significantly
> amplifies the response of aggregate output and unemployment to productivity
> shocks in any sector and changes the relative importance of sectors to
> aggregate output and unemployment. Our model nearly doubles the output response
> compared to an efficient production network and triples the unemployment
> response compared to a multi-sector search model following a productivity shock
> to durable manufacturing."

Structure: DMP-style matching between **occupations** (job seekers) and
**sectors** (employers) — a bipartite labour network layered on the production
network. State-dependence is central: responses are larger when the labour market
is slack.

## Why it matters for NK-IOSOE

1. **It is the natural extension of our labour block.** We have a single
   representative household, undifferentiated labour, and quadratic reallocation
   costs (Ψ_l aggregate, ψ_l sectoral). They replace that with genuine search
   frictions and get amplification that a reallocation cost cannot produce —
   notably *state dependence*, which our linear-in-levels solution cannot
   generate at all.
2. **The amplification numbers are a benchmark.** "Nearly doubles" the output
   response relative to an efficient network and "triples" the unemployment
   response. If we ever want to argue our labour adjustment cost is sufficient,
   this is the paper that says it is not.
3. **Sector rankings change.** They find the labour block "changes the relative
   importance of sectors to aggregate output and unemployment" — i.e. it moves
   exactly the cross-sectional object our rank-correlation moments (moments 7–9)
   target. Relevant to the sectoral-volatility fit problem: our current model
   asks TFP shocks and reallocation costs alone to reproduce the sectoral
   ranking, and a labour-market channel is a mechanism we have switched off.
4. **Positioning.** It is a US working paper with no SOE dimension, no imported
   inputs and no exchange rate. Our contribution — networks + SOE + nominal
   rigidity heterogeneity — does not overlap. Cite it as the frontier on the
   labour side and as the natural next step, not as a competitor.

## Caveats before citing

- Working paper, not published as of August 2026; check for a newer version and
  for a journal placement before it goes in the bibliography.
- The site's own citation block gives year 2023 while the page metadata says the
  paper was modified November 2024. Verify the version you cite.
- I have not read the body: the abstract does not say whether wages are
  Nash-bargained or rigid, nor whether prices are sticky. Both matter for whether
  the mechanism is complementary to ours or substitutes for our Phillips-curve
  heterogeneity. **Read before citing in the text.**

## Related notes

- [[fgi_jme_2023]] — convex hiring costs, the reduced-form version of this
- [[cardoza_grigoli_pierri_ruane_2025_restud]] — the empirical counterpart: workers really do move along the network
- [[baqaee-farhi-2022-supply-and-demand-in-disaggregated-keynesian-economies-with-an-application-to-the-covid-19-crisis]] — sectoral labour with downward nominal wage rigidity and factor unemployment
- `xmas.pdf` (García, Guarda, Kirchner & Tranamil 2019, BCCh WP 833) — full DMP block with endogenous separation, no network; the in-house donor if we ever add unemployment
