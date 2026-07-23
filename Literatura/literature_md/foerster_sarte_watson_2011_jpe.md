# Foerster, Sarte & Watson (2011) — Sectoral versus Aggregate Shocks

**Citation (verified):** Foerster, Andrew T., Pierre-Daniel G. Sarte, and Mark W. Watson (2011). "Sectoral versus Aggregate Shocks: A Structural Factor Analysis of Industrial Production." *Journal of Political Economy* 119(1), 1–38.

**Public copies:** JPE https://www.journals.uchicago.edu/doi/abs/10.1086/659311 ; FRB Richmond WP 08-07 https://ideas.repec.org/p/fip/fedrwp/08-07.html ; non-technical summary (Sarte, Richmond Fed EQ 2011) https://www.richmondfed.org/-/media/richmondfedorg/publications/research/economic_quarterly/2011/q2/pdf/sarte.pdf

**PDF status:** not yet in `Literatura/` — download from the JPE/Richmond Fed links above (binary fetch not available in the drafting environment).

```bibtex
@article{foerster2011sectoral,
  author  = {Foerster, Andrew T. and Sarte, Pierre-Daniel G. and Watson, Mark W.},
  title   = {Sectoral versus Aggregate Shocks: A Structural Factor Analysis of Industrial Production},
  journal = {Journal of Political Economy},
  year    = {2011},
  volume  = {119},
  number  = {1},
  pages   = {1--38}
}
```

## Summary

Decomposes U.S. industrial production growth into **aggregate (common) factors** and **sector-specific idiosyncratic shocks** using a *structural* factor model that embeds the input-output structure. Up to **117 sectors** — the finest disaggregation at which BEA IO tables can be matched to sectoral output. Two stages: (1) extract common factors from sectoral output growth; (2) use the IO matrix to attribute the residual cross-sectional comovement to network propagation, isolating the truly sector-specific disturbances. Data requirement is light: **sectoral output growth series + an IO matrix** — no KLEMS capital/intermediate time series needed for the core decomposition.

Headline: over **1984–2007**, sector-specific disturbances explain about **half** of the variation in U.S. IP, and their importance rose over time. Evolving IO linkages did **not** by themselves raise shock propagation.

## Relevance to NK-IOSOE (Chile)

The **data-feasible template** for our redesign. We have exactly its inputs — sectoral real output (`pib_sectorial_bc`, 12 sectors, 2006–2023 quarterly) and the 2021 IO table (`IO_2021_chile` / `mip_12x12`) — and we lack the KLEMS series Atalay needs. FSW's factor decomposition lets us **measure the sector-specific component** of each sector's output volatility from the data and **calibrate `isigma_tfp_i` (and `rho_A`) to it**, rather than SMM-estimating 12 shock variances that the network + common shocks otherwise homogenize and swamp (the source of our uniform ~0.3 sectoral vol, rank-corr ≈ 0.2, and corner solutions). Feed the shocks; let the model transmit them.
