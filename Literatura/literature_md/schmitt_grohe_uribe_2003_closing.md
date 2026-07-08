# Schmitt-Grohé & Uribe (2003) — Closing Small Open Economy Models

**Citation (verified):** Schmitt-Grohé, Stephanie, and Martín Uribe (2003). "Closing Small Open Economy Models." *Journal of International Economics* 61(1), 163–185.

**Public draft:** https://www.columbia.edu/~mu2166/closing.pdf (working-paper version, Nov 2001)

**PDF status:** not yet in `literature pdf/` — download published version.

```bibtex
@article{schmittgrohe2003closing,
  author  = {Schmitt-Groh{\'e}, Stephanie and Uribe, Mart{\'i}n},
  title   = {Closing Small Open Economy Models},
  journal = {Journal of International Economics},
  year    = {2003},
  volume  = {61},
  number  = {1},
  pages   = {163--185}
}
```

## Summary (from WP draft)

Compares five stationarity-inducing devices for the incomplete-markets SOE model: (1) endogenous discount factor (Uzawa), (1a) same without internalization, (2) **debt-elastic interest-rate premium**, (3) convex portfolio adjustment costs, (4) complete markets, (5) no closure (unit root). Main finding: all incomplete-markets variants deliver **virtually identical business-cycle dynamics** (second moments and IRFs); complete markets only smooths consumption. Without a closure, the model has a random-walk component in NFA/consumption and infinite unconditional variances.

Their DEIR spec: `p(d) = ψ₂(exp(d − d̄) − 1)` on the **current debt stock in units of goods** — no valuation ratio, no GDP denominator. Calibrated ψ₂ = 0.000742 (tiny; set to match Canadian CA/GDP volatility). The "innocuousness" result is conditional on small ψ₂.

## Relevance to NK-IOSOE

- Citation for our risk-premium closure after the 2026-07-08 change: premium on `Q_ss·Bstar/GDP_ss` (constant-price debt stock) = SGU form up to a scaling constant.
- Justifies footnote that the closure is (near-)innocuous for business-cycle dynamics and that slow NFA/RER convergence is a known property of the model class (their Fig. 1 IRFs visibly do not return to SS within the plotted horizon).
- Note: our earlier XMAS-style spec (premium on contemporaneous `Q·Bstar/GDP`) broke innocuousness because bbar ≈ 1.61 made the valuation feedback first-order.
