# Atalay (2017) — How Important Are Sectoral Shocks?

**Citation (verified):** Atalay, Enghin (2017). "How Important Are Sectoral Shocks?" *American Economic Journal: Macroeconomics* 9(4), 254–280.

**Public copies:** AEA https://www.aeaweb.org/articles?id=10.1257/mac.20160353 ; author PDF https://enghinatalay.github.io/elasticities.pdf ; Census CES-WP-14-31 https://www2.census.gov/ces/wp/2014/CES-WP-14-31.pdf

**PDF status:** not yet in `Literatura/` — download from the links above (binary fetch not available in the drafting environment).

```bibtex
@article{atalay2017sectoral,
  author  = {Atalay, Enghin},
  title   = {How Important Are Sectoral Shocks?},
  journal = {American Economic Journal: Macroeconomics},
  year    = {2017},
  volume  = {9},
  number  = {4},
  pages   = {254--280}
}
```

## Summary

Estimates the **elasticities of substitution** in a multi-sector production network and shows they determine how much of aggregate volatility is sector-specific. Data: BEA IO tables + Jorgenson KLEMS (gross output, intermediates, value added, prices), 30 industries, annual; IV using military-spending shocks. Estimates:
- **ε_m (across intermediate inputs) ≈ 0.1** — point estimates −0.13 to −0.07; goods from different industries are **complements**, not substitutes.
- **ε_Q (value-added vs intermediates) ≈ 0.84** (benchmark 1.0).

Headline: with ε_m = 0.1, **83%** of aggregate output volatility is sector-specific (avg. cross-sector productivity-shock correlation 0.06); imposing ε = 1 instead attributes only **21%** to sector-specific shocks. **Mechanism:** low substitutability means downstream industries **cannot avoid** their suppliers' shocks, so idiosyncratic disturbances propagate through the network and generate comovement that unit-elasticity models misread as "common shocks."

## Relevance to NK-IOSOE (Chile)

Two uses. (1) **Calibration values, not re-estimation.** Estimating the elasticities Atalay-style needs KLEMS intermediate/price time series we do not have for Chile, so we **borrow the values**: `epsM ≈ 0.1` (= ε_m) and `epsY ≈ 0.8` (≈ ε_Q). These coincide with the model's *original* calibration — and are the opposite of where our SMM drifted (`epsM` to its upper bound 0.5, `epsY` to a corner). (2) **Theoretical justification for fixing them.** Atalay's mechanism is exactly our failure in reverse: our estimator raised `epsM` toward substitutability, which — by his 83%→21% result — collapses the sector-specific share and homogenizes sectoral volatility. Fixing the elasticities at complementarity values is what lets sector-specific shocks survive the network and reproduce heterogeneous sectoral volatility with only moderate aggregate volatility.
