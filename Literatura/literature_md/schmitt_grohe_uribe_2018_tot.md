# Schmitt-Grohé & Uribe (2018) — How Important Are Terms-of-Trade Shocks?

**Citation (verified):** Schmitt-Grohé, Stephanie, and Martín Uribe (2018). "How Important Are Terms-of-Trade Shocks?" *International Economic Review* 59(1), 85-111. (NBER WP 21253, 2015.)

**Public copies:** IER https://onlinelibrary.wiley.com/doi/abs/10.1111/iere.12263 ; NBER https://www.nber.org/papers/w21253 ; author page http://www.columbia.edu/~mu2166/

**PDF status:** not yet in `Literatura/` — download from the NBER/author links above (binary fetch not available in the drafting environment).

```bibtex
@article{schmittgrohe2018terms,
  author  = {Schmitt-Groh{\'e}, Stephanie and Uribe, Mart{\'i}n},
  title   = {How Important Are Terms-of-Trade Shocks?},
  journal = {International Economic Review},
  year    = {2018},
  volume  = {59},
  number  = {1},
  pages   = {85--111}
}
```

## Summary

Measures the terms of trade from **observed** country price data and estimates its business-cycle contribution via country-by-country SVARs across 38 poor and emerging economies. Headline: terms-of-trade shocks explain **less than 10%** of the variance of aggregate activity empirically. A calibrated/estimated three-good SOE model (importables, exportables, nontradables) predicts a larger share (~30%), so theory **over-assigns** the terms of trade relative to the data — the **"terms-of-trade disconnect."** The lesson: when the external price process is disciplined by its own observed series rather than by domestic moments, its role is modest, not dominant.

## Relevance to NK-IOSOE (Chile)

The external (import-price) shock `eps_pvstar` is currently **SMM-estimated (theta[31], theta[32]) from domestic moments only** — std(Q), autocorr(Q), corr(GDP,Q), std(TB/GDP) — with **no observed import-price / terms-of-trade series** anchoring it. It inflates to an unconditional std around 0.54 and drives `corr(GDP,Q)=0.99` (the disconnect in extreme form: the external block explains essentially everything). SGU (2018) is the reference for the fix: **measure the external process from the observed terms of trade (for Chile: the copper price) and fix rho_pvstar/sigma_pvstar out of the estimated theta**, so its contribution is data-disciplined and modest. Complements Garcia-Cicco-Pancrazi-Uribe (2010, country premium from observed macro) and the BCCh XMAS practice of multiple external shocks each anchored to an observable (copper, world rate, foreign demand).
