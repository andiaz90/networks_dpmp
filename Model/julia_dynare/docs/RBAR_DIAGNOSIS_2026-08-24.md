# The rbar miss IS the std(GDP) miss — 2026-08-24

## The corrected target

`compute_data_moments.jl` §9c, `FORCE_RECOMPUTE=1`:

| transformation | rbar | PC1 share |
|---|---|---|
| HP log levels, COVID **in** (the old target) | 0.2210 | 0.44 |
| HP log levels, COVID **out** (**targeted now**) | **0.1238** | 0.30 |
| log growth rates, COVID out (diagnostic) | 0.0645 | 0.20 |

**44% of the old target was the pandemic.** The rbar loop called `complete_cor`
instead of `mcor`, so it was the one sectoral moment never masked. The gap to the
model widens from 1.9x to **3.4x** (model 0.4166 vs data 0.1238), and the growth-rate
reading is lower still, so the direction is not a filter artefact.

For scale: FSW (2011 JPE) report ~0.90 common-factor share for US IP growth across
117 sectors. Chile's 12 broad sectors give 0.20 on the same transformation.

## The finding

Under the aggregation identity for an equicorrelated panel,

    std(GDP)^2 = [ sum_i sig_i^2 (1-r) + r (sum_i sig_i)^2 ] / n^2

Equal weights reproduce the model's actual std(GDP) to 2% (0.02175 predicted vs
0.02131 actual), so the approximation describes the model well.

| | mean sig_i | implied std(GDP) |
|---|---|---|
| DATA, r = 0.124 | 0.0277 | 0.0128 (actual **0.0148**) |
| MODEL, r = 0.417 | 0.0316 | 0.0218 (actual **0.0213**) |
| **MODEL vols, DATA comovement** | 0.0316 | **0.0144** |
| DATA vols, MODEL comovement | 0.0277 | 0.0191 |

**Give the model the data's comovement and its aggregate volatility lands at
0.0144 against the data's 0.0148.** The model is 44% too volatile at the
aggregate (0.0213 vs 0.0148), and essentially all of that is the comovement
error, not the sectoral volatilities — those are only modestly high (9 of 12
above data, mean 0.0316 vs 0.0277).

So rbar is not one moment of 84. It is the single error that also produces the
headline aggregate-volatility miss. It cannot be written off as a limitation.

## What has been ruled out

| hypothesis | test | verdict |
|---|---|---|
| input substitution in the M/V/L nest | epsY sweep 0.7-1.4 | REFUTED — rbar 0.416-0.457 |
| cross-sector substitution | epsM sweep 0.05-0.50 | REFUTED — rbar 0.387-0.448; obj minimised at 0.20 |
| Cobb-Douglas final demand as a unit-loading common factor | PC1 / loading test | REFUTED — loadings spread 0.26-0.97, PC1 only 0.505, Mining loads 0.46 on PC1 but 0.19 on C |
| common AR(1) persistence compresses dispersion | literature search | UNVERIFIED — no supporting paper found |

PC1 loadings instead track the external+copper variance share almost
monotonically (Agriculture 1.8% -> 0.263; Public Admin 24.5% -> 0.969).

## Live leads

1. **Sector-specific shock persistence.** rho_A is common at 0.985 and rho_om
   common at 0.806, while data sectoral autocorrelations run 0.11-0.89. This is
   an independent 14.3% of the objective, so it is worth fixing on its own
   merits; whether it also moves rbar is a testable side bet, NOT a promise.
2. **The external block.** 41% of GDP variance rests on a pinned sigma_pvstar =
   0.054 from a single terms-of-trade AR(1). That is a lot of weight on one
   calibration decision, and the VD gives it rbar = 0.347 alone.
3. Sectoral demand shocks are the only rbar-reducing force (rbar = -0.028) but
   carry just 0.5% of GDP variance, which is why lambda_om = 2.19 cannot fix
   this on its own.

## Immediate consequence

Every stored estimate was fit against rbar = 0.221. The target is now 0.124 with
weight 3.0, so the objective is not comparable across the change and the whole
profile must be re-run.
