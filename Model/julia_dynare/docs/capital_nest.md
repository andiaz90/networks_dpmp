# Capital in the NK-IOSOE (branch `capital-nest`)

Added 2026-08-20. **There is one version of the model.** Capital follows
Luttini, Pastén & Rubbo (2024): semi-fixed sector-specific assets with supply
elasticity ν = 1/φ. No environment switches, no defaults, no alternative
branches — every parameter is read from `Data/*.csv` and the run fails loudly if
a file is missing.

```bash
cd Model/julia_dynare
julia --project=. main_SOE_gap.jl          # solve + IRFs + figures + tables
julia --project=. oil_shock_analysis.jl    # (also agr/min/mfg/agrmin)
julia --project=. run_smm_estimation.jl    # SMM (slow; run main first)
```

First time only: `julia --project=. -e 'import Pkg; Pkg.instantiate()'`, then
**commit the generated `Manifest.toml`** — unpinned dependency versions are what
caused the cluster-only zero-decision-rule failure in job 7198.

The calibration CSVs below are already built and committed. Rebuild them only
when the IO workbooks or the convention change:
`python3 Data/build_sector_calibration.py` (add `--dry-run` to preview).


| what | where it comes from |
|---|---|
| α_K,i = EBE_i/GO_i | `Data/sector_calibration.csv`, col `alpha_K` |
| χ_I (investment bundle) | `Data/sector_calibration.csv`, col `chi_I` |
| ν = 1/φ = 0.4288 | `Data/capital_calibration.csv`, row `nu_K` |
| K̄_i (endowment) | solved in calibration, written to `params_jl.mod` |

## Where the output goes

| path | what |
|---|---|
| `mod/dynare_irfs.csv` | **raw IRF data** — every variable × every shock × 40 quarters (~83 MB) |
| `mod/dynare_ss.csv`, `dynare_endogenous_variance.csv` | steady state, unconditional variances |
| `figures/irfs/` | aggregate and sectoral IRF panels, one per shock |
| `figures/exercises/` | output-gap and shock-comparison plots |
| **`figures/oil_shock/`** | **oil-shock decompositions — the paper's headline exercise** |
| `figures/{agr,min,mfg,agrmin}_shock/` | the other sectoral shock exercises |
| `tables/` | `moment_fit_baseline.txt`, `steady_state_baseline.csv`, `gdp_va_shares_baseline.csv`, `vae_baseline.csv`, the `.tex` decomposition tables |
| `estimation_results/` | SMM θ, checkpoint, fit table, standard errors |

Since 2026-08-20 `main_SOE_gap.jl` runs `oil_shock_analysis.jl` itself, as a
subprocess at the end, so `figures/oil_shock/` is always in step with the model
just built. `RUN_OIL=0` skips it. A failure there is reported but does not abort —
everything else is already on disk by that point.

## The problem it solves

With three CES limbs (materials, imported inputs, labour) and constant returns,
labour must absorb the entirety of value added. Rubbo (2023, Econometrica),
Remark 3:

> "With constant returns to scale and labor being the only factor of production,
> labor must account for the entirety of value added."

So the baseline implies a labour share of value added of **1.000** against
**0.416** in the Chilean accounts. Excedente bruto de explotación is 56.8% of
Chilean VA, and it is concentrated exactly where it is least defensible to call
it profit: mining is 24.8% of national EBE, real estate a further 13.4%.

## The fix

A fourth CES limb: sector-specific capital services `K_i = Kbar_i * U_i` with
weight `alpha_K_i = EBE_i / GO_i`.

$$Y^i=A^i\Big[\alpha_{mi}^{1/\epsilon_Y}M^{\frac{\epsilon_Y-1}{\epsilon_Y}}+\alpha_{vi}^{1/\epsilon_Y}V^{\cdot}+\alpha_{ki}^{1/\epsilon_Y}(\bar K_i U_i)^{\cdot}+(1-\alpha_{mi}-\alpha_{vi}-\alpha_{ki})^{1/\epsilon_Y}L^{\cdot}\Big]^{\frac{\epsilon_Y}{\epsilon_Y-1}}$$

What it does **not** touch: `alpha_i`, Γ, the Leontief multiplier, GO/VA = 1.92.
It only reallocates value added from labour to capital.

| | no-capital model | this model | data |
|---|---|---|---|
| labour share of VA | 1.000 | **0.432** | 0.416 |
| GO-weighted α_K | 0 | 0.295 | EBE/GO = 0.295 |
| investment/VA | 0 | **0.17** | FBCF/VA = 0.1704 |
| GO/VA | 1.92 | 1.92 | 1.92 |

The residual 1.6pp is net taxes on production, which the model has no home for.

## Design decisions

**The endowment never enters goods-market clearing.** K̄ is not produced, so no
sector supplies it. What *does* enter is the investment demand it generates:
`chiI_i * EInv / PH_i`, with `EInv = ν/(1+ν) · Σ R_i K_i`.

**No household budget constraint for the rebate.** The retailer profit
`ν/(1+ν)·R·K` closes by Walras, exactly like Rotemberg profits — the model has no
explicit budget constraint and does not need one.

**Two solver modes** (`steady_ntwsoe_system.jl`) — not two models, just the
difference between calibrating and estimating:

- *Calibration* (`main_SOE_gap.jl`): 5×nsec system, K solved jointly under the
  unit normalisation `K_i = alpha_Ki * Y_i`. This makes the steady-state capital
  cost share equal `alpha_Ki` exactly — with A = 1 the FOC gives R_i = MC_i, so
  `alpha_Ki (MC_i/R_i)^(ε_Y−1) = alpha_Ki`. The solved K is written to the .mod
  as the endowment `Kbar_i`.
- *Estimation* (`smm_model_moments.jl`): 4×nsec system with `Kbar` held **fixed**
  at its calibrated value. The endowment is a datum, not a function of the
  estimated parameters, so the capital cost share drifts as the steady state
  moves — as in Baqaee–Farhi and Luttini–Pastén–Rubbo, where K̄_f is fixed.

Verified: at the calibrated parameters the estimation-mode solve reproduces the
calibration-mode solve to 2.1e-14.

## The LPR specification

Luttini, Pastén & Rubbo (2024) do **not** hold capital fixed. Asset f is produced
from a fixed endowment K̄_f and an investment good I_f, their eq. (4):

$$K_f=\big[(1+\phi_f)I_f\big]^{\frac{1}{1+\phi_f}}\bar K_f\quad\Longrightarrow\quad U_f^{\phi_f}=\frac{R_f}{P^I_f}\bar K_f,\qquad U_f\equiv K_f/\bar K_f$$

so 1/φ is the capital supply elasticity. Two accounting identities follow (their
eqs. 8–9): investment expenditure is `R_f K_f/(1+φ)` and the complement
`R_f K_f φ/(1+φ)` is retailer profit rebated to households.

**We parameterise by ν = 1/φ, not φ**, because φ = ∞ (the pure fixed factor) is
not representable while ν = 0 is. The supply curve is written in deviation form,
`U_i = ((RK_i/RKss_i)(PIinv_ss/PI_inv))^ν`, so U = 1 at the steady state and the
level is absorbed into the units of K̄.

### ⚠ Two things LPR do not give you

**φ is not in the paper.** Their Section IV calibrates six objects — employment
shares, labour shares, IO linkages, capital shares, consumption shares, price
rigidity — and φ is not among them. Neither are the substitution elasticities,
the number of capital assets N_f, the α^K mapping, or the ownership matrices Ξ
and Z. There is no parameter table anywhere in the paper.

**But their own identity pins it down from Chilean data.** Investment expenditure
is the fraction 1/(1+φ) of capital income, and both sides are observed:

| | | source |
|---|---:|---|
| capital income Σ α_K,i·GO_i | 121,483 | MIP sheet 1 |
| FBCF | 36,458 | Cuadro 20, col. 7 |
| ⇒ 1/(1+φ) | **0.3001** | |
| ⇒ φ | **2.332** | |
| ⇒ **ν = 1/φ** | **0.4288** | |

`build_sector_calibration.py` computes and prints this on every run.

**The investment bundle is also ours, not theirs.** LPR never calibrate their
investment aggregator. We use a Cobb–Douglas bundle with weights read off
Cuadro 20's FBCF by producing activity: Construcción 67.0%, Servicios
empresariales 15.4%, Transporte 7.5%, Manufactura 5.9%, Comercio 2.3%, Agro 1.5%,
Inmobiliario 0.5%; **zero for minería, utilities, financiera, personales and
administración pública**. This makes construction the marginal buyer of capital
goods, which is economically right and is a channel the baseline model does not
have at all.

## Mining — and why LPR are no help here

**LPR never mention mining, copper or CODELCO.** The words appear zero times in
the body; "Mining" survives only as a hand-placed annotation on two figures. More
decisively, they strip the external sector out entirely: *"As the model is a
closed economy, we adjust totals by excluding exports and imports."* Chilean
copper is therefore absent from their calibration, and there is no treatment to
borrow.

Their capital shares are still a useful external check: mean **0.302** across
111 industries (min 0.017, max ≥0.80) against our GO-weighted **0.295**.

**What actually matters for us is a feature of our own model.** Mining output is
already exogenous — `Y_2 = Y2_ss*exp(A_2)` — so the fixed factor does **not** damp
mining quantities. My earlier warning that mining would be over-damped was wrong.
What α_K,2 = 0.576 does instead is move MC_2, hence PH_2 through the NKPC, hence
the copper rent `(PH_2 − MC_2)·Y_2` that `phi_cu` repatriates.

**That is the interaction to watch.** Capital reassigns most of mining's
operating surplus from an unmodelled residual to an explicit factor rent, while
the copper block continues to treat `(PH_2 − MC_2)·Y_2` as the rent that leaks
abroad. With `SMM_PHI_CU = 0` (the default) there is no double count. **Before
setting `SMM_PHI_CU > 0`, check that the two are not both claiming the same
surplus.**

## ⚠ The thing to watch: returns to scale

With K̄ fixed, **returns to scale in the variable factors equal 1 − α_K exactly**:

| sector | α_K | RTS | | sector | α_K | RTS |
|---|---:|---:|---|---|---:|---:|
| Agro | 0.299 | 0.701 | | Transporte | 0.289 | 0.711 |
| **Minería** | **0.576** | **0.424** | | Financiera | 0.117 | 0.883 |
| Manufactura | 0.172 | 0.828 | | **Inmobiliario** | **0.704** | **0.296** |
| Utilities | 0.346 | 0.654 | | Empresariales | 0.379 | 0.621 |
| Construcción | 0.206 | 0.794 | | Personales | 0.164 | 0.836 |
| Comercio | 0.222 | 0.778 | | Adm. pública | 0.135 | 0.865 |

Housing at 0.30 and mining at 0.42 are **severe** decreasing returns — far
steeper than the mild DRS usual in this literature.

**ν = 0.4288 softens this substantially.** Because capital responds to its
rental, the *effective* returns to scale rise (verified numerically, 2% input
scaling):

| | mining | housing | manufactura | mean |
|---|---:|---:|---:|---:|
| φ → ∞ (pure fixed factor, = 1−α_K) | 0.423 | 0.295 | 0.828 | 0.699 |
| **ν = 0.4288 (this model)** | **0.530** | **0.391** | **0.881** | **0.771** |

Still a large, heterogeneous departure from constant returns, and it bears
directly on the sectoral output-volatility fit problem. **Every SMM estimate
stored before 2026-08-20 predates it and must be re-run.**

## Lineage

| Paper | Device |
|---|---|
| Baqaee & Farhi (2022) fn. 8 | McKenzie (1959) replication argument — the theoretical licence |
| **Luttini, Pastén & Rubbo (2024)** | semi-fixed capital assets, Chile — closest template, our φ→∞ parent |
| Comin, Johnson & Jones (2023) | hard output capacity constraint instead |
| Imbs, Jondeau & Pelgrin (2011) | DRS straight in the labour exponent |
| Atalay (2017) | sector-specific *accumulated* capital, Cobb–Douglas with labour |
| FGI (2023) | **no** capital; convex hiring costs stand in for it |

Note our flat four-limb CES imposes σ(K,L) = ε_Y = 0.8. Atalay imposes 1
(Cobb–Douglas inside a value-added nest), Baqaee–Farhi 0.6. Ours sits between
them; matching Atalay exactly would require a separate value-added nest and one
more elasticity.

## Files touched

| File | Change |
|---|---|
| `Data/build_sector_calibration.py` | writes `alpha_K` (= EBE/GO) and `chi_I` (FBCF bundle); writes `capital_calibration.csv` with ν and φ |
| `mod/NK_SOE_lev_gap2.mod` | params `alphaK_i`, `Kbar_i`, `chiI_i`, `RKss_i`, `nuK`, `PIinv_ss`; vars `U_i`, `RK_i`, `PI_inv`, `EInv` (+flex); fourth CES limb; labour weight net of α_K; capital FOC and supply curve; investment in market clearing; `GDP = C + EInv + TB` |
| `steady_ntwsoe_system.jl` | capital is now required, not optional; investment demand; calibration vs estimation modes |
| `main_SOE_gap.jl` | reads α_K, χ_I, ν from the CSVs and errors if absent; 5×nsec SS solve; `Kbar_ss`, `RKss_ss`, `PIinv_ss_val`, `EInv_ss`; SS identity is `C + EInv + TB ≡ ΣVA`; prints EInv/ΣVA and realised factor shares |
| `utils.jl` | writes the new params (no defaults) |
| `smm_estimation.jl` / `smm_model_moments.jl` | carry α_K, K̄, χ_I, ν; estimation-mode SS solve with K̄ fixed |
| `shock_plots_common.jl` | 7th bar "Capital" (α_K·r̂K) in the marginal-cost decomposition; `group3_components` folds it into "Otros"; new `IPOM_PURPLE` |
| `oil/agr/min/mfg/agrmin_shock_analysis.jl` | `alpha_L_vec` nets out α_K (**was a live bug** — the labour bar was overstated and the residual hid it); own SS solves extended; `GDP = C + EInv + TB` |
| `figs_SOE_gap.jl` | investment column in the aggregate IRF table; ΔU and ΔRK in the sectoral table; investment and factor-share rows in the SS summary |
| `run_model.sh` | the runner |

## What has been verified

Julia 1.10.4 (no packages), running the **real** `steady_ntwsoe_system.jl`
against the **real** calibration CSVs:

| check | result |
|---|---|
| all 12 edited `.jl` files parse | clean |
| calibration-mode SS residual | 4.4e-15 |
| estimation-mode SS residual (K̄ fixed) | 2.8e-14 |
| estimation reproduces calibration | 2.1e-14 |
| R_ss = MC_ss, K = α_K·Y | 2.2e-16 |
| EInv/VA | 0.1665 (data 0.1704; toy prices, not the model SS) |

**Dynare has still never parsed the `.mod`.** The steady-state block is
verified, the model block is not. `julia --project=. main_SOE_gap.jl` is the
first time it will be, and it is the step most likely to fail: investment is a
final-demand block the model did not have before, so if Dynare's residuals are
non-zero suspect the market-clearing term `chiI_i*EInv/PH_i` or the GDP identity
before suspecting the capital algebra.

Checks printed by that run:

| line | expect |
|---|---|
| `[SS-ACCT] max|goods-mkt clearing resid|` | ~1e-12 |
| `[SS-ACCT] investment expenditure EInv` | EInv/ΣVA ≈ 0.17 (data 0.1704) |
| `[SS-ACCT] ΣVA vs GDP = C+EInv+TB` | gap ~0 |
| `[SS-ACCT] factor shares of ΣVA` | labour ≈0.43, capital ≈0.57 |
| `capital SS check: max\|RKss - MC_ss\|` | ~1e-12 |
