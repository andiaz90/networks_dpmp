# Estimating sectoral price rigidities κ_i from the data — design and identification

**Project:** NK-IOSOE (Networks-DPMP) · **Date:** 2026-09-02
**Status:** IMPLEMENTED. See §10 for what was built, what the data said, and what changed.

---

## 0. Why this is worth doing

The κ_i currently in the model are **not Chilean**. `Data/build_fpa_vector_chile.py` takes
Pastén–Schoenle–Weber (2020) **US PPI** monthly frequencies, and rescales all twelve durations
by the single constant factor that maps US manufacturing duration onto the one Chilean number
that is observed — Albagli, Grigoli, Luttini, Quevedo & Rojas (2026) regular-PPI manufacturing
median frequency, 0.252/month. So the *level* is anchored to Chile in one sector and the entire
*cross-sectional pattern* — the object the paper's heterogeneity claim rests on — is American.

That is the referee's first question, and the paper has no answer to it today.

Calibrated vector, for reference (θ_q = 1−(1−θ_m)³, κ = s(ε−1)/[(1−s)(1−sβ)], s = 1−θ_q):

| sector | θ_m | θ_q | duration (q) | κ_i | log κ_i |
|---|---|---|---|---|---|
| Agriculture | 0.286 | 0.637 | 1.57 | 8.0 | 2.08 |
| Mining | 0.179 | 0.446 | 2.24 | 24.6 | — (no NKPC) |
| Manufacturing | 0.252 | 0.581 | 1.72 | 11.0 | 2.40 |
| Utilities | 0.155 | 0.396 | 2.53 | 33.9 | 3.52 |
| Construction | 0.080 | 0.221 | 4.52 | 136.5 | 4.92 |
| Trade | 0.194 | 0.477 | 2.10 | 20.4 | 3.02 |
| Transport | 0.110 | 0.295 | 3.38 | 70.3 | 4.25 |
| Financial | 0.095 | 0.259 | 3.87 | 95.8 | 4.56 |
| Real estate | 0.070 | 0.196 | 5.11 | 178.8 | 5.19 |
| Business serv. | 0.101 | 0.272 | 3.67 | 85.1 | 4.44 |
| Personal serv. | 0.051 | 0.145 | 6.89 | 337.6 | 5.82 |
| Public admin. | 0.041 | 0.119 | 8.41 | 508.6 | 6.23 |

mean log κ (excl. mining) = 4.22 (κ ≈ 68), sd = 1.34.

**Mining has no Phillips curve** (`NK_SOE_lev_gap2.mod` `@#if i == 2`: PH_2 = Q·Pcstar, law of one
price). So there are **11** free κ_i at most, not 12.

---

## 1. What κ_i actually is in this model — and a normalization warning

Around π_ss = 1 the Rotemberg adjustment cost vanishes in steady state, so κ_i is a **pure
dynamics parameter**. Linearized, the sectoral NKPC is

  π^H_{i,t} = (ε/κ_i)·m̂c^r_{i,t} + β E_t π^H_{i,t+1}

Two consequences that shape the whole design:

1. **κ_i does not touch the steady state.** Verified: `mod/dynare_ss.csv` has `pi = 1.0`.
   So unlike ε_Y, freeing κ triggers **no** `recompute_ss_cached!` — the SS cache hits 100% and
   an evaluation costs what it costs today. Freeing κ is *cheap*; freeing ε_Y was not.
   → Do **not** add κ to the `need_ss` test in `smm_estimation.jl`.

2. **ε and κ_i are not separately identified from dynamics.** With `MARKUP_SUBSIDY = true`
   (LPR eq. 13) the steady-state markup is exactly 1, so ε has been stripped of its one
   non-dynamic role and now enters *only* as the slope ε/κ_i. Keep ε = 10 pinned and say in the
   paper that what is estimated is the **NKPC slope** λ_i ≡ ε/κ_i, reported as an implied Calvo
   duration. Estimating ε and κ jointly would be a rank failure, not a hard identification problem.

---

## 2. Are the current 84 moments enough? — **No.**

Not for 11 free κ_i. Enough in *count*, not in *content*. Four reasons, in order of severity.

### 2.1 There is not a single sectoral price-persistence moment

The moment vector has `autocorr(GDP)`, `autocorr(pi)`, `autocorr(Q)` and twelve `autocorr(Y_i)`
(indices 64–77) — and **zero** autocorrelations of sectoral prices or sectoral inflation.

Price stickiness *is* a persistence parameter. Identifying it from unconditional variances alone
repeats exactly the mistake diagnosed on 2026-08-21 for ρ_A: "until today this vector held exactly
ONE autocorrelation against three free persistence parameters… which is why rho_A ran to its 0.99
bound." κ_i would be identified only through how the HP filter reweights variance across
frequencies. This gap alone is disqualifying.

### 2.2 The sector-level price moments already have three other jobs

Sector-level price information in the current set is exactly 24 numbers: `std(PH_i)` (13–24) and
`corr(Y_i,PH_i)` (47–58). Both are already load-bearing for:

- **λ_A, λ_om** (θ[37:38]) — the supply/demand mix, whose stated identification is "std(Y_i) pins
  the level, corr(Y_i,PH_i) pins the mix";
- **ε_Y** (θ[2]) and the pinned ε_M — network propagation;
- **cl** and **κ_w** — both damp the price/quantity split through marginal cost.

Adding 11 κ_i puts a fourth mechanism onto the same 24 moments. Expect a ridge, not an estimate.

### 2.3 Worse: κ̂_i would absorb the shock-measurement error

The 24 sectoral shock sizes are pinned at data-measured values, with only two global scales free.
That *helps* identify κ in principle (std(PH_i) becomes informative in levels rather than only as
a ratio) — but it is fragile in exactly the way the project already knows about: the measured
σ_i vector is `compute_sectoral_shocks.jl`'s **raw first pass**, with PC1 removal and
model-inversion rescaling never applied. Any error in σ_i lands in κ̂_i one-for-one, because
std(PH_i) ≈ f(σ_i)/g(κ_i). Estimated stickiness would then be a re-labelling of a measurement
residual — indefensible in a paper whose contribution *is* the stickiness heterogeneity.

Mitigation: use **ratio** and **conditional** moments (§3), which are invariant to the σ_i scale.

### 2.4 Quarterly deflators are a blunt thermometer

The calibrated quarterly durations run **1.57 to 8.41 quarters**. Agriculture, Manufacturing and
Trade reset prices roughly **twice within a quarter**. A quarterly GDP deflator cannot distinguish
"flexible" from "very flexible" — the time aggregation has already destroyed the variation. So the
attainable object from macro data is the **sticky tail plus a level**, not twelve sharp κ's.
Expect strongly asymmetric identification: tight below, near-flat above. **Estimate log κ_i**, and
expect wide upper confidence intervals for services.

### 2.5 Precedent inside this project

`log(κ_V)` — the *import-price* Rotemberg cost, the same object one layer over — was pinned on
2026-08-21 with |dobj| = 1.8e-4 against 1.0e+1 for the λ's: flat to four orders of magnitude, not
identified by this moment set. And the sector-specific-ρ experiment (2026-08-24) is the direct
cautionary tale: freeing measured cross-sectional heterogeneity made the fit **worse** in two
parameterizations and drove λ_ρ hard onto its bound.

Design so that a negative result is a **publishable nested test**, not a failed cluster run.

---

## 3. Moments to add

Ordered by information-per-unit-of-work. Blocks A and B use data already in `Data/`.

### A. Sectoral inflation persistence — 11 moments (do this first)

  `autocorr(π^H_i)`, i ≠ 2, first-order, from `deflactor_pib_sa.csv`, HP-filtered, pandemic
  quarters dropped (same convention as `autocorr_Q` and blocks 64–77).

**Why it identifies κ.** Higher κ ⇒ the price level partially adjusts toward its flexible-price
target ⇒ positively serially correlated inflation and a persistent relative price.

**How it separates κ from shock persistence.** ρ_{A,i} raises persistence of *both* quantity and
price; κ_i raises price persistence *relative to* quantity persistence. The pair
{`autocorr(Y_i)` (already in, 66–77), `autocorr(π_i)` (new)} is what separates them. Write this
argument into the paper — it is the identification paragraph a referee will look for.

Cost: ~30 lines in `compute_data_moments.jl`, mirrored in `smm_model_moments.jl`. No new data.

### B. Price/quantity volatility ratio — 11 moments

  `std(π^H_i) / std(Δ log Y_i)` (growth rates, not HP-filtered levels).

Monotone decreasing in κ_i and — the point — **invariant to the σ_i scale**, so it survives §2.3.
Use growth rates: the HP-filtered price *level* volatility already in the set (13–24) is dominated
by low frequencies and is a weak κ signal; the inflation-rate volatility is the sharp one.

### C. Conditional pass-through to a common exogenous cost shock — 11–22 moments (highest value)

This is the moment class with genuinely exogenous variation, and the one a referee will demand.

Data: local projection, h = 0…8,

  Δ log P^H_{i,t+h} = a_i + β_i(h)·Δ log(NER·P*)_t (or Brent, or ToT) + controls + u

with `reer_chile_bis.csv`, `brent_crude_quarterly.csv`, `tot_chile_quarterly.csv` — all present.

Model counterpart: run the **identical** local projection on data simulated from (T, R), so the
estimator is the same object in both. (An analytic IRF ratio is cheaper but is not the same
statistic, and the difference matters at short horizons.)

**Use the shape, not the level.** β_i(0) also depends on α_{V,i} (imported-input intensity), which
is calibrated from the IO table, not estimated. The *timing* statistic — β_i(4)/β_i(0), or the
half-life of the cumulative response — nets that out and is stickiness almost purely. Recommend
one timing moment per sector, and optionally β_i(0) as a second, down-weighted.

This is also the moment that will discipline the sticky tail where §2.4 says variances cannot.

### D. Cross-sector inflation comovement — 1 moment

  `rbar_π` = average pairwise corr(π_i, π_j), the price analogue of the existing `rbar` (78).

Heterogeneous κ plus the network jointly determine it (Rubbo 2023). One number, near-zero cost,
and it speaks to the paper's aggregate-Phillips-curve section.

### E. Relative-price dispersion — 1 moment (optional)

  std over time of the cross-sectional sd of sectoral relative prices. Welfare-relevant in the
  Baqaee–Farhi / Rubbo sense; cheap; supports the counterfactual section.

**Something must come down in the weighting matrix.** Project rule, from the 2026-08-21 pass:
this is a *replacement*, not an addition. Candidate: cut the weight on `std(PH_i)` levels
(13–24) in favour of A + C, which carry the same information in a form that is not confounded
with σ_i.

---

## 4. Parameterization — do **not** free 11 κ_i in one go

A ladder, each rung a nested restriction of the next. This mirrors the λ_A/λ_om and λ_ρ designs
already in the codebase and makes every possible outcome reportable.

**S0 — global level sweep (already implemented, 0 new parameters).**
`SMM_KAPPA_SCALE` exists (2026-08-25) and applies κ_i = s·κ_i^cal at build time with a read-back
assertion. Profile s over a grid *before* writing any new code. It answers "is the LEVEL right?"
for the price of a few runs. Note the existing diagnostic comment already warns that errors go
both ways (Construction and Business Services are ~2× too volatile), which is itself evidence
that a level shift alone will not do it — and therefore an argument for S1.

**S1 — level + dispersion, 2 free parameters (recommended headline).**

  log κ_i(θ) = θ₃₉ + θ₄₀·(log κ_i^cal − mean_j log κ_j^cal),  i ≠ 2

- θ₃₉ = log κ̄, the overall stickiness level (Chilean anchor);
- θ₄₀ = λ_κ, the dispersion scale.

Two nested tests, both about the paper's actual contribution:

| restriction | meaning | what it tests |
|---|---|---|
| θ₄₀ = 0 | homogeneous stickiness | does κ heterogeneity matter at all? |
| θ₄₀ = 1, θ₃₉ = mean log κ^cal | the current calibration | do Chilean macro data agree with the US cross-section? |

Suggested bounds: θ₃₉ ∈ [log 2, log 800] = [0.69, 6.68]; θ₄₀ ∈ [0, 2.5].

**S2 — group-wise κ, 3–4 parameters.** Goods vs services, or terciles of κ^cal. Run only if S1's
λ̂_κ is interior and precisely estimated.

**S3 — 11 free log κ_i (robustness table only).** With the augmented moment set *and* a shrinkage
penalty toward the micro-calibrated vector, `obj + ω‖log κ − log κ^cal‖²`, ω reported. Never the
headline: 11 free parameters against quarterly deflators will produce estimates the data cannot
support (§2.4).

---

## 5. The alternative that may dominate all of this: measure κ_i directly

Frequency of price adjustment is **observable**, and this is a central bank.

- `Data/diccionario-de-datos-bbdd-ipp-anonimizadas.xlsx` is the data dictionary for the
  **anonymized IPP microdata** — the microdata itself is not in the repo. Request access.
- Albagli, Grigoli, Luttini, Quevedo & Rojas (2026) already computed Chilean regular-PPI
  frequencies; the manufacturing median (0.252) is the anchor the calibration currently uses.
  **They may already have the sectoral breakdown.** Luttini is a co-author on LPR, which this
  model already uses for the markup subsidy — this is one email, and it could replace §4 entirely.
- IPC microdata (INE) would cover the service sectors that PPI does not.

If a Chilean sectoral frequency vector can be obtained, the paper's structure improves sharply:
κ_i becomes a **calibrated, externally measured** object (Table 2, from data), the S1 estimation
becomes a **validation exercise** (does the macro-implied λ_κ agree with the micro-measured
cross-section?), and the referee question in §0 disappears. **Start the access request in
parallel with §6 — it has the longest lead time and the highest payoff.**

---

## 6. Implementation checklist

Ordered. Steps 1–3 are the ones with known failure modes in this codebase.

1. **`utils.jl`**
   - `const KAPPA_CAL` — the calibrated vector, single source of truth, read the same way
     `RHO_A_MEASURED` is. Do not recompute it in two places.
   - Append θ[39] = `log_kappa_bar`, θ[40] = `lambda_kappa` to `PARAM_LABELS`, `CSV_PARAM_NAMES`,
     `LB`, `UB`. `N_THETA` 38 → 40. **Append, never insert** — the block slices in the
     decomposition printers are hardcoded.
   - `FREE_THETA` += [39, 40].
   - Append new moment blocks to `MOMENT_NAMES` / `MOMENT_BLOCKS`, again at the **end**.
     Update `N_MOMENTS`.
   - `build_weighting_matrix`: weight the new blocks, and cut something (§3).

2. **`smm_model_moments.jl` — `_KLEIN_STRUCT_IDX = collect(1:40)`. Mandatory.**
   κ_i enters the pricing-equation Jacobian directly. Leaving the cache key at 1:38 reproduces
   the stale-decision-rule bug exactly (see `networks-dpmp-klein-cache-stale-R`): the objective
   would be silently flat in κ and the run would "converge" having never moved it.

3. **`smm_estimation.jl`, the per-evaluation param block (~line 794)**
   - Compute κ_i(θ) and `set_param!(context, "kappa_$(i)", κ_i)` for i ≠ 2, in the same loop as
     `cl_$(i)` and `rho_tfp1_$(i)`.
   - Guard: κ_i finite and > 0 ⇒ else `return NAN58, false`.
   - **Read back** `pvec("kappa_$(i)")` once at startup and assert it equals what was written.
     The 2026-08-25 κ sweep returned byte-identical moments at 0.1× and 4× because a "safety"
     reload wiped every `set_param!`. Reuse that read-back pattern verbatim.
   - Do **not** add κ to `need_ss` (§1).
   - `--threads=1` for any κ run, for the same per-thread `_SMM_PARAMS` reason as
     `SMM_KAPPA_SCALE`.

4. **`compute_data_moments.jl`** — add blocks A–E as new *columns* in `sectoral_moments.csv`;
   never overwrite existing columns. Keep the SA source (`*_sa.csv`), the 2006Q1–2023Q4 window
   and the pandemic exclusion identical to the existing blocks — the sample-mismatch bug
   (`networks-dpmp-sectoral-sample-mismatch`) came from exactly this.

5. **Pre-flight sensitivity scan** (`probe_free_dims.jl`) on θ[39], θ[40] **before** any cluster
   time. Accept only if |dobj| is within ~one order of magnitude of cl (4.5e-1). If it comes back
   at the log κ_V level (1e-4), stop and report flatness — that is a result, and it is cheap.

6. **Profile before searching.** Grid λ_κ ∈ [0, 2] at fixed θ₃₉, then grid θ₃₉. The objective
   profile is a figure worth putting in the paper and shows whether κ is identified at all —
   the same lesson as ε_Y. Only then run joint CMA-ES on the cluster.

7. **Inference.** `smm_inference.jl` sandwich; report s.e. on log κ̄ and λ_κ, and the implied
   Calvo durations with delta-method CIs. Report the two nested tests of §4 as a table.

---

## 7. Verification (do all five before believing any estimate)

1. **Nesting.** At θ₄₀ = 1, θ₃₉ = mean log κ^cal, the objective must reproduce today's baseline
   **to machine precision**. If it does not, the map is wrong.
2. **SS invariance.** Solve the steady state at 0.1× and 10× κ; every SS value must be identical
   (π_ss = 1). If not, κ has leaked into the SS block and §1's cost argument fails.
3. **Read-back.** κ_i as seen by the solver equals κ_i(θ) at three random θ draws (§6.3).
4. **Sweep non-degeneracy.** Moments must actually *move* across the κ grid. A byte-identical
   sweep means the parameter never arrived — it does **not** mean κ does not matter.
5. **Blanchard–Kahn** at both bounds of θ₃₉ and θ₄₀, and at the implied κ_i extremes.

---

## 8. Recommended sequence

| # | step | new params | new moments | cost |
|---|---|---|---|---|
| 0 | Email Albagli/Luttini re: sectoral Chilean PPI frequencies | — | — | one email, highest payoff |
| 1 | Add block A (`autocorr(π_i)`) to data + model moments | 0 | 11 | ~1 day |
| 2 | Re-run baseline; confirm current κ^cal fits the new block | 0 | — | 1 cluster run |
| 3 | `SMM_KAPPA_SCALE` profile over s | 0 | — | grid of runs |
| 4 | Add blocks B–D | 0 | ~23 | ~2 days |
| 5 | Pre-flight sensitivity scan on θ[39:40] | 2 | — | minutes |
| 6 | Profile λ_κ, then joint CMA-ES (S1) | 2 | — | cluster |
| 7 | S2 / S3 as robustness | 3–11 | — | cluster |

Step 2 is not optional. If the calibrated κ vector already misses `autocorr(π_i)` badly, that
tells you the answer before a single parameter is freed — and it tells you *which* sectors are
driving it.

---

## 9. Open risks to state explicitly in the paper

- κ̂ absorbs σ_i measurement error (§2.3) — mitigated but not eliminated by ratio/conditional
  moments. The shock-inversion refinement is still outstanding independently.
- Quarterly time aggregation caps what is learnable about the flexible sectors (§2.4).
- κ_i, cl, and κ_w all damp the price/quantity split. Report the correlation of their estimates,
  and check whether freeing κ moves ĉl away from the FGI value (19.1) — if it does, they are
  trading off and the identification claim needs more work.

---

## 10. What was built, and how the data changed the design (2026-09-02)

### 10.1 Code

| file | change |
|---|---|
| `utils.jl` | `KAPPA_MODE`/`KAPPA_SHRINK`/`KAPPA_NKPC_SECTORS`, `kappa_from_theta`, `kappa_penalty`; θ 38→**51**; moments 84→**109**; new blocks + weights |
| `smm_estimation.jl` | κ written per evaluation with guards; new data + model moments; 4-check verification suite; `SMM_MAX_EVALS`, `SMM_LAMBDA_KAPPA`, `SMM_LOG_KAPPA_BAR`; hard error if `SMM_KAPPA_SCALE` and `SMM_KAPPA_MODE` are both set |
| `smm_model_moments.jl` | `_KLEIN_STRUCT_IDX` 1:38 → **1:51** |
| `compute_data_moments.jl` | section 9d + `_dmask`/`mdstd`/`mdcor`; new CSV columns `autocorr_PH`, `ratio_dPH_dY`, row `rbar_dPH` |
| `main_SOE_gap.jl` | rebuilds `modkappa` from θ[39:51] so figures and estimation cannot disagree |
| `run_kappa_ladder.sh` | rungs 0–6, local, `--threads=1` |
| `cluster/run_kappa_smm.sh` | sbatch, single-threaded, mode-stamped outputs |
| `compare_kappa_micro.jl` | κ̂ vs κ^cal table (CSV + booktabs) in κ and in implied durations |
| `Data/check_price_moments.py` | independent numpy re-implementation of the new data moments |

θ map, in every mode:
`log κ_i = θ₃₉ + θ₄₀·(log κ_i^cal − μ_cal) + δ_i`, δ = θ[41:51], over the 11 NKPC sectors.
`SMM_KAPPA_MODE` ∈ {off (default), scale {39}, affine {39,40}, free {41:51}}. Every mode nests
the calibration at the seed; `off` reproduces prior runs bit for bit.

### 10.2 Moments actually added (blocks A, B, D)

| # | moment | source |
|---|---|---|
| 85–96 | `autocorr(PH_i)` | Γ_v, Γ1_v |
| 97–108 | `std(ΔPH_i)/std(ΔY_i)` | `var(Δx) = 2(Γ0 − Γ1)` |
| 109 | avg pairwise `corr(ΔPH_i, ΔPH_j)` | `cov(Δx,Δy) = 2(Γ0 − Γ1)` |

All three fall out of the covariance matrices the objective already forms — the returned Γ1 is
symmetrised, which is *exactly* what the differenced-series identity needs. No new frequency
loop, no `.mod` change, no extra steady-state solve.

### 10.3 THE DATA CHANGED THE DESIGN — block A is untargeted by default

Before targeting anything, the two new sectoral blocks were sign-checked against the calibrated
cross-section (2013Q1–2023Q4 common window, pandemic excluded, independent numpy implementation):

| rank corr with log κ_i^cal | all 11 NKPC | 9 market | theory |
|---|---|---|---|
| `autocorr(PH_i)` | **−0.118** | **−0.183** | should be **positive** |
| `std(ΔPH_i)/std(ΔY_i)` | −0.182 | **−0.533** | should be **negative** ✓ |

Leave-one-out, market sectors: the ratio stays in [−0.86, −0.38] and never flips; `autocorr(PH_i)`
stays in [−0.48, +0.12] and is essentially never positive. No single sector drives either.

So the persistence block — the one this plan opened by calling the disqualifying gap — **has the
wrong sign in Chilean data**. The most likely reason is measurement, not economics: a quarterly
sectoral GDP deflator is a *value-added* deflator computed as a residual, noisiest exactly where
output is most imputed (services), and classical measurement error biases an autocorrelation
**down**. The noise is therefore correlated with κ^cal in precisely the wrong direction, and
targeting the block would push κ̂ toward *more flexible services* — backwards.

**Decision:** block 85–96 ships with weight **0** — computed and printed as a diagnostic, the same
treatment already given to `std(L_i)` and the non-market sectors. `SMM_W_AC_PH=1.0` turns it on for
a robustness column; if it is ever turned on, the paper must say it disagrees with the micro
calibration's ordering and why. Identification rests on 97–108 and 109.

Also worth reporting on its own: **`rbar_ΔPH` = −0.017** in the data. Sectoral price changes are
essentially uncorrelated across Chilean sectors. Given the model's known excess output comovement
(rbar 0.41 model vs 0.22 data), expect a large miss here, and expect it to be informative.

### 10.4 What this does to the plan

- §2.1 said the missing price-persistence moment was the disqualifying gap. It was — but the
  gap **cannot be filled with GDP deflators**. That strengthens, rather than weakens, the two
  recommendations that do not depend on deflator quality: block C (pass-through, conditional on an
  exogenous common cost shock) and Tier 0 (micro frequencies from the IPP).
- The identification is now thinner than designed: 12 ratio moments + 1 comovement moment against
  2 free parameters in the headline spec. That is workable for `affine`. It is **not** workable for
  `free` (11 parameters, ~38 usable quarters) except as a shrinkage robustness column.
- The §6 pre-flight sensitivity scan (rung 3) therefore matters more, not less. Run it before
  spending queue time.

### 10.5 First estimate — local affine run, 2026-09-02

260 evaluations, 1344 ms/eval, epsY pinned, `--threads=1`, warm-started.
**Exploratory, not reportable** — CMA-ES self-terminated after 260 evals, which is
too few for 10 free dimensions.

| | seed | estimate | |
|---|---|---|---|
| obj | 9.161 | **8.621** | −5.9% |
| log κ̄ | 4.222 | **3.820** | κ̄ 68.2 → 45.6; mean duration ≈3.5q → ≈2.9q |
| λ_κ | 1.000 | **0.560** | dispersion **compressed**, not widened |
| cl | 7.42 | **19.24** | FGI (2023) estimate is 19.1 (s.e. 12.6) |
| κ_w | 1500 (UB) | **1192** | came **off** the bound |
| ρ_om | 0.990 (UB) | 0.836 | off the bound |

**Read the two κ parameters very differently.**

*The level is real.* log κ̄ falls, |dobj| = 1.23, and the direction was predicted from
the rung-2 fit (model prices too sluggish in every targeted sector).

*The dispersion is not yet an estimate.* λ̂_κ = 0.56 is a large move bought for almost
nothing: |dobj| for θ[40] is 2.2e−2, about 34× less responsive per unit than θ[39], and
**the block it was supposed to be identified from did not improve** — the ratio block went
2.744 → 2.698 while rising from 30.0% to 31.3% of the objective. What actually improved
was corr(Y_i,PH_i) (2.458 → 2.174), autocorr(Y_i) (1.320 → 1.156), rbar_ΔPH (0.180 → 0.081)
and std(L) g/s (0.172 → 0.032). κ is being pulled by everything except its own moment.
That is the signature of a flat direction. **Profile λ_κ (rung 4) before quoting a number.**

Note this also reverses the §10.4 prediction that λ_κ > 1. The prediction came from the
miss factor being larger in flexible sectors; the estimator went the other way. Either
reading needs the profile, not a point estimate.

**Two things the run surfaced that are worth acting on.**

1. **Finance is 9.4% of the total objective and is unreachable.** Data ratio 5.68 against
   0.25–2.88 everywhere else; the residual moved only 5.30 → 5.10 across the whole
   estimation. Financial-services value added is FISIM, so its deflator is a residual of a
   residual. An unreachable moment with a large weighted residual does not identify κ — it
   biases the *level* down, because the only way to shrink it at all is to make every price
   more flexible. `SMM_W_RATIO_DROP=8` untargets it; run both and report both.
2. **κ_w came off its 1500 bound the moment κ was free.** Wage stickiness had been standing
   in for price flexibility. Good news for the labour block, and it means κ and κ_w are
   partly substitutes — report their joint movement, not just κ̂.

Also: the sensitivity scan is **not reliable for parameters sitting on a bound**. The test
point is θ + 0.25·(UB − θ), so κ_w at its UB shows |dobj| = 7.2e−5 and ρ_om at 0.9896 shows
2.7e−4 — neither is evidence of flatness, just of no room to move. Read that column only for
interior parameters.

### 10.6 Run B — Finance untargeted. λ_κ flips, and hits a corner.

3792 evals, 475 ms/eval, flat over the last ~1700 — effectively converged, unlike Run A
(260 evals, CMA-ES self-terminated). Objectives are not comparable across the two: B's is
0.807 lower by construction, that being exactly Finance's contribution.

| | κ^cal | Run A (all moments) | Run B (Finance dropped) |
|---|---|---|---|
| obj | — | 8.621 | 7.376 |
| evals | — | 260 | 3792 |
| log κ̄ | 4.222 | 3.820 | **2.974** |
| λ_κ | 1.000 | 0.560 | **1.713** |
| sd(log κ) | 1.281 | 0.717 | **2.193** |
| mean duration | 3.98q | 2.99q | 3.17q |
| cl | — | 19.24 | 6.11 |
| κ_w | — | 1192 | 1499 (back at UB) |

**λ_κ flipped sides, 0.56 → 1.71.** Compression became amplification. Two things changed at
once, though — the Finance moment *and* 15× more search — so this is not a clean decomposition.
What can be said: in B the ratio block finally *improved* (1.873 → 1.470, −22%), which it never
did in A. With Finance out, κ moves the moment it was built for. That is the first evidence the
identification works as designed.

**But λ̂_κ = 1.713 is a corner solution against a numerical guard, not an interior estimate.**
At that dispersion the implied κ are:

| | Agri | Manuf | Trade | Utilit | Transp | BusServ | Fin | Constr | RealEst | PersServ | PubAdmin |
|---|---|---|---|---|---|---|---|---|---|---|---|
| κ^cal | 8.0 | 11.0 | 20.4 | 33.9 | 70.3 | 85.1 | 95.8 | 136.5 | 178.8 | 337.6 | 508.6 |
| κ̂ | **0.50** | 0.86 | 2.48 | 5.92 | 20.7 | 28.7 | 35.1 | 64.4 | 102.2 | 303.5 | 612.4 |
| duration (q) | 1.05 | 1.09 | 1.23 | 1.46 | 2.11 | 2.37 | 2.56 | 3.26 | 3.98 | 6.55 | 9.22 |

Agriculture lands on **exactly** the `kappa < 0.5 → reject` guard in `smm_model_moments`. The
8% failure rate (312 of 3792) is the search repeatedly bumping into it. Raising λ_κ further
would push κ_1 below the floor and be rejected, so **the floor caps λ_κ from above**: the data
want at least this much dispersion and possibly more. λ̂ = 1.713 is a lower bound, and must be
reported as one.

The floor is not economically silly — κ = 0.5 is a 1.05-quarter duration, i.e. full flexibility
at quarterly frequency, which is the finest a quarterly model can resolve. But that is the
substantive statement: **the data want Agriculture and Manufacturing fully flexible.**

**This is the argument for the "free" mode.** The affine map has one dial for eleven sectors, so
buying flexibility for Agriculture forces Public Admin to κ = 612 — *above* its calibrated 509,
which nothing in the data asked for. Per-sector deviations (θ[41:51], with the ridge) can make
the flexible end flexible without dragging the sticky end further out.

Two other movements to keep an eye on, both of which say the runs are landing in different
basins rather than converging on one answer: **cl** went 19.24 (A) → 6.11 (B), and **κ_w** went
back to its 1500 upper bound in B after coming off it in A. cl and κ_w and κ are all damping the
same price/quantity margin.

### 10.7 Next actions, in order

1. `bash run_kappa_ladder.sh 0 1 2` — regenerate moments, verify the wiring, score the baseline.
   **Rung 1 must show all four checks OK before anything below is meaningful.**
2. `bash run_kappa_ladder.sh 3` — the sensitivity scan. If |dobj| for θ[39:40] is at the log κ_V
   level (1e-4), stop and report flatness.
3. `bash run_kappa_ladder.sh 4` — profile λ_κ. This is the paper figure.
4. `sbatch cluster/run_kappa_smm.sh affine` — the headline estimate.
5. `sbatch cluster/run_kappa_smm.sh free 2.0` — the robustness column, then
   `julia --project=. compare_kappa_micro.jl`.
6. In parallel, and independently of all of the above: **email Albagli/Luttini** about sectoral
   Chilean PPI frequencies (§5). It is still the highest-payoff item in this document.
