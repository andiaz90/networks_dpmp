# SMM Estimation of the NK-IOSOE Chile Model — Method Note

*Agustín Díaz · Central Bank of Chile · 2026-07-10*
*Describes the current Julia implementation (`run_smm_estimation.jl`, `smm_estimation.jl`, `smm_model_moments.jl`, shared definitions in `utils.jl`). Supersedes the 15-moment MATLAB description in the paper draft — the paper's estimation section is pending an update.*

---

## 1. What is being estimated, in one paragraph

We estimate 36 structural and shock parameters of the 12-sector NK-IOSOE model by minimizing the weighted distance between 60 second moments of HP-filtered Chilean quarterly data and their model counterparts. Despite the "SMM" label, no simulation is involved: for every candidate parameter vector θ the model is re-solved to first order (Klein 2000), and the model moments are computed *analytically* from the state-space representation — unconditional HP-filtered covariances obtained from the Lyapunov equation with a spectral HP filter (λ = 1600, 256 frequencies). This makes each evaluation exact (no simulation noise) and fast (~0.2 s), which is what allows a global optimizer (CMA-ES) over a 36-dimensional parameter space.

Formally: θ̂ = argmin ψ(θ)′ W ψ(θ), where ψ(θ) = m_data − m_model(θ), and W is a fixed diagonal weighting matrix (Section 5).

## 2. The 36 estimated parameters

Layout of θ (single source of truth: `PARAM_LABELS`, `LB`, `UB` in `utils.jl`):

| θ index | Parameter | Model object | Bounds | Notes |
|---|---|---|---|---|
| 1 | `ilabcosts` | Sectoral labor reallocation cost scale (multiplies cl_i in the labor agency problem) | [0.001, 50] | The "c" of the quadratic switching cost |
| 2 | `epsY` | ε_Y — elasticity of substitution across production inputs (materials / imported input / labor) | [0.30, 1.50] | < 1 ⇒ gross complements |
| 3 | `epsM` | ε_m — elasticity across sectoral materials within the M bundle | [0.05, 0.50] | Governs IO network propagation |
| 4 | `log(kappaV)` | κ_v — Rotemberg adjustment cost on import prices (estimated in logs) | [log 10³, log 10⁸] | Large κ_v ⇒ sticky import prices, incomplete ERPT |
| 5 | `rho_om` | Common AR(1) persistence of the 12 sectoral demand shocks ω_i | [−0.95, 0.95] | |
| 6 | `rho_A` | Common AR(1) persistence of the 12 sectoral TFP shocks A_i | [0.10, 0.95] | |
| 7–18 | `isigma_tfp_1..12` | Sectoral TFP innovation std devs σ_{A,i} | [10⁻⁴, 0.10] | One per sector |
| 19–30 | `sigma_om_1..12` | Sectoral demand innovation std devs σ_{ω,i} | [10⁻⁵, 0.20] | One per sector |
| 31 | `rho_pvstar` | AR(1) of world imported-input price PV* | [0.50, 0.99] | |
| 32 | `sigma_pvstar` | Innovation std of PV* | [0.005, 0.20] | |
| 33 | `rho_xi` | AR(1) of intertemporal preference shock ξ (Euler equation) | [0.00, 0.95] | Aggregate demand shock |
| 34 | `sigma_xi` | Innovation std of ξ | [0, 0.05] | |
| 35 | `etastar` | η* — foreign demand elasticity for exports | [0.50, 6.00] | |
| 36 | `kappaw` | Rotemberg wage adjustment cost (EHL-style wage PC, added 2026-07-08) | [0, 400] | 0 nests flexible wages; 115 ≈ 4-quarter Calvo at ε_w = 10 |

Everything else is calibrated outside the estimation: preferences (γ = 2, ψ, β = 0.986), ε = 10, sectoral Rotemberg price costs κ_i (Pastén et al.), the IO matrix Γ and factor shares α_mi, α_vi (2021 IO tables), consumption weights (Ω_G = 0.57, ω_i), home bias ϱ_i, Taylor rule, ε_w = 10, and the steady-state external targets (TB/GDP = 2.3%). Monetary (eps_i), PV*, ξ, TFP, and demand shocks are active in the moment computation (27 of 29 exogenous shocks; `epschi` and `eps_postar` are excluded).

## 3. The data targets (60 moments)

All data moments are computed by `compute_data_moments.jl` from BCCh/INE/BIS series, HP-filtered with λ = 1600, and stored in `Data/sectoral_moments.csv` and `Data/aggregate_moments.csv`. Nominal sample 2006Q1–2023Q4; the effective common window is shorter for some sectors (Finance and Business Services output start 2013Q1; stds are computed on valid windows only). ⚠ The paper still says 2006Q1 everywhere — pending fix.

**Block A — sectoral volatilities (36 moments).** std of HP-filtered log real output std(Y_i), producer price deflator std(PH_i), and employment std(L_i) for each of the 12 sectors. Sources: `pib_sectorial_bc` (output), `deflactor_pib` (prices), `count_workers_by_sector` (employment, monthly → quarterly).

**Block B — aggregates (7 moments).** std(GDP) = 0.041, std(π) = 0.022 (QoQ deflator), corr(GDP, π) = −0.02, std(TB/GDP) = 0.035, std(Q) = 0.043 (BIS REER), autocorr(Q) = 0.72, corr(GDP, Q) = 0.14.

**Block C — cross-sector rank correlations (3 moments).** Spearman rank correlation between the model's and the data's cross-sectoral ordering of output / price / labor volatilities, each with target 1. These discipline *which sectors* are volatile, not just average volatility — the model must reproduce that Agriculture is the most volatile sector, Public Administration the least, etc.

**Block D — sectoral price–quantity comovement (12 moments).** corr(Y_i, PH_i) per sector. This is the supply-vs-demand identifier: sectoral TFP shocks push output and prices in opposite directions (negative corr), sectoral demand shocks push them together (positive corr). Data values range from −0.97 (Personal Services) to +0.21 (Construction).

**Block E — labor comovement (2 moments, added 2026-07-08).** corr(N, GDP) = +0.71 and corr(N, GDP/N) = +0.08. Motivated by the Galí facts: employment is strongly procyclical while employment–productivity correlation is ≈ 0. Discipline the aggregate demand/supply mix and the wage stickiness κ_w.

## 4. Identification: which moments pin down which parameters

Heuristic mapping (all parameters affect all moments in general equilibrium; this is where the *first-order* information comes from):

| Parameter(s) | Primary identifying moments | Mechanism |
|---|---|---|
| σ_{A,i} (12) | std(Y_i), output rank corr | Sectoral TFP variance scales own-sector output volatility directly |
| σ_{ω,i} (12) | corr(Y_i,PH_i), std(PH_i), price rank corr | Demand shocks generate *positive* price–quantity comovement; their size relative to σ_{A,i} sets the sign and magnitude of corr(Y_i,PH_i) |
| ρ_A, ρ_ω | Relative volatility of HP-filtered stds across blocks | Persistence shifts variance toward low frequencies that the HP filter removes; affects all stds jointly |
| ilabcosts | std(L_i) relative to std(Y_i), labor rank corr | Higher reallocation cost dampens sectoral employment responses |
| κ_w | std(L_i) levels, corr(N,GDP), corr(N,GDP/N), labor rank corr | Sticky wages damp the MRS channel: labor stops overreacting to TFP; raises corr(N,GDP) toward the data |
| ε_Y, ε_m | Cross-sector comovement, price rank corr | Substitutability governs how shocks propagate through the IO network to other sectors' MC and prices |
| κ_v | std(PH_i) of import-intensive sectors, corr(GDP,π) | Import price stickiness controls ERPT into domestic marginal costs |
| ρ_pv*, σ_pv* | std(Q), autocorr(Q), corr(GDP,Q), std(TB/GDP) | The PV* shock is the main driver of RER dynamics and the trade balance |
| η* | std(TB/GDP) | Export demand elasticity scales the TB response to relative prices |
| ρ_ξ, σ_ξ | std(GDP), corr(GDP,π), corr(N,GDP) | Aggregate demand shock: raises output-inflation comovement and employment procyclicality |

Sanity checks worth running on any estimate: parameters at bounds (the report flags `<< at LB` / `>> at UB`) indicate either weak identification or a binding economic constraint; a σ_{ω,i} at its lower bound says the data want sector i to be purely supply-driven.

## 5. Weighting matrix

W is diagonal and fixed (not the efficient two-step GMM matrix):

| Moments | Weight | Rationale |
|---|---|---|
| All 40 std devs (blocks A + std moments of B) | 1/d² (d = data value; floor at d = 10⁻⁴) | Converts errors to *percent* deviations — scale-free across moments of very different magnitude |
| corr(GDP,π) | 0.2 × 1/d² | Structurally hard for a supply-shock-driven model; down-weighted |
| autocorr(Q) | 2.0 | Fixed |
| Rank correlations (3) | 2.0 each | Fixed |
| corr(Y_i,PH_i) (12) | 3.0 each | Fixed since 2026-07-09. The old proportional 1/d² weights exploded for near-zero data correlations (Business Services: d = 0.007 ⇒ w ≈ 30,800) and made this block 99.5% of the objective. Correlations live on [−1,1], so absolute errors are already comparable |
| corr(N,GDP), corr(N,GDP/N) | 2.0 each | Fixed |

⚠ Objective values are **not comparable across the 2026-07-09 reweighting** (or across any change to the moment set). The current objective at the pre-reweighting checkpoint θ is the baseline to beat.

Standard errors: `smm_inference.jl` computes them from a numerical Jacobian of the moment function at θ̂ (delta method with the fixed W). With diagonal W these are valid but not efficient; the J-statistic is not chi-squared distributed under the fixed W and should not be reported as an overidentification test.

## 6. Computational pipeline

```
compute_data_moments.jl      → Data/{sectoral,aggregate}_moments.csv
main_SOE_gap.jl              → steady state, params_jl.mod, Dynare solve (subprocess),
                               mod/nk_iosoe_context.jls  +  full fit report
run_smm_estimation.jl        → CMA-ES estimation
       ├── pre-flight: re-solve at θ₀, print parameter + fit tables, abort if solve fails
       ├── CMA-ES on [0,1]-rescaled θ, multithreaded (per-thread Dynare contexts)
       │     failed solves → penalty 1e8; checkpoint written on EVERY improvement
       └── outputs: Data/smm_estimates.csv (final θ̂), Data/smm_results.csv (60-moment fit),
                    Data/smm_checkpoint.csv (live best; warm-start source)
```

Practical notes:

- **Warm start:** any run automatically restarts from `smm_checkpoint.csv` if it has the current 36-parameter layout.
- **Fit report without estimating:** `SMM_REPORT_ONLY=1 julia --threads=auto run_smm_estimation.jl` prints the parameter table (with bound flags) and the 60-moment fit at the current checkpoint θ, touching nothing.
- **`main_SOE_gap.jl`** now auto-loads the newer of checkpoint/estimates and prints the identical 60-moment table with the true SMM weights, so its TOTAL is the SMM objective.
- **Wall clock:** set `SMM_MAX_HOURS` slightly under the SLURM limit; the run exits gracefully and saves final state.
- **Cluster:** 48 threads, ~0.18 s/eval locally (4 threads); production runs target 300,000 evaluations max.

## 7. Known caveats / open items (as of 2026-07-10)

1. The paper's estimation section describes the old 15-moment MATLAB implementation — must be rewritten to match this note.
2. Effective data sample is 2009Q1–2023Q4 for the common window (paper says 2006Q1) — fix in draft.
3. Diagonal W: efficient (optimal) weighting and a proper J-test are not implemented; SEs are delta-method under fixed W.
4. No estimation results are final yet: the only complete runs predate the sticky-wage extension and the 2026-07-09 reweighting. Current checkpoint θ is an interrupted, pre-reweighting point with `ilabcosts` near its upper bound and ρ_ω < 0 — treat as a starting value, not as estimates.
5. Excluded shocks: `epschi` (χ, labor disutility) and `eps_postar` (world oil price) are off during estimation; the oil block is used for the policy exercises only.
