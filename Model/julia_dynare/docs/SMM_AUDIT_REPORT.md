# SMM Methodology Audit — NK-IOSOE Chile Model

**Date:** April 15, 2026

---

## Overall Assessment

The SMM implementation is **mostly correct** but has **three issues** that affect the quality of the estimates. Two are methodologically important (Issues 1 and 2) and one is a potential bug (Issue 3).

---

## Issue 1: HP Filter Mismatch (IMPORTANT)

**Data moments:** `y_d[i] = std(hp_cycle(log(Y_i), λ=1600))`
**Model moments:** `pstd("Y_i") = sqrt(Var(Y_i)) / Y_i_ss`

The data moments are computed on **HP-filtered log deviations**. The model moments are **unconditional analytical second moments** from the discrete Lyapunov equation, which are the population analogs of `std(log(Y_i))` — there is *no* HP filter applied to the model.

**Why this matters:** HP filtering removes low-frequency variation. In a DSGE model with persistent shocks (ρ_A = 0.5–0.9), a substantial fraction of output variance is at low frequencies. The HP filter removes this variance from the data but the model's Lyapunov moments include it. The result: the model systematically **overstates** the volatility that should be matched to the data, leading the optimizer to **underestimate** shock standard deviations (isigma_tfp_i, sigma_om, sigma_pvstar).

The bias is proportional to shock persistence. For ρ_A = 0.8, the ratio of HP-filtered to unfiltered variance is approximately:

```
Var_HP / Var_unfilt ≈ 1 - (1-ρ²)/(1-ρ²+λ(1-ρ)⁴) ≈ 0.60
```

So the model's unfiltered std is ~√(1/0.60) ≈ 1.29× larger than the HP-filtered data std. The optimizer compensates by shrinking shock sizes by ~23%.

**Fix:** Apply the HP filter to the model's theoretical spectrum. For a first-order solution with state-space `y_t = T·y_{t-1} + R·ε_t`, the HP-filtered autocovariance can be computed analytically by filtering the spectral density:

```
S_HP(ω) = |1 - H(ω)|² · S_y(ω)
```

where `H(ω) = 1/(1 + λ(2-2cos(ω))²)` is the HP filter's transfer function. The filtered variance is then `∫ S_HP(ω) dω / 2π`, which can be evaluated numerically on a grid of 512 frequencies in ~1ms.

Alternatively (simpler): simulate a long time series (T=10,000) from the linearized model, apply the HP filter, compute sample moments. This is slower but trivial to implement.

**References:** Burnside (1998, JMCB); Cogley and Nason (1995, AER) — standard treatment of HP-filtered model moments in DSGE estimation.

---

## Issue 2: Monetary Policy Shock Excluded from Σe (MODERATE)

The model has a monetary policy shock `eps_i` (varexo index 2) with AR(1) process:

```
vi = rhoi*vi(-1) + sigma_i*eps_i
```

But in `smm_model_moments.jl` (line 722), the Σe matrix is constructed with:

```julia
Σe[1, 1]  = 1.0    # eps_om
Σe[4, 4]  = 1.0    # eps_pvstar
Σe[5:16]  = 1.0    # epsA_1:12
Σe[17,17] = 1.0    # eps_xi
```

**`Σe[2, 2] = 0.0` — the monetary policy shock `eps_i` is never activated.**

Similarly, `Σe[3, 3] = 0.0` — the labor supply shock `epschi` is inactive.

The monetary policy shock `eps_i` is not estimated (its parameters `rhoi`, `sigma_i` are not in θ), but `sigma_i` is set to a calibrated value in `params_jl.mod`. If `sigma_i > 0` in the calibration, the model *should* include this shock in Σe so it contributes to aggregate volatility. Otherwise, all aggregate fluctuations must come from TFP, preference, and import price shocks — forcing the TFP shocks to be too large to compensate for the missing monetary policy channel.

**Impact:** If `sigma_i = 0.001` (the calibrated value from the paper's Table 3), the monetary policy shock contributes modestly to GDP and inflation volatility. Excluding it biases isigma_tfp upward and distorts corr(GDP, π) because monetary shocks generate positive GDP-inflation comovement (demand shock via the Taylor rule).

**Fix:** Add `Σe[2, 2] = 1.0` to activate the monetary policy shock. If you don't want to estimate (rhoi, sigma_i), keep them fixed at their calibrated values — but they should still enter the model's moment computation.

---

## Issue 3: Lag-1 Covariance Formula May Be Incorrect (CHECK)

The autocorrelation of Q (moment 42) is computed as:

```julia
Γ1 = T * T[sr,:] * (P * T' + R[sr,:] * Σe * R')    # Julia
Γ_1 = T * A * (P_state * T' + B * Qe * R')           # MATLAB
```

The correct formula for the lag-1 cross-covariance `E[y_t · y_{t-1}']` in the state-space `y_t = T·s_t, s_t = A·s_{t-1} + B·ε_t` is:

```
E[y_t · y_{t-1}'] = T · E[s_t · y_{t-1}']
                   = T · E[(A·s_{t-1} + B·ε_t)(s_{t-1}'·T' + ε_{t-1}'·R')]
                   = T · (A · P · T' + 0)           [since E[ε_t · ε_{t-1}'] = 0]
                   = T · A · P · T'
```

But the code computes `T * A * (P * T' + B * Qe * R')`. The extra term `B * Qe * R'` equals `E[ε_t · ε_t'] * R'` which would contribute only if `ε_t` and `y_t` are contemporaneously correlated at lag 0, not lag 1. This term should be zero because `E[B·ε_t · (R·ε_{t-1})'] = B · E[ε_t · ε_{t-1}'] · R' = 0`.

**However**, the issue depends on the timing convention. If `y_t` includes contemporaneous `R·ε_t` and we're computing `E[y_t · s_{t-1}']`:

```
E[y_t · s_{t-1}'] = E[(T·s_t + R·ε_t) · s_{t-1}']
                   = T·A·P + R·E[ε_t · s_{t-1}'] = T·A·P
```

This confirms the extra term is spurious. **But** it only matters for `autocorr(Q)`, which is a single moment (moment 42). The contamination would make autocorr(Q) slightly too high, causing rho_pvstar to be slightly underestimated. The quantitative effect is likely small (since `B*Qe*R'` is typically much smaller than `A*P*T'`), but it should be corrected.

**Fix:**
```julia
Γ1 = T * T[sr,:] * P * T'
```

---

## Things That Are Correct

### ✓ Parameter → Moment Mapping
The 23 estimated parameters map logically to the 46 moments:
- `isigma_tfp_1:12` → `std(Y_1:12)` (12 output volatilities, well-identified)
- `epsY, epsM` → `std(PH_1:12)` (price volatilities, over-identified through IO)
- `ilabcosts` → `std(L_1:12)` (labor volatilities, single aggregate cost)
- `rho_om, sigma_om` → sectoral reallocation + goods expenditure share
- `rho_pvstar, sigma_pvstar` → `std(Q)`, `autocorr(Q)`, `corr(GDP,Q)`
- `rho_xi, sigma_xi` → `corr(GDP,pi)` (demand shock for comovement)
- `log(kappaV)` → import price pass-through → exchange rate moments
- `rho_A` → persistence of all sectoral TFP → autocorrelation structure

This is standard and appropriate for a 12-sector SOE model.

### ✓ Lyapunov Approach (not simulation)
Using the analytical Lyapunov equation `P = A·P·A' + Q` to compute unconditional moments is correct for a first-order perturbation. This avoids simulation error and makes the objective smooth (no Monte Carlo noise). This is standard in the SMM literature.

### ✓ SS Recomputation When epsY/epsM Change
The production function CES elasticities (epsY, epsM) enter the steady state nonlinearly — changing them requires re-solving the SS system. The code correctly detects when these change and re-solves via `nlsolve`/`fsolve`. The SS quantities (prices, quantities, aggregates) are then written back to the Dynare parameters before calling `resol()`. This is necessary and correct.

### ✓ Rank Correlations as Moments
Including `Spearman_rho(std_Y_model, std_Y_data)` as a target moment is a good econometric choice — it ensures the model reproduces the *cross-sectional ordering* of sectoral volatilities, not just the average level. This is important for the paper's narrative about heterogeneous sectoral responses.

### ✓ Percentage Std Dev Normalization
`pstd("Y_i") = sqrt(Var(Y_i)) / |Y_i_ss|` correctly computes the percentage standard deviation of a level variable from the linearized variance. For the first-order approximation `Y_i,t ≈ Y_i_ss·(1 + ŷ_i,t)`, this gives `std(ŷ_i)` ≈ `std(log(Y_i))`, which is what the data moments measure.

### ✓ CMA-ES Optimizer Choice
CMA-ES is appropriate for this problem: 23 parameters, non-differentiable objective (due to Schur decomposition and potential BK violations), and no analytical gradient. The population-based search with covariance adaptation is the standard in structural macro estimation.

### ✓ Weighting Matrix
The diagonal weighting with `W_ii = 1/d_i²` (proportional fit) is standard. The downweighting of `corr(GDP,π)` is justified by the structural supply-shock limitation (discussed in the paper).

---

## Priority Ranking

| # | Issue | Impact | Difficulty |
|---|-------|--------|------------|
| 1 | HP filter mismatch | HIGH — biases all shock σ's downward | Medium (spectral filter or simulation) |
| 2 | Monetary shock excluded | MODERATE — biases TFP σ's upward | Trivial (one line) |
| 3 | Lag-1 covariance formula | LOW — affects only autocorr(Q) | Trivial (remove extra term) |
