# SMM Estimation Diagnostic Report: NK-IOSOE Chile Model

**Date:** April 15, 2026  
**Author:** Analysis for Agustín Díaz

---

## 1. Why the Model Struggles to Match the Data

### 1.1 Structural Identification Problem: 23 Parameters → 46 Moments

The estimation targets 46 moments with 23 parameters. While formal overidentification is desirable, the *effective* identification is much weaker than it appears:

- **Moments 1–12 (std(Y_i))** are primarily identified by **isigma_tfp_1–12** (12 free σ's map to 12 targets). This is well-identified.
- **Moments 13–24 (std(PH_i))** are *over-identified*: sectoral price volatilities depend on the IO structure, κ_i (fixed), epsY, epsM, and the TFP shocks. With κ_i fixed from Pastén et al., the model can only adjust epsY/epsM globally — it cannot match sector-specific price dispersion independently.
- **Moments 25–36 (std(L_i))** depend on ilabcosts (one scalar) and the IO structure. A single aggregate labor adjustment cost cannot match 12 heterogeneous labor volatilities.
- **Moment 39 (corr(GDP,π))** is structurally hard: a supply-shock-only model generates negative corr(GDP,π), but the data may show near-zero or positive correlation. The preference shock (ξ) was added to address this, but it enters through the Euler equation only — its ability to generate positive comovement is limited by the Taylor rule's strong inflation response (φ_π = 2.5).

**Diagnosis:** The model is *structurally unable* to match all 46 moments simultaneously. The main tensions are: (a) sector-specific price volatilities vs. two common elasticities (epsY, epsM); (b) sector-specific labor volatilities vs. one aggregate adjustment cost (ilabcosts); (c) the corr(GDP,π) sign problem.

### 1.2 The Steady-State Recomputation Bottleneck

Every time epsY or epsM change (which is *every* CMA-ES iteration, since they're free parameters), `recompute_ss()` is called. This function:

1. Calls `fsolve` on the 15-variable outer system (12 pH + w + Q + C) — **~0.3s per call**
2. Inside that, calls `fsolve` on the 48-variable production system (M, L, V, Y × 12) — **~0.1s per call**
3. Then updates all Dynare SS parameters and calls `resol()` — **~0.05s per call**

With CMA-ES population size λ ≈ 4 + 3*ln(23) ≈ 13, each generation evaluates ~13 candidates. At ~0.5s per evaluation, that's ~6.5s per generation. With `MaxFunEvals = 115,000`, total wall time is approximately **115,000 × 0.5s ≈ 16 hours** — and that's if every evaluation succeeds.

### 1.3 The Weighting Matrix Distortion

The current diagonal weighting `W_ii = 1/max(|d_i|, 0.01)^2` with manual fixes creates an *irregular* objective surface:

- Small moments (std(PH_10) ≈ 0.017) still get weight ~870 even after the 0.25× cap, while large moments (std(Y_2) for Mining ≈ 0.10) get weight ~100. The 9:1 ratio means price moments of business services contribute 9× more to the loss than mining output volatility.
- The 0.02× scaling on corr(GDP,π) effectively removes this moment from the objective — the optimizer ignores it entirely.
- Rank correlation moments (44–46) target 1.0, so their natural weight is 1/1² = 1, then halved to 0.5. This is 1700× less than the weight on std(PH_10). The rank correlations are essentially ignored.

### 1.4 The epsY–epsM Flat Ridge

The production function CES with epsY (gross elasticity) and epsM (materials elasticity) creates a near-flat ridge in the objective. Both parameters control how intermediate inputs substitute with labor and imports. In a 12-sector IO network, their effects are highly collinear: increasing epsY while decreasing epsM (or vice versa) produces similar aggregate moments. CMA-ES explores this ridge inefficiently because the covariance matrix adaptation takes many generations to align with a narrow valley.

### 1.5 Missing Demand Shocks for GDP-Inflation Comovement

The ξ (preference) shock enters the Euler equation: `ξ_t C_t^{-γ} = β E_t[ξ_{t+1} C_{t+1}^{-γ} r_t/π_{t+1}]`. This creates a demand impulse (ξ↑ → C↑ → π↑), but the Taylor rule's strong response (φ_π = 2.5) dampens the inflation effect quickly. The result: ξ generates mild positive corr(GDP,π), but not enough to offset the dominant negative correlation from TFP shocks. The model needs either a government spending shock or a stronger demand channel.

---

## 2. Why Convergence Is Slow

### 2.1 Redundant SS Recomputation

`need_ss_update` in `smm_model_moments.m` (line 114) compares epsY/epsM against `M_in` (the *original* compiled model), not against the previous iteration. Since CMA-ES draws epsY/epsM from a continuous distribution, **every single evaluation triggers a full SS recomputation** even when epsY/epsM barely change. This is the single largest performance bottleneck.

### 2.2 No Parallelism in CMA-ES

The standard `cmaes.m` (Hansen 2012) evaluates candidates sequentially. Each generation's λ ≈ 13 candidates are independent — they could be evaluated in parallel using MATLAB's `parfor`. This alone would give a 4–8× speedup on a modern workstation (4–8 cores).

### 2.3 Oversized Search Volume

`insigma = (ub - lb)/6` covers the full feasible range. For parameters like `log(kappaV)` with range [log(1e3), log(1e16)] = [6.9, 36.8], the initial σ = 5.0. This means CMA-ES spends many early generations exploring extreme kappaV values (e.g., 10^15) that always return `resol` failures (penalty 1e8). The optimizer needs hundreds of generations just to shrink the covariance to a sensible region.

### 2.4 High Penalty Rate

With 23 parameters and wide bounds, a large fraction of CMA-ES candidates violate Blanchard-Kahn conditions or produce numerical instabilities. Based on the code's `smm_fail_count` mechanism, failure rates of 30–50% are typical. Each failed evaluation still costs the `resol()` call time (~0.05s) but provides no gradient information to CMA-ES.

---

## 3. Recommended Improvements

### 3.1 Cache the SS to Avoid Redundant Recomputation (HIGH IMPACT)

Store the last (epsY, epsM) and corresponding SS solution. Only recompute when the change exceeds a tolerance (e.g., |Δ| > 1e-4). Use the previous SS as the initial guess for fsolve — this warm-start typically converges in 2–3 iterations instead of 15–20.

### 3.2 Parallelize CMA-ES Population Evaluation (HIGH IMPACT)

Replace the sequential objective evaluation in cmaes.m with parfor-based parallel evaluation. Expected speedup: 4–8× on a standard workstation, more on the BCCh cluster.

### 3.3 Tighten Parameter Bounds (MEDIUM IMPACT)

Based on economic priors and the literature:
- `epsY ∈ [0.3, 1.5]` (not [0.1, 3.0]) — production CES elasticities above 1.5 are implausible
- `epsM ∈ [0.05, 0.5]` (not [0.01, 1.5]) — materials are complements
- `log(kappaV) ∈ [log(1e3), log(1e8)]` (not log(1e16)) — extremely high import adj. costs are degenerate
- `isigma_tfp_i ∈ [1e-4, 0.10]` (not 0.50) — TFP shocks > 10% per quarter are unrealistic

### 3.4 Use Two-Stage Weighting (MEDIUM IMPACT)

Stage 1: Identity weighting (all moments contribute equally in absolute terms). Run CMA-ES for ~20,000 evaluations to get θ̂₁.
Stage 2: Estimate optimal weighting from model moments at θ̂₁ (bootstrap or analytical), then re-estimate.

### 3.5 Profile-Out epsY/epsM via Concentrated Objective (HIGH IMPACT)

For each candidate (epsY, epsM), the SS can be solved once and the remaining 21 parameters optimized conditional on that SS. This separates the "expensive" SS-dependent parameters from the "cheap" shock-process parameters. Implement as a two-level optimizer: outer loop over (epsY, epsM) with grid or Nelder-Mead; inner loop over remaining 21 parameters with CMA-ES (no SS recomputation needed).

### 3.6 Add a Government Spending Shock (STRUCTURAL)

To match corr(GDP,π) > 0, the model needs a proper demand shock. A government spending shock `G_t` that enters the resource constraint (`GDP_t = C_t + G_t + TB_t`) with AR(1) process would directly generate positive GDP-inflation comovement. Two new parameters: (ρ_G, σ_G).

---

## 4. Implementation Priority

| Priority | Change | Expected Impact |
|----------|--------|-----------------|
| 1 | SS caching + warm-start | 3–5× faster per evaluation |
| 2 | Parallel CMA-ES | 4–8× faster total |
| 3 | Tighter bounds | 30–50% fewer failed evaluations |
| 4 | Two-level optimizer (profile-out epsY/epsM) | Removes SS bottleneck entirely |
| 5 | Better weighting matrix | Improved moment fit |
| 6 | Gov spending shock | Fixes corr(GDP,π) structural mismatch |
