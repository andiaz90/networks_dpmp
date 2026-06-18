# NK-IOSOE — GMM/SMM Implementation & Model-Solution Audit

**Date:** 2026-06-18 · **Scope:** `Model/julia_dynare/` (estimation + model solve) · **Branch:** `agustin_development`

This audit answers two questions:

1. How does your estimation differ from the method-of-moments machinery Dynare ships?
2. Is the model being solved correctly?

The short answers: (1) your routine is an *analytical-moment GMM* (despite the "SMM" name), more flexible than Dynare's `method_of_moments` on the moment side but missing its inference machinery, and it re-implements the model solution and moment computation by hand — which is where the problems are; (2) **the model itself is solved correctly** (steady state exact, Blanchard–Kahn satisfied), but the **estimator is wired to a different model version than the one it actually loads**, so the moments it minimises are computed from the wrong structural shocks and several estimated parameters do nothing.

---

## 0. Executive summary of findings

| # | Severity | Finding | Always active? |
|---|----------|---------|----------------|
| **C1** | 🔴 Critical | Shock layout mismatch: code assumes the 28-shock `_smm.mod`; the loaded context is the 18-shock `NK_SOE_lev_gap2`. `Σe` activates the wrong shocks. | Yes (every evaluation) |
| **C2** | 🔴 Critical | The lev model's demand shock is a single `eps_om` scaled by scalar `sigma_om`. The estimator writes `sigma_om_1..12` (which don't exist → silent no-op) and never sets `sigma_om` (also not written by `params_jl.mod`). The whole demand-shock channel is effectively dead/uncontrolled. | Yes |
| **C3** | 🟠 High | Hand-rolled Klein solver hardcodes stale dimensions `491/78/56/17`; the real model is `570/80/58/18`. | Only if Klein fallback runs |
| **C4** | 🟠 High | `_load_jacobian_structure!` zips `g1_v` (column-major, from the compiled file) with `dynamic.json` entries (row-major) → scrambled Jacobian. | Only if Klein fallback runs |
| **C5** | 🟡 Medium | Which context is loaded is decided by file existence (`nk_iosoe_context.jls` vs `nk_iosoe_smm_context.jls`); silently determines whether C1/C2 bite. | Yes |
| **C6** | 🟢 Verified-OK | The model solve is correct: SS residual 3.3e-13; all 80 state-transition eigenvalues < 1 (max 0.9951). | — |

**Bottom line:** the linear solution Dynare produces is correct, but the SMM layer around it is matching moments from a scrambled shock set and estimating ~12 parameters that have no effect. Estimates produced in this configuration should not be trusted until C1/C2 are fixed.

---

## 1. Your estimator vs. Dynare's `method_of_moments`

### 1.1 What you built

Your routine, despite the "SMM" label, does **not simulate**. For each θ it:

1. solves the model to first order (Dynare.jl's native solver, with a hand-rolled Klein fallback),
2. forms the **analytical** unconditional second moments from the discrete Lyapunov equation `P = A P A' + BΣB'`,
3. **HP-filters the model's spectral density analytically** (`hp_filtered_cov_fast`, the `|1−H(ω)|²·S_y(ω)` weighting), so model moments are comparable to your HP-filtered data moments,
4. minimises a diagonally-weighted distance with **CMA-ES**.

So methodologically it is **GMM with model-implied analytical moments**, plus an HP-filter correction.

### 1.2 What Dynare offers

Dynare (MATLAB) has a mature `method_of_moments` command:

- **GMM**: analytical moments from the (optionally pruned) perturbation solution (Andreasen, Fernández-Villaverde & Rubio-Ramírez 2018), orders 1–3; matches means, (co)variances and autocovariances declared in a `matched_moments` block.
- **SMM**: simulates long samples and matches simulated moments.
- **Inference machinery you currently don't have**: optimal weighting matrix (Newey–West/Bartlett HAC), one-step / two-step / iterated GMM, analytical or numerical **standard errors**, and the **J-test of overidentifying restrictions**. Multiple optimizers (including CMA-ES, `mode_compute=9`), parameter bounds.

**Important:** Dynare.jl (the Julia port you use) does **not** ship a complete `method_of_moments`. That is the legitimate reason a hand-rolled routine exists at all.

### 1.3 Head-to-head

| Dimension | Your code | Dynare `method_of_moments` |
|---|---|---|
| Moment source | Analytical (Lyapunov), order 1 | GMM analytical (orders 1–3) or SMM simulated |
| Model filtering | **HP-filters the model spectrum** (matches HP-filtered data) | Matches *raw* model moments; no HP filter |
| Custom moments | Sectoral std(Y/P/L), cross-sectional **rank correlations**, `corr(Y_i,P_i)`, `std(TB/GDP)`, `autocorr(Q)` | Only product moments in `matched_moments`; rank-corr / ratios **not expressible** |
| Weighting | Diagonal, inverse data-variance (+ ad-hoc multipliers) | Optimal HAC, 1-/2-step/iterated |
| Standard errors / J-test | **None** | Yes (analytical or numerical) |
| Model solution | Dynare.jl native, **hand-rolled Klein fallback** | Dynare core (battle-tested) |
| Optimizer | CMA-ES | CMA-ES + many others |

**Take-away.** Your approach is genuinely *more flexible on the moment side* (the HP-filter correction and the rank-correlation / cross-sectional moments are real value-adds you could not get from Dynare's `matched_moments`). The cost is that you re-implement the model solve and the moment algebra by hand, and you give up Dynare's inference (optimal weighting, standard errors, the overid J-test). For a top-journal submission you will eventually need **standard errors and a J-test**; budget for adding them (sensitivity-based SEs à la Andrews–Gentzkow–Shapiro are cheapest given you already have analytical moments).

---

## 2. Is the model solved correctly? — verification

Two independent, standard checks were run against the *actual compiled model* (`mod/NK_SOE_lev_gap2/model/julia/`, evaluated in a clean Julia 1.11.6):

1. **Steady state is exact.** Evaluating the compiled residual `f` at `dynare_ss.csv` with the `params_jl.mod` parameters gives `max|f| = 3.3e-13`. The steady state is a true fixed point of the model.

2. **Blanchard–Kahn / determinacy holds.** From Dynare's exported decision rule `dynare_g1_1.csv` (= `ghx`), the implied state transition `S = ghx[state_rows, :]` (80×80) has **all eigenvalues inside the unit circle** (max |λ| = 0.9951, zero explosive roots). The solution is unique and stable.

Together these are the conventional evidence that the linear solution is correct. (An independent reconstruction of the equilibrium residual `A·yₜ₋₁ + B·yₜ + C·yₜ₊₁ ≈ 0` was attempted but was inconclusive due to (a) the sparse-index ordering bug in C4 and (b) the extreme scale heterogeneity of the level model — `κᵥ`-bearing Jacobian entries ≈ 1e13 alongside IRF components spanning 1e3…1e-300. This is a limitation of the offline harness, **not** evidence of a solve error; the SS and BK checks above are dispositive.)

**Verdict:** Dynare solves the model correctly. The problems below are in the estimation layer wrapped around it.

---

## 3. Correctness findings (detail)

### C1 — Shock-layout mismatch (Critical)

`run_smm_estimation.jl` loads `mod/nk_iosoe_context.jls`, which is compiled from **`NK_SOE_lev_gap2.mod`**. Its `varexo` order is:

```
1 eps_om   2 eps_i   3 epschi   4 eps_pvstar   5 eps_postar
6..17 epsA_1..epsA_12   18 eps_xi                      (18 shocks)
```

But `smm_model_moments`/`smm_estimation` build `Σe` for the **28-shock `_smm.mod`** Option-A order (`# 1=eps_i, 2=epschi, 3=eps_pvstar, 4:15=epsA_1:12, 16=eps_xi, 17:28=eps_om_1:12`). Applied to the loaded 18-shock model, the unit variances land on the **wrong shocks**:

| `Σe` index set to 1 | code *thinks* it is | in the loaded model it is |
|---|---|---|
| 1 | eps_i | **eps_om** |
| 3 | eps_pvstar | **epschi** |
| 4–15 | epsA_1..12 | **eps_pvstar, eps_postar, epsA_1..10** |
| 16 | eps_xi | **epsA_11** |
| 17–18 | eps_om_1, eps_om_2 | **epsA_12, eps_xi** |

So the analytical moments are computed from a **scrambled shock set** (oil and labour-supply shocks switched on, the monetary shock switched off, sectoral TFP shocks mislabelled). Every moment evaluation is affected.

### C2 — Dead demand-shock channel (Critical)

The lev model's preference shock is a **single aggregate** `eps_om`, entering as `sigma_om*eps_om` (lines 546–547 of `NK_SOE_lev_gap2.mod`), governed by the **scalar** parameter `sigma_om`. The estimator, however:

- estimates `θ[19:30] = sigma_om_1..12` and calls `set_param!(ctx, "sigma_om_$(i)", …)`. Those parameters **do not exist** in the lev model, so `param_idx` returns `nothing` and `set_param!` **silently no-ops** → 12 estimated parameters have **zero effect** on the model.
- never sets the scalar `sigma_om` the model actually uses — and `sigma_om` is **not even written by `params_jl.mod`** (confirmed: `sigma_om` and `rho_om2` are the only two declared parameters absent from `params_jl.mod`). So the demand shock's amplitude is left at its uninitialised default.

Net effect: the entire sectoral-demand dimension of the estimation is non-operative, and 12 of the 35 parameters are unidentified by construction.

### C3 — Klein fallback: stale hardcoded dimensions (High)

`_klein_solve` hardcodes `n_bk=78; n_endo=491; n_fw=56; n_exo=17` (642 Jacobian columns). The actual compiled model is `80 / 570 / 58 / 18` (726 columns: 80+570+58+18). Every slice `A=G[:,1:78]`, `B=G[:,79:569]`, `C=G[:,570:625]`, `D=G[:,626:642]` is therefore wrong. Masked today only because Dynare.jl's native solver runs first; if the fallback ever activates (the ARM-Mac scenario it was written for), it returns garbage or spurious BK failures.

### C4 — Klein fallback: Jacobian sparse-index order mismatch (High)

`_load_jacobian_structure!` parses `(eq, col)` from `dynamic.json` **in file order** (which is **row-major**: eq1's columns 87,92,93,…) and zips them with `g1_v` from the compiled `SparseDynamicG1!` (which is **column-major** — its first assignment is `g1_v[1] = -params[6]`, not the `1.0` that row-major entry 1 would imply). The resulting `sparse(rows, cols, g1_v)` places values at the wrong positions. Same caveat as C3: latent unless the fallback runs.

### C5 — Implicit context selection (Medium)

`run_smm_estimation.jl` picks `nk_iosoe_context.jls` if present else `nk_iosoe_smm_context.jls`. Both exist in `mod/`. Which one loads silently determines whether C1/C2 apply. There is no assertion that the loaded context's shock count/order matches what the moment code assumes.

### C6 — Methodology items (from the earlier `SMM_AUDIT_REPORT.md`)

- HP-filter-of-model: **now implemented** (`hp_filtered_cov_fast`) — good, this was the main earlier concern.
- Monetary shock in `Σe`: intended to be on, but see C1 — under the lev layout index 1 is `eps_om`, so the *monetary* shock is actually off. Fixing C1 resolves this.
- Lag-1 autocovariance spectral formula: looks correct in the current spectral implementation.

---

## 4. Recommended fixes (in priority order)

1. **Decide which model is the estimation target and make it explicit.** Either
   (a) estimate on the **`_smm.mod`** (28-shock, 12 sectoral demand shocks) — then `run_smm_estimation.jl` must load `nk_iosoe_smm_context.jls` and the compiled files under the `_smm` build, and `main_SOE_gap.jl` must (re)build that context; or
   (b) estimate on the **lev model** — then rewrite the parameter vector and `Σe` to the 18-shock layout (drop `sigma_om_1..12`, estimate the scalar `sigma_om`; map TFP via `isigma_tfp_i`).
   Recommendation: pick (a) if the 12 sectoral demand shocks are part of the identification story; otherwise (b) is simpler.

2. **Map `Σe` shocks by NAME, not hardcoded index.** Build `Σe` from `context.models[1].exogenous` names (`findfirst(==("epsA_$i"), exo_names)` etc.) so a model/order change can never silently scramble the shocks. Add this to `validate_smm_setup`.

3. **Make `set_param!` fail loudly on unknown names during setup.** A one-time check that every parameter the estimator writes (`sigma_om_$(i)`, `isigma_tfp_$(i)`, …) actually exists in the loaded context would have caught C2 immediately. (Keep the silent no-op only inside the hot loop, if at all.)

4. **De-hardcode the Klein solver.** Read `n_eq` from `dynamic.json` nrows, `n_exo` from `exogenous_nbr`, `n_bk = length(i_bkwrd_b)`, `n_fw = length(i_fwrd_b)`; and fix the sparse mapping to use the compiled column-major order (or read Dynare's CSC `colptr/rowval`). Even though the fallback is dormant, leaving it silently wrong is a trap.

5. **Assert context consistency at load.** In `run_smm_estimation.jl`, after deserialising, assert `exogenous_nbr == expected` and that the named shocks exist; abort with a clear message otherwise.

6. **(Journal-readiness) Add inference.** Standard errors (sensitivity-based, Andrews–Gentzkow–Shapiro, are cheap given analytical moments) and an overidentification J-test.

---

## 5. Verification harness

The offline checks above were produced with a clean Julia 1.11.6 evaluating the compiled model files directly (no Dynare needed for SS/Jacobian). To reproduce on a machine with the model built:

- **SS check:** evaluate `SparseDynamicResid!` at `dynare_ss.csv` + `params_jl.mod`; expect `max|residual| < 1e-10`.
- **BK check:** load `dynare_g1_1.csv` and `dynare_state_rows.csv`; compute `eigvals(G[state_rows, :])`; expect all `|λ| < 1`.
- **Shock-map check (the important one):** print `context.models[1].exogenous` names in order and compare to the `Σe` index comments in `smm_model_moments`. They currently disagree (C1).

When Dynare is available, the cleanest full check is Dynare's own `model_diagnostics` + `check` (eigenvalue/BK report) and, for the moments, a one-off comparison of your Lyapunov `std(GDP)` against `stoch_simul`'s `oo_.var` at the same θ.
